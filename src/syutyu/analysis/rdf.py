"""種別とサイトを分けた RDF。力場に依存しない。
"""
import numpy as np

from ..io.dump import frame_indices
from ..io.store import Result
from ..kernels import minimum_image


def radial_distribution(trajectory, pairs, *, rmax, bins=200, start_step=None, stop_step=None,
                        stride=1, intermolecular=False, normalization):
    """中心と隣接の組ごとに、フレーム平均の RDF と積算配位数を返す。

    引数:
        - `trajectory` — `Trajectory`。
        - `pairs` — `(中心, 隣接)` の選択式の並び。例 `('type:1', 'type:2')`。
        - `rmax` — g(r) の上限。Å。対象全フレームの最短辺の半分以下。
        - `bins` — ビン数。既定 200。
        - `start_step`、`stop_step`、`stride` — 使うフレームの範囲と間引き。
        - `intermolecular` — 同じ分子（mol が同じで 0 でない）の原子対を除く。
        - `normalization` — `finite`（分母 `N_A*N_B-N_overlap`）か `population`（`N_A*N_B`）。必須。
    戻り値:
        `Result`。`tables['rdf']` に `r_center_A`、`r_upper_A`、`g_<i>`、`CN_<i>`、`ordered_counts_<i>`。
    例外:
        `ValueError` — `normalization` が2つ以外、`rmax` が範囲外、`pairs` が空、`intermolecular` なのに mol 列がない。
    注意:
        各フレームの体積で規格化し、CN は中心原子1個あたりの直接カウント。可変体積では平均 g に平均密度を掛け直さない。LAMMPS の `special_bonds` の除外は再現しない。
    """
    if normalization not in ('population', 'finite'):
        raise ValueError('normalization must be population or finite')
    frames = frame_indices(trajectory, start_step=start_step, stop_step=stop_step, stride=stride)
    if (not np.isfinite(rmax) or rmax <= 0 or bins < 1
            or rmax > np.min(trajectory.boxes[frames]) / 2):
        raise ValueError('RDF requires 0 < rmax <= half the shortest box in every selected frame')
    if not pairs:
        raise ValueError('Specify at least one center/neighbor pair')
    if intermolecular and 'mol' not in trajectory.labels:
        raise ValueError('Intermolecular RDF requires a mol column')
    edges = np.linspace(0, rmax, bins + 1)
    shell = 4 * np.pi / 3 * np.diff(edges**3)
    table = {'r_center_A': (edges[:-1] + edges[1:]) / 2, 'r_upper_A': edges[1:]}
    pair_info = []
    for pair_index, (center, neighbor) in enumerate(pairs):
        a = np.flatnonzero(trajectory.select(center))
        b = np.flatnonzero(trajectory.select(neighbor))
        denominator = len(a) * len(b)
        if normalization == 'finite':
            denominator -= len(np.intersect1d(trajectory.ids[a], trajectory.ids[b]))
        if denominator == 0:
            raise ValueError('Finite RDF normalization requires at least one non-self atom pair')
        counts = np.zeros(bins, dtype=np.int64)
        g_sum = np.zeros(bins)
        same_selection = np.array_equal(a, b)
        for frame in frames:
            xyz, box = trajectory.positions[frame], trajectory.boxes[frame]
            hist = np.zeros(bins, dtype=np.int64)
            for offset in range(0, len(a), 64):
                batch = a[offset:offset + 64]
                delta = minimum_image(xyz[batch, None, :] - xyz[b][None, :, :], box)
                keep = (trajectory.ids[batch, None] < trajectory.ids[b][None, :] if same_selection
                        else trajectory.ids[batch, None] != trajectory.ids[b][None, :])
                if intermolecular:
                    mol = trajectory.labels['mol']
                    same_molecule = (mol[batch, None] == mol[b][None, :]) & (mol[batch, None] != '0')
                    keep &= ~same_molecule
                distances = np.sqrt(np.sum(delta**2, axis=-1))[keep]
                hist += np.histogram(distances[distances < rmax], edges)[0]
            if same_selection:
                hist *= 2  # unique i<j pairs, then account for both center atoms
            counts += hist
            g_sum += hist * np.prod(box) / (denominator * shell)
        g_key, cn_key = f'g_{pair_index}', f'CN_{pair_index}'
        table[g_key] = g_sum / len(frames)
        table[cn_key] = np.cumsum(counts) / (len(frames) * len(a))
        table[f'ordered_counts_{pair_index}'] = counts
        pair_info.append({'center': center, 'neighbor': neighbor, 'g': g_key, 'cn': cn_key,
                          'centers': len(a), 'neighbors': len(b)})
    return Result({'rdf': table}, {
        'kind': 'rdf', 'title': 'Partial radial distribution', 'atoms': len(trajectory.ids),
        'source_dump': trajectory.source_path,
        'selected_frames': len(frames), 'selected_steps': [int(trajectory.steps[frames[0]]), int(trajectory.steps[frames[-1]])],
        'frame_stride': stride, 'bin_width_A': float(edges[1]), 'rdf_pairs': pair_info,
        'intermolecular': intermolecular, 'variable_box': not trajectory.fixed_box,
        'method': 'custom', 'rdf_normalization': normalization,
        'normalization': ('Per-frame (N_A*N_B-N_overlap)/V; CN per center atom' if normalization == 'finite'
                          else 'Per-frame N_A*N_B/V; excludes self pairs; CN per center atom'),
        'limitations': ['Equilibrium not confirmed', 'All pairs unless intermolecular requested; no special_bonds exclusions',
                        'With changing volume, CN is averaged directly, not reconstructed using one average density'],
    })
