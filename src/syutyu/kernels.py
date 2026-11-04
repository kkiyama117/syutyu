"""解析と教材用の実装が共有する、NumPy だけの数値カーネル。
"""
import numpy as np


def minimum_image(delta, box):
    """相対ベクトルを各方向 [-L/2, L/2] に戻す。

    引数:
        - `delta` — 相対ベクトル。Å。末尾の軸が x, y, z。
        - `box` — 直方体セルの辺の長さ。Å。
    戻り値:
        最小像へ移した相対ベクトル。Å。
    注意:
        直方体セルのみ。三斜晶の tilt は扱わない。
    """
    return delta - box * np.rint(delta / box)


def pair_rdf(positions, box, edges, *, normalization):
    """全フレームの i<j を数え、g(r) と積算配位数を返す。

    引数:
        - `positions` — フレーム×原子×3 の座標。Å。
        - `box` — 直方体セルの辺の長さ。Å。
        - `edges` — ビンの境界。0 から始まり、最短辺の半分以下で終わること。
        - `normalization` — `population`（密度 `N/V`）か `finite`（`(N-1)/V`、Colab の有限N補正）。必須。
    戻り値:
        `(g, cn, counts)`。`counts` はビンごとの原子対数、`cn` はビン上端までの平均隣接原子数（自分自身を含まない）。
    例外:
        `ValueError` — `normalization` が2つ以外、フレームがないか原子が2個未満、ビンの範囲が不正。
    注意:
        有限Nの理想気体で `population` は `1-1/N`、`finite` は 1 になる。規格化を既定値で継承させないため、キーワード必須にしている。
    """
    if normalization not in ('population', 'finite'):
        raise ValueError('normalization must be population or finite')
    frames, n, _ = positions.shape
    if frames == 0 or n < 2:
        raise ValueError('Need frames and at least two atoms')
    if edges[0] != 0 or np.any(np.diff(edges) <= 0) or edges[-1] > np.min(box) / 2:
        raise ValueError('Require bins from zero to at most half the shortest box length')
    i, j = np.triu_indices(n, k=1)  # 自己対と二重計数を除く
    counts = np.zeros(len(edges) - 1, dtype=np.int64)
    for xyz in positions:
        delta = minimum_image(xyz[i] - xyz[j], box)
        distances = np.sqrt(np.sum(delta**2, axis=1))
        counts += np.histogram(distances[distances < edges[-1]], bins=edges)[0]
    shell = 4 * np.pi / 3 * np.diff(edges**3)  # 球殻の厳密な体積
    rho = (n - 1 if normalization == 'finite' else n) / np.prod(box)
    neighbors = 2 * counts / (frames * n)  # i<jの各対を両原子から見直す
    g = neighbors / (rho * shell)
    cn = np.cumsum(neighbors)
    return g, cn, counts


def unwrap(positions, box):
    """隣接フレームの最小像変位を積算して非折り返し座標を作る。

    引数:
        - `positions` — フレーム×原子×3 の折り返し座標。Å。
        - `box` — 直方体セルの辺の長さ。Å。固定セルを前提とする。
    戻り値:
        先頭フレームを基準に積算した座標。Å。
    注意:
        各成分の真の変位が L/2 未満であることが前提。image flag のない dump では、この前提そのものを確かめられない。
    """
    increments = minimum_image(np.diff(positions, axis=0), box)
    return np.concatenate([positions[:1], positions[:1] + np.cumsum(increments, axis=0)])


def msd_multi_origin(unwrapped):
    """lag ごとに可能な全開始時刻と原子を平均した MSD を返す。

    引数:
        - `unwrapped` — フレーム×原子×3 の非折り返し座標。Å。
    戻り値:
        `(msd, samples)`。`msd[lag]` は Å²、`samples[lag]` は使った (開始時刻, 原子) の組の数。
    注意:
        重なった標本は独立でない。lag 0 は 0 として返す。
    """
    frames, n, _ = unwrapped.shape
    msd = np.zeros(frames)
    for lag in range(1, frames):
        displacement = unwrapped[lag:] - unwrapped[:-lag]
        msd[lag] = np.mean(np.sum(displacement**2, axis=-1))
    return msd, (frames - np.arange(frames)) * n
