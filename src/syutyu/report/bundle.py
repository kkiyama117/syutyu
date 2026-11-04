"""VMD 用の持ち出し ZIP を作る。network・認証・GUI を使わない。
"""
from hashlib import sha256
import json
import os
from pathlib import Path
import tempfile
from zipfile import ZipFile, ZIP_DEFLATED

VMD_SCRIPT = '''# Run in VMD: source /path/to/extracted/load-vmd.tcl
set bundle_dir [file dirname [file normalize [info script]]]
mol new [file join $bundle_dir trajectory.lammpstrj] type lammpstrj waitfor all
mol modstyle 0 top VDW 0.4 12.0
mol modcolor 0 top Name
display resetview
animate goto start
'''


def iter_dump_blocks(path):
    """1フレーム分の生の dump を、ファイル全体を読み込まずに順に返す。

    引数:
        - `path` — custom dump。
    戻り値:
        フレームごとのバイト列。原子列は問わずそのまま保つ。
    例外:
        `ValueError` — 見出しが9行でそろわない、原子数が正でない、セル境界が2列か3列でない、原子行が足りない。
    注意:
        持ち出しのための検査用で、数値の解析には使わない。一般三斜晶やベクトル形式のセルは扱わない。
    """
    with Path(path).open('rb') as stream:
        while first := stream.readline():
            header = [first] + [stream.readline() for _ in range(8)]
            if (header[0].strip() != b'ITEM: TIMESTEP'
                    or header[2].strip() != b'ITEM: NUMBER OF ATOMS'
                    or not header[4].startswith(b'ITEM: BOX BOUNDS')
                    or not header[8].startswith(b'ITEM: ATOMS ')):
                raise ValueError('Incomplete or unsupported LAMMPS dump header')
            int(header[1])
            count = int(header[3])
            if count <= 0:
                raise ValueError('Atom count must be positive')
            for line in header[5:8]:
                if len(line.split()) not in (2, 3):
                    raise ValueError('Unsupported box bounds')
                list(map(float, line.split()))
            columns = len(header[8].split()) - 2
            rows = [stream.readline() for _ in range(count)]
            if any(len(row.split()) != columns for row in rows):
                raise ValueError('Incomplete atom records')
            yield b''.join(header + rows)


def export_vmd(trajectory, output, *, stride=1):
    """指定フレームと VMD スクリプト・manifest・説明を1つの ZIP にまとめる。

    引数:
        - `trajectory` — custom dump。
        - `output` — 作る ZIP。既存のパスは拒否する。
        - `stride` — 何フレームおきに保存するか（1以上）。1 なら全バイトをそのまま残す。
    戻り値:
        manifest の辞書。元ファイルと保存フレームの SHA-256、枚数を含む。
    例外:
        `FileExistsError` — `output` が既にある。
        `FileNotFoundError` — 軌跡がない。
        `ValueError` — `stride` が正の整数でない、保存フレームが0。
    注意:
        間引きは表示用で、動力学の解析に使うと元の軌跡と違う結果になる。完成した ZIP だけを公開し、途中のファイルは残さない。
    """
    trajectory, output = Path(trajectory), Path(output)
    if not isinstance(stride, int) or isinstance(stride, bool) or stride < 1:
        raise ValueError('stride must be a positive integer')
    if output.exists():
        raise FileExistsError(f'Refusing to overwrite {output}')
    if not trajectory.is_file():
        raise FileNotFoundError(trajectory)
    output.parent.mkdir(parents=True, exist_ok=True)
    # Publish only a complete ZIP; the temporary file is on the same filesystem.
    with tempfile.TemporaryDirectory(prefix='.vmd-bundle-', dir=output.parent) as temporary:
        archive_path = Path(temporary) / 'bundle.zip'
        source_hash, saved_hash = sha256(), sha256()
        source_frames = saved_frames = 0
        with ZipFile(archive_path, 'w', compression=ZIP_DEFLATED) as archive:
            with archive.open('trajectory.lammpstrj', 'w', force_zip64=True) as destination:
                for index, block in enumerate(iter_dump_blocks(trajectory)):
                    source_hash.update(block)
                    source_frames += 1
                    if index % stride == 0:
                        destination.write(block)
                        saved_hash.update(block)
                        saved_frames += 1
            if not saved_frames:
                raise ValueError('Empty trajectory')
            manifest = {
                'schema_version': 1, 'source_name': trajectory.name,
                'source_sha256': source_hash.hexdigest(), 'trajectory_sha256': saved_hash.hexdigest(),
                'source_frames': source_frames, 'exported_frames': saved_frames, 'stride': stride,
                'purpose': 'Visualization only; scientific interpretation requires original run conditions',
            }
            archive.writestr('manifest.json', json.dumps(manifest, indent=2) + '\n')
            archive.writestr('load-vmd.tcl', VMD_SCRIPT)
            archive.writestr('README.txt',
                'VMD用の持ち出しセット\n\n'
                '1. このZIPをノートPCで展開してください（ZIP内から直接実行しない）。\n'
                '2. VMDのTk Consoleで、展開したスクリプトを読み込みます。\n'
                '   source {/path/to/extracted/load-vmd.tcl}\n'
                '   Windows例: source {C:/Users/you/Downloads/ar-vmd/load-vmd.tcl}\n'
                '   波括弧は空白を含むパスにも対応します。実際の展開先へ置き換えてください。\n'
                '   コマンドからなら: vmd -e load-vmd.tcl\n\n'
                f'全{source_frames}フレーム中{saved_frames}枚を保存（stride={stride}）。\n'
                '軌跡の間引きは表示用です。MSD等の解析には元の軌跡を使ってください。\n'
                '計算条件・密度・平衡化は元の計算記録で確認してください。\n'
                '配布資料・パスワード・他の作業ファイルはこのZIPへ含めません。\n')
        # Hard-link publication is atomic and refuses an existing destination.
        os.link(archive_path, output)
    return manifest
