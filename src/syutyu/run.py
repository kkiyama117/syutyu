"""明示された LAMMPS 入力だけを実行する。入力スクリプトは信頼済みコードとして扱う。

名前を挙げた補助ファイルだけをコピーする。LAMMPS の shell コマンドや include、ファイルアクセスを隔離する sandbox ではない。物理単位や系固有の設定は推測しない。
"""
from datetime import datetime, timezone
from hashlib import sha256
import json
import math
import os
from pathlib import Path
import shutil
import signal
import subprocess
import tempfile
import time


_RESERVED = {'screen.txt', 'lammps-help.txt', 'run.json', 'log.lammps'}


def _run_process(command, cwd, output, timeout):
    """On POSIX, stop the local MPI process group on timeout or interruption."""
    process = subprocess.Popen(command, cwd=cwd, env=dict(os.environ, OMP_NUM_THREADS='1'),
                               stdout=output, stderr=subprocess.STDOUT, start_new_session=(os.name == 'posix'))
    try:
        return process.wait(timeout=timeout)
    except (subprocess.TimeoutExpired, KeyboardInterrupt):
        if os.name == 'posix':
            for sig in (signal.SIGTERM, signal.SIGKILL):
                try:
                    os.killpg(process.pid, sig)
                except ProcessLookupError:
                    pass
                if sig == signal.SIGTERM:
                    try:
                        process.wait(timeout=2)
                    except subprocess.TimeoutExpired:
                        pass
        else:
            process.kill()
        process.wait()
        raise


def _input_files(input_file, assets):
    """Validate every source before reading or copying any of them."""
    input_file = Path(input_file).absolute()
    root = input_file.parent.resolve()
    files = {}
    sources = set()
    for requested in (root / input_file.name, *assets):
        requested = Path(requested)
        if '..' in requested.parts:
            raise ValueError(f'Input path traversal is not allowed: {requested}')
        source = requested if requested.is_absolute() else root / requested
        try:
            name = source.relative_to(root)
            resolved = source.resolve(strict=True)
            resolved.relative_to(root)
        except ValueError as error:
            raise ValueError(f'Input must stay within {root}: {requested}') from error
        if not resolved.is_file():
            raise ValueError(f'Input must be a regular file: {requested}')
        if name.name in _RESERVED or name.parts[0] in _RESERVED:
            raise ValueError(f'Reserved output filename: {name}')
        if name in files or resolved in sources:
            raise ValueError(f'Duplicate input or main-input collision: {name}')
        files[name] = resolved
        sources.add(resolved)
    return files


def _resolve_executable(executable):
    resolved = shutil.which(executable)
    if resolved is None:
        raise FileNotFoundError(f'Executable not found: {executable}')
    return str(Path(resolved).resolve())


def run_lammps(input_file, output_root, *, assets=(), ranks=2, timeout=180,
               executable='lmp') -> Path:
    """入力と明示した補助ファイルを新しいディレクトリへ写し、実行してそのパスを返す。

    引数:
        - `input_file` — LAMMPS の入力。親フォルダを基準に補助ファイルを解決する。
        - `output_root` — run フォルダを作る場所。既存でもよい。
        - `assets` — 入力の親フォルダからの相対パス（その配下の絶対パスも可）。サブフォルダは保たれる。
        - `ranks` — MPI プロセス数。1 なら `lmp` を直接、2 以上なら `mpirun -np` で起動する。
        - `timeout` — 1回の実行の上限秒。`lmp -help` は最大30秒。
        - `executable` — LAMMPS の実行ファイル名。PATH から解決する。
    戻り値:
        作った run フォルダの `Path`。`run.json` に入力ハッシュ・コマンド・成否・時間・警告・出力を残す。
    例外:
        `ValueError` — `ranks` や `timeout` が不正、`..` による移動、入力の外を指すパス、`run.json` などの予約名、重複した入力。
        `FileNotFoundError` — 実行ファイルや `mpirun` が見つからない。
        `RuntimeError` — `-help` の失敗、非ゼロ終了、出力の `ERROR` 行、時間切れ、その他の実行失敗。
    注意:
        失敗した run もフォルダと `run.json` を残す。入力と過去の run は書き換えない。POSIX では時間切れや中断時にプロセスグループを止める。
    """
    if not isinstance(ranks, int) or isinstance(ranks, bool) or ranks < 1:
        raise ValueError('ranks must be a positive integer')
    if not math.isfinite(timeout) or timeout <= 0:
        raise ValueError('timeout must be finite and positive')
    files = _input_files(input_file, assets)
    main_name = next(iter(files)).as_posix()
    executable = _resolve_executable(executable)
    command = [executable, '-in', main_name, '-log', 'log.lammps']
    if ranks > 1:
        command = [_resolve_executable('mpirun'), '-np', str(ranks)] + command
    output_root = Path(output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    run_dir = Path(tempfile.mkdtemp(prefix='run-', dir=output_root)).resolve()
    started = time.monotonic()
    record = {
        'created_at': datetime.now(timezone.utc).isoformat(),
        'input_file': main_name,
        'assets': [name.as_posix() for name in list(files)[1:]],
        'input_sha256': {}, 'executable': executable, 'command': command,
        'mpi_ranks': ranks, 'omp_threads': 1, 'timeout_seconds': timeout,
        'status': 'running', 'phase': 'copy', 'returncode': None,
        'warnings': [], 'outputs': [],
    }
    try:
        for name, source in files.items():
            destination = run_dir / name
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, destination)
            record['input_sha256'][name.as_posix()] = sha256(destination.read_bytes()).hexdigest()
        record['phase'] = 'help'
        with (run_dir / 'lammps-help.txt').open('w', encoding='utf-8') as output:
            help_code = _run_process([executable, '-help'], run_dir, output, min(timeout, 30))
        record['help_returncode'] = help_code
        if help_code:
            raise RuntimeError(f'LAMMPS -help exited with code {help_code}')
        record['phase'] = 'simulation'
        with (run_dir / 'screen.txt').open('w', encoding='utf-8') as output:
            record['returncode'] = _run_process(command, run_dir, output, timeout)
        if record['returncode']:
            raise RuntimeError(f'LAMMPS exited with code {record["returncode"]}')
        screen = (run_dir / 'screen.txt').read_text(encoding='utf-8', errors='replace')
        if any(line.lstrip().startswith('ERROR') for line in screen.splitlines()):
            raise RuntimeError('LAMMPS reported an ERROR')
        record['status'] = 'completed'
    except subprocess.TimeoutExpired as error:
        record['status'] = 'timeout'
        record['error'] = str(error)
        raise RuntimeError(f'LAMMPS timed out during {record["phase"]}; inspect {run_dir}') from error
    except KeyboardInterrupt:
        record['status'] = 'interrupted'
        record['error'] = 'KeyboardInterrupt'
        raise
    except Exception as error:
        record['status'] = 'failed'
        record['error'] = str(error)
        raise RuntimeError(f'LAMMPS run failed: {error}; inspect {run_dir}') from error
    finally:
        record['finished_at'] = datetime.now(timezone.utc).isoformat()
        record['elapsed_seconds'] = time.monotonic() - started
        for name in ('screen.txt', 'log.lammps'):
            path = run_dir / name
            if path.is_file():
                for line in path.read_text(encoding='utf-8', errors='replace').splitlines():
                    if 'WARNING' in line and line not in record['warnings']:
                        record['warnings'].append(line)
        input_names = {name.as_posix() for name in files}
        record['outputs'] = sorted({
            path.relative_to(run_dir).as_posix() for path in run_dir.rglob('*')
            if path.is_file() and path.relative_to(run_dir).as_posix() not in input_names
        } | {'run.json'})
        (run_dir / 'run.json').write_text(json.dumps(record, indent=2) + '\n', encoding='utf-8')
    return run_dir
