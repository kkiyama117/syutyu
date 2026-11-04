from hashlib import sha256
import json
import os
from pathlib import Path
import sys
import tempfile
import threading
import time
import unittest
from unittest.mock import patch

from syutyu import run_lammps

HELP_OK = 'if [ "$1" = "-help" ]; then echo "fake LAMMPS help"; exit 0; fi\n'


def alive(pid):
    try:
        return Path(f'/proc/{pid}/stat').read_text().split()[2] != 'Z'
    except (FileNotFoundError, ProcessLookupError):
        return False


@unittest.skipUnless(os.name == 'posix', 'POSIX fake executables')
class RunnerTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.inputs = self.root / 'inputs'
        self.inputs.mkdir()
        self.main = self.inputs / 'salt.in'
        self.main.write_text('units real\nread_data data/ions.data\ninclude settings/pair.in\n')
        self.output = self.root / 'runs'
        self.lmp = self.fake('lmp', HELP_OK)

    def asset(self, name, text='fixture'):
        path = self.inputs / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
        return path

    def fake(self, name, script):
        path = self.root / 'bin' / name
        path.parent.mkdir(exist_ok=True)
        path.write_text('#!/bin/sh\n' + script)
        path.chmod(0o755)
        return path

    def launch(self, input_file=None, **options):
        options = {'ranks': 1, 'executable': str(self.lmp), **options}
        return run_lammps(input_file or self.main, self.output, **options)

    def test_arbitrary_input_assets_and_relative_executable(self):
        first = self.asset('data/ions.data', 'ionic data')
        second = self.asset('settings/pair.in', 'pair settings')
        self.asset('.private-fixture', 'unlisted secret fixture; never copied')
        exe = self.root / 'fake-lmp'
        exe.write_text(
            '#!/bin/sh\n'
            'if [ "$1" = "-help" ]; then echo "fake LAMMPS help"; exit 0; fi\n'
            'test "$1" = "-in" && test "$2" = "salt.in" || exit 3\n'
            'test "$3" = "-log" && test "$4" = "log.lammps" || exit 4\n'
            'test -f data/ions.data && test -f settings/pair.in || exit 5\n'
            'echo "WARNING: synthetic warning"\n'
            'printf "Step Temp\\n0 300\\n1 301\\nLoop time of 0.1\\n" > log.lammps\n'
            'mkdir results; echo "arbitrary output" > results/ions.dump\n'
        )
        exe.chmod(0o755)
        old = Path.cwd()
        try:
            os.chdir(self.root)
            run = run_lammps(self.main, self.output, assets=['data/ions.data', second],
                             ranks=1, executable='./fake-lmp')
            next_run = run_lammps(self.main, self.output, assets=[first, second],
                                  ranks=1, executable='./fake-lmp')
        finally:
            os.chdir(old)
        self.assertNotEqual(run, next_run)
        record = json.loads((run / 'run.json').read_text())
        self.assertEqual(record['status'], 'completed')
        self.assertEqual(record['executable'], str(exe))
        self.assertEqual(record['command'], [str(exe), '-in', 'salt.in', '-log', 'log.lammps'])
        self.assertEqual(record['input_file'], 'salt.in')
        self.assertEqual(record['assets'], ['data/ions.data', 'settings/pair.in'])
        self.assertEqual(record['input_sha256'], {
            'salt.in': sha256(self.main.read_bytes()).hexdigest(),
            'data/ions.data': sha256(b'ionic data').hexdigest(),
            'settings/pair.in': sha256(b'pair settings').hexdigest(),
        })
        self.assertEqual((run / 'data/ions.data').read_text(), 'ionic data')
        self.assertEqual(record['warnings'], ['WARNING: synthetic warning'])
        self.assertEqual(record['outputs'], ['lammps-help.txt', 'log.lammps',
                                             'results/ions.dump', 'run.json', 'screen.txt'])
        self.assertGreaterEqual(record['elapsed_seconds'], 0)
        self.assertIn('finished_at', record)
        self.assertFalse((run / '.private-fixture').exists())
        self.assertFalse((run / 'thermo.csv').exists())
        self.assertNotIn('conditions', record)

    def test_missing_main_and_assets_refused_before_launch(self):
        with self.assertRaises(FileNotFoundError):
            self.launch(self.inputs / 'absent.in')
        with self.assertRaises(FileNotFoundError):
            self.launch(assets=['absent.data'])
        self.assertFalse(self.output.exists())

    def test_path_traversal_and_outside_assets_refused(self):
        outside = self.root / 'outside.data'
        outside.write_text('not an input')
        for name in ('../outside.data', outside, 'nested/../salt.in'):
            with self.subTest(name=name), self.assertRaises(ValueError):
                self.launch(assets=[name])
        self.assertFalse(self.output.exists())

    def test_duplicate_main_and_assets_refused(self):
        asset = self.asset('data/ions.data')
        for assets in ([self.main], ['salt.in'], [asset, 'data/ions.data']):
            with self.subTest(assets=assets), self.assertRaises(ValueError):
                self.launch(assets=assets)

    def test_reserved_input_and_asset_names_refused(self):
        for name in ('screen.txt', 'lammps-help.txt', 'run.json', 'log.lammps'):
            source = self.asset(name)
            with self.subTest(name=name):
                with self.assertRaisesRegex(ValueError, 'Reserved'):
                    self.launch(assets=[name])
                with self.assertRaisesRegex(ValueError, 'Reserved'):
                    self.launch(source)
                self.asset('nested/' + name)
                with self.assertRaisesRegex(ValueError, 'Reserved'):
                    self.launch(assets=['nested/' + name])

    def test_reserved_output_cannot_be_an_asset_directory(self):
        self.asset('screen.txt/data')
        with self.assertRaisesRegex(ValueError, 'Reserved'):
            self.launch(assets=['screen.txt/data'])
        self.assertFalse(self.output.exists())

    def test_symlink_escape_and_duplicate_alias_refused(self):
        outside = self.root / 'outside.data'
        outside.write_text('outside fixture')
        (self.inputs / 'escape.data').symlink_to(outside)
        (self.inputs / 'external').symlink_to(self.root, target_is_directory=True)
        for assets in (['escape.data'], ['external/outside.data']):
            with self.subTest(assets=assets), self.assertRaises(ValueError):
                self.launch(assets=assets)
        with self.assertRaises(ValueError):
            self.launch(self.inputs / 'escape.data')
        (self.inputs / 'alias.in').symlink_to(self.main)
        with self.assertRaises(ValueError):
            self.launch(assets=['alias.in'])
        self.assertFalse(self.output.exists())

    def test_directory_asset_refused(self):
        (self.inputs / 'data').mkdir()
        with self.assertRaisesRegex(ValueError, 'regular file'):
            self.launch(assets=['data'])

    def test_missing_executable_creates_no_run(self):
        with patch.dict(os.environ, PATH=str(self.root / 'empty-path')):
            with self.assertRaisesRegex(FileNotFoundError, 'Executable not found'):
                run_lammps(self.main, self.output, ranks=1)
        self.assertFalse(self.output.exists())

    def test_mpi_launcher_is_resolved(self):
        mpirun = self.fake('mpirun', 'test "$1" = "-np" || exit 9\nshift 2\nexec "$@"\n')
        with patch.dict(os.environ, PATH=f'{mpirun.parent}{os.pathsep}{os.environ["PATH"]}'):
            run = self.launch(ranks=3)
        record = json.loads((run / 'run.json').read_text())
        self.assertEqual(record['status'], 'completed')
        self.assertEqual(record['command'][:4],
                         [str(mpirun.resolve()), '-np', '3', str(self.lmp.resolve())])

    def test_failed_and_timed_out_runs_keep_metadata(self):
        cases = [('exit-7', 'exit 7\n', {}, 'failed', 'code 7'),
                 ('screen-error', 'echo "ERROR: synthetic failure"\n', {}, 'failed', 'reported an ERROR'),
                 ('sleeper', 'exec sleep 30\n', {'timeout': .5}, 'timeout', 'timed out')]
        for name, script, options, status, message in cases:
            with self.subTest(status=status, message=message):
                before = set(self.output.glob('run-*'))
                with self.assertRaisesRegex(RuntimeError, message) as raised:
                    self.launch(executable=str(self.fake(name, HELP_OK + script)), **options)
                run, = set(self.output.glob('run-*')) - before
                record = json.loads((run / 'run.json').read_text())
                self.assertIn(str(run), str(raised.exception))
                self.assertEqual(record['status'], status)
                self.assertIn('error', record)
                self.assertTrue((run / 'salt.in').is_file())
                self.assertTrue((run / 'screen.txt').is_file())
                self.assertIn('run.json', record['outputs'])

    def test_help_failure_retains_diagnostics(self):
        broken = self.fake('broken-help', 'echo "broken help"\nexit 2\n')
        with self.assertRaisesRegex(RuntimeError, 'help exited with code 2'):
            self.launch(executable=str(broken))
        run, = self.output.iterdir()
        record = json.loads((run / 'run.json').read_text())
        self.assertEqual(record['phase'], 'help')
        self.assertEqual(record['status'], 'failed')
        self.assertIn('broken help', (run / 'lammps-help.txt').read_text())

    def test_invalid_execution_parameters(self):
        for ranks in (0, -1, 1.5, True):
            with self.subTest(ranks=ranks), self.assertRaises(ValueError):
                run_lammps(self.main, self.output, ranks=ranks)
        for timeout in (0, -1, float('nan'), float('inf')):
            with self.subTest(timeout=timeout), self.assertRaises(ValueError):
                run_lammps(self.main, self.output, timeout=timeout)

    @unittest.skipUnless(Path('/proc').is_dir(), 'Linux process state check')
    def test_timeout_stops_worker_ignoring_term(self):
        worker = self.root / 'worker.py'
        worker.write_text(
            'import subprocess, sys, time\n'
            'from pathlib import Path\n'
            "child = subprocess.Popen([sys.executable, '-c', 'import signal, time; "
            "signal.signal(signal.SIGTERM, signal.SIG_IGN); time.sleep(30)'])\n"
            "Path('child.pid').write_text(str(child.pid))\n"
            'time.sleep(30)\n')
        slow = self.fake('slow-lmp', HELP_OK + f'exec "{sys.executable}" "{worker}"\n')
        with self.assertRaisesRegex(RuntimeError, 'timed out'):
            self.launch(executable=str(slow), timeout=1)
        run, = self.output.iterdir()
        self.assertEqual(json.loads((run / 'run.json').read_text())['status'], 'timeout')
        pid = int((run / 'child.pid').read_text())
        # Bounded kernel scheduling allowance; no live worker may remain.
        deadline = time.monotonic() + 1
        while alive(pid) and time.monotonic() < deadline:
            threading.Event().wait(.01)
        self.assertFalse(alive(pid), 'TERM-ignoring worker survived the timeout')


if __name__ == '__main__':
    unittest.main()
