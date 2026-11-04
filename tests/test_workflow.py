import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from syutyu import rdf, save_result
from support import run_cli, tiny_dump


class WorkflowTests(unittest.TestCase):
    def test_headless_cli_from_another_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            dump = tiny_dump(tmp)
            results = save_result(rdf(dump, method='custom', pairs=[('all', 'all')], rmax=4, bins=8),
                                  tmp / 'rdf')
            out = tmp / 'rdf.png'
            env = dict(os.environ, DISPLAY='', WAYLAND_DISPLAY='', MPLBACKEND='Agg')
            args = ('plot', results, '--output', out)
            completed = run_cli(*args, cwd=tmp, env=env)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertTrue(out.read_bytes().startswith(b'\x89PNG'))
            again = run_cli(*args, cwd=tmp, env=env)
            self.assertEqual(again.returncode, 2)
            self.assertIn('overwrite', again.stderr)

    def test_import_does_not_import_pyplot(self):
        check = subprocess.run([sys.executable, '-c',
                                "import syutyu,sys; assert 'matplotlib.pyplot' not in sys.modules"],
                               capture_output=True, text=True, timeout=10)
        self.assertEqual(check.returncode, 0, check.stderr)


if __name__ == '__main__':
    unittest.main()
