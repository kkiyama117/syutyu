import subprocess
import sys

CUBE = ((0, 10),) * 3


def write_dump(path, columns, frames, steps=None, boxes=None):
    """Write a LAMMPS text dump; frames are lists of atom rows, boxes per-frame (lo, hi) triples."""
    blocks = []
    for index, rows in enumerate(frames):
        step = index * 10 if steps is None else steps[index]
        bounds = CUBE if boxes is None else boxes[index]
        blocks.append(f'ITEM: TIMESTEP\n{step}\nITEM: NUMBER OF ATOMS\n{len(rows)}\n'
                      'ITEM: BOX BOUNDS pp pp pp\n' + ''.join(f'{lo} {hi}\n' for lo, hi in bounds)
                      + f'ITEM: ATOMS {columns}\n' + ''.join(row + '\n' for row in rows))
    path.write_text(''.join(blocks))
    return path


def tiny_dump(directory):
    """Three two-atom Ar frames in a fixed 10 Å cube."""
    return write_dump(directory / 'tiny.lammpstrj', 'id element x y z',
                      [['1 Ar 1 1 1', '2 Ar 3 1 1']] * 3, steps=[0, 1000, 2000])


def run_cli(*args, **kwargs):
    """Run `python -m syutyu` with text capture; kwargs pass through to subprocess.run."""
    return subprocess.run([sys.executable, '-m', 'syutyu', *map(str, args)],
                          capture_output=True, text=True, timeout=30, **kwargs)
