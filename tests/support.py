import subprocess
import sys

import MDAnalysis
import numpy as np
from MDAnalysis.coordinates.memory import MemoryReader

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


def memory_universe(steps, positions, boxes, types, velocities=None, angles=90.0):
    """Wrap explicit arrays in a Universe whose ids run 1..N and whose atoms share resid 0."""
    steps = np.asarray(steps)
    n_atoms = len(types)
    universe = MDAnalysis.Universe.empty(n_atoms, n_residues=1, atom_resindex=[0] * n_atoms,
                                         trajectory=True)
    universe.add_TopologyAttr('id', list(range(1, n_atoms + 1)))
    universe.add_TopologyAttr('type', list(types))
    universe.add_TopologyAttr('element', [''] * n_atoms)
    universe.add_TopologyAttr('resid', [0])
    dimensions = np.column_stack([np.asarray(boxes, dtype=float),
                                  np.full((len(steps), 3), angles, dtype=float)])
    universe.load_new(np.asarray(positions, dtype=float), format=MemoryReader,
                      dimensions=dimensions, velocities=velocities,
                      dt=steps[1] - steps[0] if len(steps) > 1 else 1, time_offset=steps[0])
    return universe


def frames_of(universe):
    """Return (steps, positions, boxes, velocities) over every frame, as float64 except steps."""
    steps, positions, boxes, velocities = [], [], [], []
    has_velocities = universe.trajectory.ts.has_velocities
    for timestep in universe.trajectory:
        steps.append(int(timestep.data['step']) if 'step' in timestep.data
                     else int(round(timestep.time)))
        positions.append(np.asarray(timestep.positions, dtype=float))
        boxes.append(np.asarray(timestep.dimensions[:3], dtype=float))
        if has_velocities:
            velocities.append(np.asarray(timestep.velocities, dtype=float))
    universe.trajectory[0]
    return (np.array(steps, dtype=np.int64), np.array(positions), np.array(boxes),
            np.array(velocities) if has_velocities else None)


def run_cli(*args, **kwargs):
    """Run `python -m syutyu` with text capture; kwargs pass through to subprocess.run."""
    return subprocess.run([sys.executable, '-m', 'syutyu', *map(str, args)],
                          capture_output=True, text=True, timeout=30, **kwargs)
