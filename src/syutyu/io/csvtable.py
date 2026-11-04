"""Named numeric CSV I/O. Trajectories and thermo logs have dedicated readers."""
import numpy as np


def read_csv(path):
    data = np.genfromtxt(path, delimiter=',', names=True, ndmin=1)
    if not data.dtype.names or data.size == 0:
        raise ValueError(f'Empty or unnamed CSV: {path}')
    return {key: data[key] for key in data.dtype.names}


def write_csv(path, table):
    np.savetxt(path, np.column_stack(list(table.values())), delimiter=',',
               header=','.join(table), comments='')
