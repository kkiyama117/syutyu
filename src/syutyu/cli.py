"""Thin CLI over current APIs: no backward-compatible command aliases."""
import argparse
import json
from pathlib import Path
import subprocess


def make_parser():
    parser = argparse.ArgumentParser(prog='syutyu', description='MD tools: generic LAMMPS runs, partial RDF, atomic MSD, plots and VMD export.')
    commands = parser.add_subparsers(dest='command', required=True)
    run = commands.add_parser('run', help='Run an explicit LAMMPS input and listed assets')
    run.add_argument('input', type=Path)
    run.add_argument('--asset', action='append', default=[], help='File relative to input parent; repeat for data/force fields/includes')
    run.add_argument('--output-root', type=Path, required=True)
    run.add_argument('--ranks', type=int, default=2)
    run.add_argument('--timeout', type=float, default=180.)
    run.add_argument('--executable', default='lmp')
    for name in ('rdf', 'msd', 'vacf'):
        cmd = commands.add_parser(name, help=f'{name.upper()}: select custom dump analysis or native LAMMPS output')
        cmd.add_argument('trajectory', type=Path, help='Raw dump (custom) or native compute output (lammps)')
        cmd.add_argument('--method', choices=('custom', 'lammps'), default='custom')
        cmd.add_argument('--output', type=Path, required=True)
        cmd.add_argument('--start-step', type=int)
        cmd.add_argument('--stride', type=int, default=1, help='Use every Nth frame (MSD unwraps before thinning)')
        cmd.add_argument('--overwrite', action='store_true')
        if name == 'rdf':
            cmd.add_argument('--pair', nargs=2, action='append', metavar=('CENTER', 'NEIGHBOR'))
            cmd.add_argument('--rmax', type=float, help='Custom only: Å, at most half the shortest box length')
            cmd.add_argument('--bins', type=int, help='Custom only; default 200')
            cmd.add_argument('--normalization', choices=('finite', 'population'), default='finite')
            cmd.add_argument('--intermolecular', action='store_true', help='Custom only: exclude same nonzero mol ID')
        else:
            cmd.add_argument('--select', action='append', help='all, type:3, element:Na or mol:12; repeatable')
            cmd.add_argument('--dt-ps', type=float, required=True, help='ps per LAMMPS timestep, NOT dump interval')
            cmd.add_argument('--origins', choices=('multiple', 'single'), help='Custom default multiple; LAMMPS is single only')
            cmd.add_argument('--origin-step', type=int, help='Required for native files: when the compute reference was created')
            if name == 'msd':
                cmd.add_argument('--laboratory', action='store_true', help='Custom NPT with explicit unwrapped coordinates; includes box motion')
            else:
                cmd.add_argument('--velocity-unit', choices=('A/fs', 'A/ps'), required=True)
                cmd.add_argument('--max-lag-ps', type=float)
    native_msd = commands.add_parser('native-msd', help='Import LAMMPS compute msd tables (useful for wrapped NPT dumps)')
    native_msd.add_argument('--series', nargs=2, action='append', required=True, metavar=('NAME', 'FILE'))
    native_msd.add_argument('--dt-ps', type=float, required=True)
    native_msd.add_argument('--origin-step', type=int, default=0)
    native_msd.add_argument('--column', type=int, default=-1, help='Zero-based total-MSD column; -1 is last')
    native_msd.add_argument('--output', type=Path, required=True)
    native_msd.add_argument('--overwrite', action='store_true')
    thermo = commands.add_parser('thermo', help='Extract one thermo block; never merge timestep resets')
    thermo.add_argument('log', type=Path)
    thermo.add_argument('--block', type=int, default=-1, help='0-based; -1 is last block')
    thermo.add_argument('--dt-ps', type=float, required=True)
    thermo.add_argument('--output', type=Path, required=True)
    thermo.add_argument('--overwrite', action='store_true')
    plot = commands.add_parser('plot', help='Save a headless figure')
    plot.add_argument('results', type=Path)
    plot.add_argument('--view', choices=('rdf', 'cn', 'msd', 'vacf', 'vacf-normalized', 'temperature'))
    plot.add_argument('--output', type=Path, required=True)
    plot.add_argument('--no-compare', action='store_true')
    plot.add_argument('--overwrite', action='store_true')
    export = commands.add_parser('export-vmd', help='Create portable trajectory + VMD loader ZIP')
    export.add_argument('trajectory', type=Path)
    export.add_argument('--output', type=Path, required=True)
    export.add_argument('--stride', type=int, default=1)
    return parser


def main(argv=None):
    parser = make_parser()
    args = parser.parse_args(argv)
    try:
        from .io.store import Result, load_result, save_result
        if args.command == 'run':
            from .run import run_lammps
            directory = run_lammps(args.input, args.output_root, assets=args.asset, ranks=args.ranks,
                                   timeout=args.timeout, executable=args.executable)
            print(json.dumps({'run_dir': str(directory)}))
        elif args.command in ('rdf', 'msd', 'vacf'):
            from .analysis.dispatch import msd, rdf, vacf
            custom = args.method == 'custom'
            options = {'start_step': args.start_step, 'stride': args.stride} if custom else {}
            if not custom and (args.stride != 1 or (args.command != 'rdf' and args.start_step is not None)):
                raise ValueError('Native kinetic files cannot change the compute origin using --start-step/--stride')
            if args.command == 'rdf':
                options['pairs'] = args.pair
                if custom:
                    if not args.pair or args.rmax is None:
                        raise ValueError('Custom RDF requires --pair and --rmax')
                    options.update(rmax=args.rmax, bins=200 if args.bins is None else args.bins,
                                   intermolecular=args.intermolecular)
                else:
                    options['start_step'] = args.start_step
                    if args.rmax is not None or args.bins is not None or args.intermolecular:
                        raise ValueError('--rmax/--bins/--intermolecular apply only to custom RDF')
                result = rdf(args.trajectory, method=args.method, normalization=args.normalization, **options)
            else:
                if args.select is not None or custom:
                    options['selections'] = args.select or ['all']
                if args.command == 'msd':
                    if args.laboratory or custom:
                        options['laboratory'] = args.laboratory
                    result = msd(args.trajectory, method=args.method, dt_ps=args.dt_ps,
                                 origin_step=args.origin_step, origin_mode=args.origins, **options)
                else:
                    result = vacf(args.trajectory, method=args.method, dt_ps=args.dt_ps,
                                  velocity_unit=args.velocity_unit, origin_step=args.origin_step,
                                  origin_mode=args.origins, max_lag_ps=args.max_lag_ps, **options)
            save_result(result, args.output, overwrite=args.overwrite)
            print(json.dumps({'results': str(args.output.resolve()), 'summary': result.summary}, indent=2))
        elif args.command == 'native-msd':
            from .io.native import read_msd
            if len(dict(args.series)) != len(args.series):
                raise ValueError('MSD series names must be unique')
            result = read_msd(dict(args.series), dt_ps=args.dt_ps,
                              origin_step=args.origin_step, column=args.column)
            save_result(result, args.output, overwrite=args.overwrite)
            print(json.dumps({'results': str(args.output.resolve())}))
        elif args.command == 'thermo':
            from .io.thermo import read_thermo
            table = read_thermo(args.log, block=args.block, dt_ps=args.dt_ps)
            result = Result({'thermo': table},
                            {'kind': 'thermo', 'title': 'LAMMPS thermodynamic log',
                             'block': args.block, 'dt_ps': args.dt_ps})
            save_result(result, args.output, overwrite=args.overwrite)
            print(json.dumps({'results': str(args.output.resolve()), 'records': table.height}))
        elif args.command == 'plot':
            from .report.figures import save_figure
            save_figure(load_result(args.results), args.output, args.view, not args.no_compare,
                        overwrite=args.overwrite)
            print(json.dumps({'plot': str(args.output.resolve())}))
        else:
            from .report.bundle import export_vmd
            record = export_vmd(args.trajectory, args.output, stride=args.stride)
            print(json.dumps({'bundle': str(args.output.resolve()), **record}, indent=2))
    except (OSError, ValueError, IndexError, RuntimeError, subprocess.SubprocessError) as error:
        parser.exit(2, f'syutyu: {error}\n')
    return 0
