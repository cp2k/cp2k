"""Prepare a common CP2K relaxation; execute only with --run or --submit."""

# SPDX-License-Identifier: GPL-2.0-or-later

import argparse
import json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "structure", help="ASE-readable structure file (fully periodic)"
    )
    parser.add_argument("--profile", required=True, help="Existing AiiDA profile")
    parser.add_argument(
        "--code", required=True, help="Configured CP2K code label or UUID"
    )
    parser.add_argument(
        "--electronic-type", required=True, choices=("metal", "insulator")
    )
    parser.add_argument(
        "--protocol", default="moderate", choices=("fast", "moderate", "precise")
    )
    parser.add_argument(
        "--relax-type",
        default="positions",
        choices=("none", "positions", "positions_cell"),
    )
    parser.add_argument("--mpi", action="store_true")
    parser.add_argument("--ranks", type=int, default=1)
    parser.add_argument(
        "--walltime", type=int, default=3600, help="Walltime in seconds"
    )
    execution = parser.add_mutually_exclusive_group()
    execution.add_argument(
        "--run", action="store_true", help="Run synchronously (no daemon required)"
    )
    execution.add_argument(
        "--submit", action="store_true", help="Submit to the configured AiiDA daemon"
    )
    args = parser.parse_args()
    if args.ranks < 1 or args.walltime < 1 or (args.ranks != 1 and not args.mpi):
        parser.error("Use positive resources and --mpi when requesting multiple ranks")

    from aiida import engine, load_profile, orm
    from ase.io import read
    from cp2k.aiida import build_relax_builder

    atoms = read(args.structure)
    if atoms.get_initial_charges().any() or atoms.get_initial_magnetic_moments().any():
        parser.error(
            "This example is neutral/nonmagnetic; configure charge and spin explicitly for other systems"
        )
    if atoms.info.get("charge", 0) or atoms.info.get("spin", 0):
        parser.error("This example does not infer charge or spin from file metadata")
    load_profile(args.profile)
    builder = build_relax_builder(
        orm.StructureData(ase=atoms),
        code=orm.load_code(args.code),
        electronic_type=args.electronic_type,
        protocol=args.protocol,
        relax_type=args.relax_type,
        options={
            "resources": {"num_machines": 1, "num_mpiprocs_per_machine": args.ranks},
            "withmpi": args.mpi,
            "max_wallclock_seconds": args.walltime,
            "environment_variables": {"OMP_NUM_THREADS": "1"},
        },
    )
    builder.metadata.label = "CP2K common relaxation"
    if args.submit:
        print(f"Submitted workflow UUID: {engine.submit(builder).uuid}")
    elif args.run:
        results, node = engine.run_get_node(builder)
        print(f"Workflow UUID: {node.uuid}; exit status: {node.exit_status}")
        if not node.is_finished_ok:
            raise SystemExit(node.exit_message or "CP2K workflow failed")
        print(f"Total energy: {results['total_energy'].value:.12f} eV")
    else:
        print(json.dumps(builder.cp2k.parameters.get_dict(), indent=2))
        print("Prepared only; no job submitted. Use --run or --submit to execute.")


if __name__ == "__main__":
    main()
