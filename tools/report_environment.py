"""Write reproducibility metadata without starting any optimization job."""
import argparse
import importlib.metadata
import json
import os
import platform
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('environment.json'))
    args = parser.parse_args()
    versions = {}
    for package in ('gurobipy', 'pyomo', 'numpy', 'matplotlib'):
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = None
    report = {'python': sys.version, 'platform': platform.platform(),
              'logical_cpu_count': os.cpu_count(), 'packages': versions,
              'external_solver_versions': 'Record the executable version headers from the solver logs.'}
    args.output.write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(args.output)


if __name__ == '__main__':
    main()

