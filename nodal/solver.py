import argparse
import os

import nodal as n

parser = argparse.ArgumentParser(
    description="Solve electrical circuits using nodal analysis"
)
parser.add_argument(
    "netlist_path", metavar="FILE", help="csv file describing the netlist"
)
parser.add_argument("-s", "--sparse", action="store_true", help="use a sparse matrix")


def main():
    args = parser.parse_args()
    path = args.netlist_path

    if os.path.isfile(path):
        netlist = n.Netlist(path)
    else:
        print(f"File not found: {path}")
        exit(1)

    circuit = n.Circuit(netlist, sparse=args.sparse)
    solution = circuit.solve()
    print(solution)


if __name__ == "__main__":
    main()
