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

    if os.path.isfile(args.netlist_path):
        netlist = n.Netlist(args.netlist_path)
    else:
        print("File not found: {}".format(args.netlist_path))
        exit(1)
    circuit = n.Circuit(netlist, sparse=args.sparse)
    solution = circuit.solve()
    print(solution)


if __name__ == "__main__":
    main()
