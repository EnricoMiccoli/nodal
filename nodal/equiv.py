"""Script meant for command line usage, exported as `nodal-resistance`.

Also provides the equivalent_resistance() function.
"""

import argparse
from copy import deepcopy

import nodal as n


parser = argparse.ArgumentParser(
    description="Calculate equivalent resistance using nodal analysis"
    "\n"
    "Label nodes as '1' and 'g' to mark where to connect to the network."
)
parser.add_argument(
    "netlist_path", metavar="FILE", help="csv file describing the resistive network"
)
parser.add_argument("-s", "--sparse", action="store_true", help="use a sparse matrix")


def check_resistive(netlist):
    """Check if all components in netlist are resistors"""

    for component in netlist.components.values():
        if component.type != "R":
            return False
    return True


def equivalent_resistance(netlist, a, b, sparse=False):
    """Calculate equivalent resistance as seen through nodes a and b.

    Works by adding a 1 ampere current generator between a and b, then
    solving the circuit. The resistance returned then is
        R = [e(a) - e(b)] / 1

    Raises
        * ValueError: the original netlist has a non-resistor component
        * KeyError: either node a or b not found in netlist
    """

    if not check_resistive(netlist):
        raise ValueError("Network is not resistive")

    for node in (a, b):
        if node not in netlist.nodenum and node != netlist.ground:
            raise KeyError(f"Node `{node}` not found in netlist")

    resistive_network = deepcopy(netlist)
    resistive_network.process_component(["a1", "A", "1", a, b])
    circuit = n.Circuit(resistive_network, sparse=sparse)
    solution = circuit.solve()

    e = [0, 0]
    for i, node in enumerate((a, b)):
        if node != "g":
            j = solution.nodenum[node]
            e[i] = solution.result[j]

    return e[0] - e[1]


def main():
    args = parser.parse_args()
    A = "1"
    B = "g"

    try:
        netlist = n.Netlist(args.netlist_path)
    except FileNotFoundError:
        exit(1)

    try:
        r = equivalent_resistance(netlist, A, B, sparse=args.sparse)
    except ValueError:
        print("Invalid netlist\n")
        print("Resistors are the only component allowed in the circuit")
        exit(1)
    except KeyError:
        print("Invalid netlist\n")
        print("Node '1' or 'g' not found.")
        exit(1)

    print(f"R = {r}")


if __name__ == "__main__":
    main()
