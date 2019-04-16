import nodal as n
import sys

if __name__ == '__main__':
    if len(sys.argv) >= 2:
        netlist_path = sys.argv[1]
    else:
        print("Missing argument: netlist file")
        exit(1)
    if len(sys.argv) == 3:
        sparse = (sys.argv[2] == "sparse")
    else:
        sparse = False

    netlist = n.Netlist(netlist_path)
    circuit = n.Circuit(netlist, sparse=sparse)
    solution = circuit.solve()
    print(solution)
