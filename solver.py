import nodal as n
import sys

if __name__ == '__main__':
    if len(sys.argv) == 2:
        netlist_path = sys.argv[1]
    else:
        print("Missing argument: netlist file")
        exit(1)

    netlist = n.Netlist(netlist_path)
    circuit = n.Circuit(netlist)
    solution = circuit.solve()
    print(solution)
