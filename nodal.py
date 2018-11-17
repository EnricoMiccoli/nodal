import numpy as np
import sys
import csv

NETLIST = sys.argv[1]
# CSV parsing
NCOL = 0 # component number
TCOL = 1 # type of component
VCOL = 2 # value of component, eg resistance
ACOL = 3 # node connected to first lead, currents enter here
BCOL = 4 # node connected to second lead

# Compute problem dimensions and find highest degree node
with open(NETLIST, 'r') as infile:
    netlist = csv.reader(infile)
    nodes = set()
    ncomps = 0
    degrees = {}
    for component in netlist:
        ncomps += 1
        curnodes = component[ACOL:BCOL+1]
        newnodes = [key for key in curnodes if key not in degrees]
        for node in newnodes:
            degrees[node] = 0
        for node in curnodes:
            degrees[node] += 1

    f = (len(degrees), ncomps)
    ground = max(degrees.keys(), key=(lambda x: degrees[x]))

    i = 0
    nodenum = {}
    for node in [k for k in degrees.keys() if k != ground]:
        nodenum[node] = i
        i += 1
    assert len(nodenum) == len(degrees) - 1

# Build coefficient matrix 
with open(NETLIST, 'r') as infile:
    G = np.zeros(shape=((f[0]-1,)*2))
    A = np.zeros(shape=((f[0]-1,1)))
    netlist = csv.reader(infile)
    for component in netlist:
        assert type(component[TCOL]) == str
        assert component[TCOL] in ["A", "R"]
        anode = component[ACOL]
        bnode = component[BCOL]
        if anode != ground:
            i = nodenum[anode]
            assert i >= 0 and i <= len(nodenum)
        if bnode != ground:
            j = nodenum[bnode]
            assert j >= 0 and j <= len(nodenum)
        if component[TCOL] == "R":
            conductance = 1 / float(component[VCOL])
            if anode != ground:
                G[i,i] += conductance
            if bnode != ground:
                G[j,j] += conductance
            if bnode != ground and anode != ground:
                G[i,j] -= conductance
                G[j,i] -= conductance
        elif component[TCOL] == "A":
            current = float(component[VCOL])
            if anode != ground:
                A[i] += current
            if bnode != ground:
                A[j] -= current
        else:
            exit(1)

e = np.linalg.solve(G, A)

nodelabel = dict(reversed(x) for x in nodenum.items())
print("Ground node: {}".format(ground))
i = 0
for potential in e:
    print("e({}) = {}".format(nodelabel[i], potential))
    i += 1
