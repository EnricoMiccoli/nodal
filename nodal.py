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

NODE_TYPES = ["A", "R", "E"]
NODE_TYPES_ANOM = ["E"]

# Compute problem dimensions and find highest degree node
with open(NETLIST, 'r') as infile:
    netlist = csv.reader(infile)
    ncomps = 0
    degrees = {}
    anomalies = 0
    anomnum = {}
    for component in netlist:
        ncomps += 1
        curnodes = component[ACOL:BCOL+1]
        newnodes = [key for key in curnodes if key not in degrees]
        if component[TCOL] in NODE_TYPES_ANOM:
            anomnum[component[NCOL]] = anomalies
            anomalies += 1
        for node in newnodes:
            degrees[node] = 0
        for node in curnodes:
            degrees[node] += 1

    ground = max(degrees.keys(), key=(lambda x: degrees[x]))

    i = 0
    nodenum = {}
    for node in [k for k in degrees.keys() if k != ground]:
        nodenum[node] = i
        i += 1
    assert len(nodenum) == len(degrees) - 1
    

# We are going to solve Gx = A
# Build coefficient matrix 
with open(NETLIST, 'r') as infile:
    n = len(degrees) - 1 + anomalies # number of unknowns
    G = np.zeros(shape=(n, n))
    A = np.zeros(shape=(n, 1))
    netlist = csv.reader(infile)
    currents = []
    for component in netlist:
        assert type(component[TCOL]) == str
        assert component[TCOL] in NODE_TYPES
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
        elif component[TCOL] == "E":
            tension = float(component[VCOL])
            k = anomnum[component[NCOL]]
            i = len(nodenum) + k
            currents.append(component[NCOL])
            A[i] += tension
            if anode != ground:
                j = nodenum[anode]
                G[i,j] = 1
                G[j,i] = -1
            if bnode != ground:
                j = nodenum[bnode]
                G[i,j] = -1
                G[j,i] = 1
        else:
            exit(1)

try:
    e = np.linalg.solve(G, A)
except np.linalg.linalg.LinAlgError:
    print("Model error: matrix is singular")
    print(G)
    exit(1)

nodelabel = dict(reversed(x) for x in nodenum.items())
print("Ground node: {}".format(ground))
i = 0
for potential in e[0:len(nodenum)]:
    print("e({}) = {}".format(nodelabel[i], potential))
    i += 1
i = 0
for current in e[len(nodenum):len(nodenum)+anomalies]:
    print("i({}) = {}".format(currents[i], current))
    i += 1
