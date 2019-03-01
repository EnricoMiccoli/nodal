import numpy as np
import sys
import csv
import logging
logging.basicConfig(level=logging.DEBUG)

NETLIST = sys.argv[1]
# CSV parsing
NCOL = 0 # component name
TCOL = 1 # type of component
VCOL = 2 # value of component, eg resistance
ACOL = 3 # node connected to first lead, currents enter here
BCOL = 4 # node connected to second lead
# for dependent sources:
CCOL = 5 # first node of controlling variable
DCOL = 6 # second node of controlling variable
PCOL = 7 # name of the driving component

NODE_TYPES_CC = ["CCCS", "CCVS"]
NODE_TYPES_DEP = ["VCVS", "VCCS"] + NODE_TYPES_CC
NODE_TYPES_ANOM = ["E"] + NODE_TYPES_DEP
NODE_TYPES = ["A", "R"] + NODE_TYPES_ANOM

components = {}
# We will need to iterate over components twice
# in the same order, so we save keys
# TODO make this more memory efficient
component_keys = []

# Compute problem dimensions and find highest degree node
with open(NETLIST, 'r') as infile:
    netlist = csv.reader(infile)
    nums = {}
    nums["components"] = 0
    nums["anomalies"] = 0
    nums["be"] = 0  # number of branch equations
    nums["kcl"] = 0 # number of non-ground nodes
    degrees = {}
    anomnum = {}
    for component in netlist:
        key = component[NCOL]
        assert type(key) == str
        component_keys.append(key)
        components[key] = [None] * 8

        newcomp = components[key]
        newcomp[NCOL] = key
        assert component[TCOL] in NODE_TYPES
        newcomp[TCOL] = component[TCOL]

        try:
            newcomp[VCOL] = float(component[VCOL])
        except ValueError:
            print((
                "Bad input: expected a number"
                "for component value of {},"
                "got {} instead.").format(
                    component[NCOL],
                    component[VCOL])
                )
            exit(1)

        newcomp[ACOL] = component[ACOL]
        newcomp[BCOL] = component[BCOL]
        if component[TCOL] in NODE_TYPES_DEP:
            newcomp[CCOL] = component[CCOL]
            newcomp[DCOL] = component[DCOL]
            if component[TCOL] in NODE_TYPES_CC:
                newcomp[PCOL] = component[PCOL]
            else:
                assert len(component) == 7
        else:
            assert len(component) == 5

        nums["components"] += 1
        curnodes = component[ACOL:BCOL+1]
        newnodes = [key for key in curnodes
                if key not in degrees]
        if component[TCOL] in NODE_TYPES_ANOM:
            anomnum[component[NCOL]] = nums["anomalies"]
            nums["anomalies"] += 1
        for node in newnodes:
            degrees[node] = 0
        for node in curnodes:
            degrees[node] += 1

ground = max(degrees.keys(), key=(lambda x: degrees[x]))
#ground = sys.argv[2]
logging.debug("ground node-> {}".format(ground))

i = 0
nodenum = {}
for node in [k for k in degrees.keys() if k != ground]:
    nodenum[node] = i
    i += 1
assert len(nodenum) == len(degrees) - 1
logging.debug("nodenum={}".format(nodenum))
nums["kcl"] = len(nodenum)
nums["be"] = nums["anomalies"]
logging.debug("nums={}".format(nums))
logging.debug("anomnum={}".format(anomnum))
# From now on nums shall become immutable

# We are going to solve Ge = A
# Build coefficient matrix 
n = nums["kcl"] + nums["be"] # number of unknowns
G = np.zeros(shape=(n, n))
A = np.zeros(shape=(n, 1))
currents = []
for key in component_keys: # preserve order of iteration
    component = components[key]
    anode = component[ACOL]
    bnode = component[BCOL]
    if anode != ground:
        i = nodenum[anode]
        assert i >= 0 and i <= nums["kcl"]
    if bnode != ground:
        j = nodenum[bnode]
        assert j >= 0 and j <= nums["kcl"]
    if component[TCOL] == "R":
        try:
            conductance = 1 / component[VCOL]
        except ZeroDivisionError:
            print("Model error: resistors can't have null resistance")
            exit(1)
        if anode != ground:
            G[i,i] += conductance
        if bnode != ground:
            G[j,j] += conductance
        if bnode != ground and anode != ground:
            G[i,j] -= conductance
            G[j,i] -= conductance
    elif component[TCOL] == "A":
        current = component[VCOL]
        if anode != ground:
            A[i] += current
        if bnode != ground:
            A[j] -= current
    elif component[TCOL] == "E":
        tension = component[VCOL]
        k = anomnum[component[NCOL]]
        i = nums["kcl"] + k
        currents.append(component[NCOL])
        A[i] += tension
        if anode != ground:
            j = nodenum[anode]
            assert G[i,j] == 0
            G[i,j] = 1
            G[j,i] = -1
        if bnode != ground:
            j = nodenum[bnode]
            assert G[i,j] == 0
            G[i,j] = -1
            G[j,i] = 1
    elif component[TCOL] == "VCVS":
        currents.append(component[NCOL])
        r = component[VCOL]
        k = anomnum[component[NCOL]]
        i = nums["kcl"] + k
        cnode = component[CCOL]
        dnode = component[DCOL]
        # we need to write this BE:
        # ea - eb = r (ec - ed)
        # ea - eb - r ec + r ed = 0
        if anode != ground:
            j = nodenum[anode]
            assert G[i,j] == 0
            G[i,j] = 1
            G[j,i] = -1
        if bnode != ground:
            j = nodenum[bnode]
            assert G[i,j] == 0
            G[i,j] = -1
            G[j,i] = 1
        if cnode != ground:
            j = nodenum[cnode]
            G[i,j] += -r
        if dnode != ground:
            j = nodenum[dnode]
            G[i,j] += r
    elif component[TCOL] == "CCCS":
        currents.append(component[NCOL])
        g = component[VCOL]
        k = anomnum[component[NCOL]]
        j = nums["kcl"] + k
        if anode != ground:
            i = nodenum[anode]
            assert G[i,j] == 0
            G[i,j] = -g
        if bnode != ground:
            i = nodenum[bnode]
            assert G[i,j] == 0
            G[i,j] = g
        # we write the branch equation:
        # i_cccs = i_driver
        i = j
        assert G[i,i] == 0
        G[i,i] = 1
        driver = components[component[PCOL]]
        # case 1: i_driver is unknown 
        if driver[TCOL] == 'R':
            cnode = component[CCOL]
            dnode = component[DCOL]
            assert (cnode == driver[ACOL] and dnode == driver[BCOL]) or (cnode == driver[BCOL] and dnode == driver[ACOL])
            assert cnode != None
            assert dnode != None
            if cnode != ground:
                j = nodenum[cnode]
                assert G[i,j] == 0
                G[i,j] = +1 / driver[VCOL]
            if dnode != ground:
                j = nodenum[dnode]
                assert G[i,j] == 0
                G[i,j] = -1 / driver[VCOL]
        elif driver[TCOL] in NODE_TYPES_ANOM:
            j = anomnum[driver[NCOL]]
            if driver[ACOL] == component[CCOL]:
                assert driver[BCOL] == component[DCOL]
                G[i,j] = -1
            else:
                assert driver[ACOL] == component[DCOL]
                assert driver[BCOL] == component[CCOL]
                G[i,j] = 1
        # case 2: i_driver is known 
        elif driver[TCOL] == 'A':
            assert A[i] == 0
            A[i] = driver[VCOL]
        else:
            exit(1)
    else:
        exit(1)

logging.debug("currents={}".format(currents))
logging.debug("G=\n{}".format(G))
logging.debug("A=\n{}".format(A))
try:
    e = np.linalg.solve(G, A)
except np.linalg.linalg.LinAlgError:
    print("Model error: matrix is singular")
    logging.error(G)
    exit(1)

nodelabel = dict(reversed(x) for x in nodenum.items())
print("Ground node: {}".format(ground))
i = 0
for potential in e[0:nums["kcl"]]:
    print("e({}) = {}".format(nodelabel[i], potential))
    i += 1
i = 0
for current in e[nums["kcl"]:nums["kcl"]+nums["be"]]:
    print("i({}) = {}".format(currents[i], current))
    i += 1
