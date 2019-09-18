import csv
import logging

import numpy as np
import scipy as sp
import scipy.sparse as spsp
import scipy.sparse.linalg as spspla

logging.basicConfig(level=logging.ERROR)

# CSV parsing
NCOL = 0  # component name
TCOL = 1  # type of component
VCOL = 2  # value of component, eg resistance
ACOL = 3  # node connected to first lead, currents enter here
BCOL = 4  # node connected to second lead
# for dependent sources:
CCOL = 5  # first node of controlling variable
DCOL = 6  # second node of controlling variable
PCOL = 7  # name of the driving component

NODE_TYPES_CC = ["CCCS", "CCVS"]
NODE_TYPES_DEP = ["VCVS", "VCCS"] + NODE_TYPES_CC
NODE_TYPES_ANOM = ["E"] + NODE_TYPES_DEP
NODE_TYPES = ["A", "R"] + NODE_TYPES_ANOM + ["OPAMP"]

NODE_ARGS_NUMBER = {
    "OPAMP": 7,
    "R": 5,
    "A": 5,
    "E": 5,
    "VCCS": 7,
    "VCVS": 7,
    "CCCS": 8,
    "CCVS": 8,
}


def find_ground_node(degrees):
    if "g" in degrees:
        ground = "g"
    else:
        ground = max(degrees.keys(), key=(lambda x: degrees[x]))
    logging.debug("ground node-> {}".format(ground))
    return ground


def check_input_component(data):
    s = len(data)
    if s == 0 or data[0][0] == "#":
        return

    key = data[NCOL]
    assert type(key) == str

    if s < 5:
        raise ValueError(f"Missing arguments for component {key}")

    ctype = data[TCOL]

    if ctype not in NODE_TYPES:
        raise ValueError(f"Unknown type {ctype} for component {key}")

    n = NODE_ARGS_NUMBER[ctype]
    if s != n:
        raise ValueError(
            f"Wrong number of arguments for component {key}: expected {n}, got {s}"
        )

    try:
        float(data[VCOL])
    except ValueError:
        raise ValueError(
            "Bad input: expected a number for component value "
            "of {}, got {} instead".format(data[NCOL], data[VCOL])
        )


def build_component(data):
    new_component = [None] * 8

    new_component[NCOL] = data[NCOL]
    new_component[TCOL] = data[TCOL]
    new_component[VCOL] = float(data[VCOL])
    new_component[ACOL] = data[ACOL]
    new_component[BCOL] = data[BCOL]

    if data[TCOL] in NODE_TYPES_DEP:
        new_component[CCOL] = data[CCOL]
        new_component[DCOL] = data[DCOL]
        if data[TCOL] in NODE_TYPES_CC:
            new_component[PCOL] = data[PCOL]

    return new_component


def read_netlist(netlist_path):
    components = {}
    # We will need to iterate over components twice
    # in the same order, so we save keys
    # TODO make this more memory efficient
    component_keys = []

    try:
        infile = open(netlist_path, "r")
    except FileNotFoundError:
        logging.error("File does not exist in specified path")
        raise
    infile.close()

    # Iterate over components in the netlist file
    with open(netlist_path, "r") as infile:
        netlist = csv.reader(infile, skipinitialspace=True)
        nums = {}
        nums["components"] = 0
        nums["anomalies"] = 0
        nums["be"] = 0  # number of branch equations
        nums["kcl"] = 0  # number of non-ground nodes
        nums["opamps"] = 0
        degrees = {}
        anomnum = {}
        for data in netlist:

            # Check for proper input
            check_input_component(data)

            # Skip comments and empty lines
            if data == [] or data[0][0] == "#":
                continue

            # Build the new component
            key = data[NCOL]
            component_keys.append(key)
            newcomp = build_component(data)
            components[key] = newcomp

            # Update the different component counts
            nums["components"] += 1
            curnodes = [data[ACOL], data[BCOL]]
            newnodes = [key for key in curnodes if key not in degrees]
            if data[TCOL] in NODE_TYPES_ANOM:
                anomnum[data[NCOL]] = nums["anomalies"]
                nums["anomalies"] += 1
            for node in newnodes:
                degrees[node] = 0
            for node in curnodes:
                degrees[node] += 1

    # Set ground node
    ground = find_ground_node(degrees)

    # Update node counts
    i = 0
    nodenum = {}
    for node in [k for k in degrees.keys() if k != ground]:
        nodenum[node] = i
        i += 1
    assert len(nodenum) == len(degrees) - 1

    # Update equations count
    logging.debug("nodenum={}".format(nodenum))
    nums["kcl"] = len(nodenum)
    nums["be"] = nums["anomalies"]
    logging.debug("nums={}".format(nums))
    logging.debug("anomnum={}".format(anomnum))
    # From now on nums shall become immutable

    state = [nums, degrees, anomnum, components, component_keys, ground, nodenum]
    # TODO these variables should become attributes of an object
    return state


def build_coefficients(state, sparse):
    [nums, _, anomnum, components, component_keys, ground, nodenum] = state
    n = nums["kcl"] + nums["be"]  # number of unknowns
    if sparse:
        G = spsp.dok_matrix((n, n), dtype=np.float64)
    else:
        G = np.zeros(shape=(n, n))
    A = np.zeros(n)
    currents = []
    for key in component_keys:  # preserve order of iteration
        component = components[key]
        anode = component[ACOL]
        bnode = component[BCOL]
        if anode != ground:
            i = nodenum[anode]
            assert 0 <= i <= nums["kcl"]
        if bnode != ground:
            j = nodenum[bnode]
            assert 0 <= j <= nums["kcl"]

        if component[TCOL] == "R":
            try:
                conductance = 1 / component[VCOL]
            except ZeroDivisionError:
                raise ValueError("Model error: resistors can't have null resistance")
            if anode != ground:
                G[i, i] += conductance
            if bnode != ground:
                G[j, j] += conductance
            if bnode != ground and anode != ground:
                G[i, j] -= conductance
                G[j, i] -= conductance

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
                assert G[i, j] == 0
                G[i, j] = 1
                G[j, i] = -1
            if bnode != ground:
                j = nodenum[bnode]
                assert G[i, j] == 0
                G[i, j] = -1
                G[j, i] = 1

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
                assert G[i, j] == 0
                G[i, j] = 1
                G[j, i] = -1
            if bnode != ground:
                j = nodenum[bnode]
                assert G[i, j] == 0
                G[i, j] = -1
                G[j, i] = 1
            if cnode != ground:
                j = nodenum[cnode]
                G[i, j] += -r
            if dnode != ground:
                j = nodenum[dnode]
                G[i, j] += r

        elif component[TCOL] == "VCCS":
            currents.append(component[NCOL])
            g = component[VCOL]
            k = anomnum[component[NCOL]]
            j = nums["kcl"] + k
            if anode != ground:
                i = nodenum[anode]
                assert G[i, j] == 0
                G[i, j] = -1
            if bnode != ground:
                i = nodenum[bnode]
                assert G[i, j] == 0
                G[i, j] = 1
            # we write the branch equation:
            # i_cccs = g (ec - ed)
            # i_cccs - g ec + g ed = 0
            i = j
            G[i, i] = +1
            cnode = component[CCOL]
            dnode = component[DCOL]
            if cnode != ground:
                j = nodenum[cnode]
                G[i, j] = -g
            if dnode != ground:
                j = nodenum[dnode]
                G[i, j] = +g

        elif component[TCOL] == "CCVS":
            r = component[VCOL]
            k = anomnum[component[NCOL]]
            i = nums["kcl"] + k
            currents.append(component[NCOL])
            cnode = component[CCOL]
            dnode = component[DCOL]
            try:
                driver = components[component[PCOL]]
            except KeyError:
                raise KeyError(f"Driving component {component[PCOL]} not found")
            assert cnode is not None
            assert dnode is not None
            assert driver is not None
            assert (cnode == driver[ACOL] and dnode == driver[BCOL]) or (
                cnode == driver[BCOL] and dnode == driver[ACOL]
            )
            if anode != ground:
                j = nodenum[anode]
                G[i, j] = 1
                G[j, i] = -1
            if bnode != ground:
                j = nodenum[bnode]
                G[i, j] = -1
                G[j, i] = 1

            # we write the branch equation:
            # v_cccv = r * i_driver
            # ea - eb - r * i_driver = 0
            if driver[TCOL] == "R":
                # i_driver = (ec - ed)/R_driver
                if cnode != ground:
                    j = nodenum[cnode]
                    G[i, j] = r / driver[VCOL]
                if dnode != ground:
                    j = nodenum[dnode]
                    G[i, j] = -r / driver[VCOL]
            elif driver[TCOL] in NODE_TYPES_ANOM:
                j = anomnum[driver[NCOL]]
                if driver[ACOL] == component[CCOL]:
                    assert driver[BCOL] == component[DCOL]
                    G[i, j] = -r
                else:
                    assert driver[ACOL] == component[DCOL]
                    assert driver[BCOL] == component[CCOL]
                    G[i, j] = r
            elif driver[TCOL] == "A":
                A[i] = r * driver[VCOL]
            else:
                raise ValueError(f"Unknown component type: {driver[TCOL]}")

        elif component[TCOL] == "CCCS":
            currents.append(component[NCOL])
            g = component[VCOL]
            k = anomnum[component[NCOL]]
            j = nums["kcl"] + k
            if anode != ground:
                i = nodenum[anode]
                assert G[i, j] == 0
                G[i, j] = -1
            if bnode != ground:
                i = nodenum[bnode]
                assert G[i, j] == 0
                G[i, j] = 1
            # we write the branch equation:
            # i_cccs = g * i_driver
            i = j
            assert G[i, i] == 0
            G[i, i] = 1
            try:
                driver = components[component[PCOL]]
            except KeyError:
                raise KeyError(f"Driving component {component[PCOL]} not found")
            # case 1: i_driver is unknown
            if driver[TCOL] == "R":
                cnode = component[CCOL]
                dnode = component[DCOL]
                assert (cnode == driver[ACOL] and dnode == driver[BCOL]) or (
                    cnode == driver[BCOL] and dnode == driver[ACOL]
                )
                assert cnode is not None
                assert dnode is not None
                if cnode != ground:
                    j = nodenum[cnode]
                    assert G[i, j] == 0
                    G[i, j] = +g / driver[VCOL]
                if dnode != ground:
                    j = nodenum[dnode]
                    assert G[i, j] == 0
                    G[i, j] = -g / driver[VCOL]
            elif driver[TCOL] in NODE_TYPES_ANOM:
                j = anomnum[driver[NCOL]]
                if driver[ACOL] == component[CCOL]:
                    assert driver[BCOL] == component[DCOL]
                    G[i, j] = -g
                else:
                    assert driver[ACOL] == component[DCOL]
                    assert driver[BCOL] == component[CCOL]
                    G[i, j] = g
            # case 2: i_driver is known
            elif driver[TCOL] == "A":
                assert A[i] == 0
                A[i] = g * driver[VCOL]
            else:
                raise ValueError(f"Unknown component type: {driver[TCOL]}")

        elif component[TCOL] == "OPAMP":
            raise NotImplementedError

        else:
            raise ValueError(f"Unknown component type: {driver[TCOL]}")

    logging.debug("currents={}".format(currents))
    logging.debug("G=\n{}".format(G))
    logging.debug("A=\n{}".format(A))
    if sparse:
        G = G.tocsr()
    return [G, A, currents]


def solve_system(G, A):
    try:
        e = np.linalg.solve(G, A)
    except np.linalg.linalg.LinAlgError:
        logging.error("Model error: matrix is singular")
        logging.debug(G)
        raise
    return e


def solve_sparse_system(G, A):
    try:
        e = spspla.spsolve(G, A)
    except sp.linalg.LinAlgError:
        logging.error("Model error: matrix is singular")
        logging.debug(G)
        raise
    return e


class Netlist:
    def __init__(self, path):
        self.state = read_netlist(path)


class Solution:
    def __init__(self, result, state, currents):
        self.result = result
        self.nodenum = state[6]
        self.nums = state[0]
        self.currents = currents
        self.ground = state[5]
        self.anomnum = state[2]

    def __str__(self):
        output = "Ground node: {}".format(self.ground)
        names = sorted(self.nodenum)
        for name in names:
            i = self.nodenum[name]
            potential = self.result[i]
            output += "\ne({}) \t= {}".format(name, potential)
        names = sorted(self.anomnum)
        for name in names:
            i = self.anomnum[name]
            current = self.result[self.nums["kcl"] + i]
            output += "\ni({}) \t= {}".format(name, current)
        return output


class Circuit:
    def __init__(self, netlist, sparse=False):
        if not isinstance(netlist, Netlist):
            raise TypeError("Input isn't a netlist")
        self.state = netlist.state
        self.sparse = sparse
        model = build_coefficients(self.state, self.sparse)
        self.G = model[0]
        self.A = model[1]
        self.currents = model[2]

    def solve(self):
        if self.sparse:
            result = solve_sparse_system(self.G, self.A)
        else:
            result = solve_system(self.G, self.A)
        solution = Solution(result, self.state, self.currents)
        return solution
