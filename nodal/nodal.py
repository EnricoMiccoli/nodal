import numpy as np
import scipy as sp
import scipy.sparse as spsp
import scipy.sparse.linalg as spspla
import csv
import logging

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


def check_input_component(component):
    key = component[NCOL]
    assert type(key) == str

    s = len(component)
    if s < 5:
        logging.error(f"Missign arguments for component {key}")
        raise ValueError
    ctype = component[TCOL]

    if ctype not in NODE_TYPES:
        logging.error(f"Unknown type {ctype} for component {key}")
        raise ValueError

    n = NODE_ARGS_NUMBER[ctype]
    if s != n:
        logging.error(
            f"Wrong number of arguments for component {key}: " f"expected {n}, got {s}"
        )
        raise ValueError


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

    # Iterate over components in the netlist file
    #TODO skip empty lines and comments
    with open(netlist_path, "r") as infile:
        netlist = csv.reader(infile)
        nums = {}
        nums["components"] = 0
        nums["anomalies"] = 0
        nums["be"] = 0  # number of branch equations
        nums["kcl"] = 0  # number of non-ground nodes
        nums["opamps"] = 0
        degrees = {}
        anomnum = {}
        for component in netlist:

            # Check for proper input
            check_input_component(component)

            # Initialize the new component
            key = component[NCOL]
            ctype = component[TCOL]
            component_keys.append(key)
            components[key] = [None] * 8

            newcomp = components[key]
            newcomp[NCOL] = key
            newcomp[TCOL] = ctype

            # Assign value
            try:
                newcomp[VCOL] = float(component[VCOL])
            except ValueError:
                logging.error(
                    "Bad input: expected a number for component value "
                    "of {} got {} instead".format(component[NCOL], component[VCOL])
                )
                raise

            # Assign positive and negative leads
            newcomp[ACOL] = component[ACOL]
            newcomp[BCOL] = component[BCOL]

            # Assign type
            if component[TCOL] in NODE_TYPES_DEP:
                newcomp[CCOL] = component[CCOL]
                newcomp[DCOL] = component[DCOL]
                if component[TCOL] in NODE_TYPES_CC:
                    newcomp[PCOL] = component[PCOL]

            # Update the different component counts
            nums["components"] += 1
            curnodes = component[ACOL : BCOL + 1]
            newnodes = [key for key in curnodes if key not in degrees]
            if component[TCOL] in NODE_TYPES_ANOM:
                anomnum[component[NCOL]] = nums["anomalies"]
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
    [nums, degrees, anomnum, components, component_keys, ground, nodenum] = state
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
            assert i >= 0 and i <= nums["kcl"]
        if bnode != ground:
            j = nodenum[bnode]
            assert j >= 0 and j <= nums["kcl"]

        if component[TCOL] == "R":
            try:
                conductance = 1 / component[VCOL]
            except ZeroDivisionError:
                logging.error("Model error: resistors can't have null resistance")
                raise ValueError
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
                logging.error("Driving component {} not found".format(component[PCOL]))
                raise
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
                logging.error("Unknown component type: {}".format(driver[TCOL]))
                raise ValueError("Unknown component type")

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
                logging.error("Driving component {} not found".format(component[PCOL]))
                raise
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
                logging.error("Unknown component type: {}".format(driver[TCOL]))
                raise ValueError("Unknown component type")

        elif component[TCOL] == "OPAMP":
            raise NotImplementedError

        else:
            logging.error("Unknown component type: {}".format(component[TCOL]))
            raise ValueError("Unknown component type")

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


class Netlist(object):
    def __init__(self, path):
        self.state = read_netlist(path)


class Solution(object):
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


class Circuit(object):
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
