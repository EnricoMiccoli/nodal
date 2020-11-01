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

# Component types
NODE_TYPES_CC = ["CCCS", "CCVS"]
NODE_TYPES_DEP = ["VCVS", "VCCS"] + NODE_TYPES_CC
NODE_TYPES_ANOM = ["E"] + NODE_TYPES_DEP
NODE_TYPES = ["A", "R"] + NODE_TYPES_ANOM + ["OPAMP", "OPMODEL"]

NODE_ARGS_NUMBER = {
    "OPAMP": 7,
    "OPMODEL": 7,
    "R": 5,
    "A": 5,
    "E": 5,
    "VCCS": 7,
    "VCVS": 7,
    "CCCS": 8,
    "CCVS": 8,
}

# OPAMP modeling
OPMODEL_RI = 1e7  # (ohm)
OPMODEL_RO = 10  # (ohm)
OPMODEL_GAIN = 1e5  # (adimensional)


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


def append_opmodel(data, netlist):
    # OPMODEL component specification:
    # [
    #   name,
    #   "OPMODEL",
    #   value of feedback resistor,
    #   output terminal,
    #   ground terminal
    #   non-inverting terminal,
    #   inverting terminal,
    # ]

    name = data[NCOL]

    # Values
    ri = str(OPMODEL_RI)
    ro = str(OPMODEL_RO)
    rf = data[VCOL]
    gain = str(OPMODEL_GAIN)

    # Nodes
    out = data[ACOL]
    ground = data[BCOL]
    pos = data[CCOL]
    neg = data[DCOL]
    phony = f"{name}_internal_node"

    input_resistor = [f"{name}_ri", "R", ri, pos, neg]
    output_resistor = [f"{name}_ro", "R", ro, phony, out]
    feedback_resistor = [f"{name}_rf", "R", rf, neg, out]
    vcvs = [f"{name}_vcvs", "VCVS", gain, phony, ground, pos, neg]

    # rf is set to 0 when there is direct feedback without a resistor
    if rf != "0":
        netlist.append(feedback_resistor)
    else:
        assert neg == out
    netlist.append(input_resistor)
    netlist.append(output_resistor)
    netlist.append(vcvs)


def read_netlist(netlist_path):
    try:
        infile = open(netlist_path, "r")
    except FileNotFoundError:
        logging.error("File does not exist in specified path")
        raise
    infile.close()

    with open(netlist_path, "r") as infile:
        netlist = csv.reader(infile, skipinitialspace=True)
        state = State()
        state.process_netlist(netlist)

    return state


def build_coefficients(state, sparse):
    nums = state.nums
    anomnum = state.anomnum
    components = state.components
    component_keys = state.component_keys
    ground = state.ground
    nodenum = state.nodenum

    n = nums["kcl"] + nums["be"]  # number of unknowns
    if sparse:
        G = spsp.dok_matrix((n, n), dtype=np.float64)
    else:
        G = np.zeros(shape=(n, n))
    A = np.zeros(n)
    currents = []
    for key in component_keys:  # preserve order of iteration
        component = components[key]
        anode = component.anode
        bnode = component.bnode
        if anode != ground:
            i = nodenum[anode]
            assert 0 <= i <= nums["kcl"]
        if bnode != ground:
            j = nodenum[bnode]
            assert 0 <= j <= nums["kcl"]

        if component.type == "R":
            try:
                conductance = 1 / component.value
            except ZeroDivisionError:
                raise ValueError("Model error: resistors can't have null resistance")
            if anode != ground:
                G[i, i] += conductance
            if bnode != ground:
                G[j, j] += conductance
            if bnode != ground and anode != ground:
                G[i, j] -= conductance
                G[j, i] -= conductance

        elif component.type == "A":
            current = component.value
            if anode != ground:
                A[i] += current
            if bnode != ground:
                A[j] -= current

        elif component.type == "E":
            tension = component.value
            k = anomnum[component.name]
            i = nums["kcl"] + k
            currents.append(component.name)
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

        elif component.type == "VCVS":
            currents.append(component.name)
            r = component.value
            k = anomnum[component.name]
            i = nums["kcl"] + k
            cnode = component.pos_control
            dnode = component.neg_control
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

        elif component.type == "VCCS":
            currents.append(component.name)
            g = component.value
            k = anomnum[component.name]
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
            cnode = component.pos_control
            dnode = component.neg_control
            if cnode != ground:
                j = nodenum[cnode]
                G[i, j] = -g
            if dnode != ground:
                j = nodenum[dnode]
                G[i, j] = +g

        elif component.type == "CCVS":
            r = component.value
            k = anomnum[component.name]
            i = nums["kcl"] + k
            currents.append(component.name)
            cnode = component.pos_control
            dnode = component.neg_control
            try:
                driver = components[component.driver]
            except KeyError:
                raise KeyError(f"Driving component {component.driver} not found")
            assert cnode is not None
            assert dnode is not None
            assert driver is not None
            assert (cnode == driver.anode and dnode == driver.bnode) or (
                cnode == driver.bnode and dnode == driver.anode
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
            if driver.type == "R":
                # i_driver = (ec - ed)/R_driver
                if cnode != ground:
                    j = nodenum[cnode]
                    G[i, j] = r / driver.value
                if dnode != ground:
                    j = nodenum[dnode]
                    G[i, j] = -r / driver.value
            elif driver.type in NODE_TYPES_ANOM:
                j = anomnum[driver.name]
                if driver.anode == component.pos_control:
                    assert driver.bnode == component.neg_control
                    G[i, j] = -r
                else:
                    assert driver.anode == component.neg_control
                    assert driver.bnode == component.pos_control
                    G[i, j] = r
            elif driver.type == "A":
                A[i] = r * driver.value
            else:
                raise ValueError(f"Unknown component type: {driver.type}")

        elif component.type == "CCCS":
            currents.append(component.name)
            g = component.value
            k = anomnum[component.name]
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
                driver = components[component.driver]
            except KeyError:
                raise KeyError(f"Driving component {component.driver} not found")
            # case 1: i_driver is unknown
            if driver.type == "R":
                cnode = component.pos_control
                dnode = component.neg_control
                assert (cnode == driver.anode and dnode == driver.bnode) or (
                    cnode == driver.bnode and dnode == driver.anode
                )
                assert cnode is not None
                assert dnode is not None
                if cnode != ground:
                    j = nodenum[cnode]
                    assert G[i, j] == 0
                    G[i, j] = +g / driver.value
                if dnode != ground:
                    j = nodenum[dnode]
                    assert G[i, j] == 0
                    G[i, j] = -g / driver.value
            elif driver.type in NODE_TYPES_ANOM:
                j = anomnum[driver.name]
                if driver.anode == component.pos_control:
                    assert driver.bnode == component.neg_control
                    G[i, j] = -g
                else:
                    assert driver.anode == component.neg_control
                    assert driver.bnode == component.pos_control
                    G[i, j] = g
            # case 2: i_driver is known
            elif driver.type == "A":
                assert A[i] == 0
                A[i] = g * driver.value
            else:
                raise ValueError(f"Unknown component type: {driver.type}")

        elif component.type == "OPAMP":
            raise NotImplementedError

        else:
            raise ValueError(f"Unknown component type: {driver.type}")

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


class Component:
    def __init__(self, data):
        check_input_component(data)

        self.name = data[NCOL]
        self.type = data[TCOL]
        self.value = float(data[VCOL])
        self.anode = data[ACOL]
        self.bnode = data[BCOL]

        if data[TCOL] in NODE_TYPES_DEP:
            self.pos_control = data[CCOL]
            self.neg_control = data[DCOL]
            if data[TCOL] in NODE_TYPES_CC:
                self.driver = data[PCOL]
            else:
                self.driver = None
        else:
            self.pos_control = None
            self.neg_control = None


class Netlist:
    def __init__(self, path):
        self.state = read_netlist(path)


class State:
    def __init__(self):
        self.nums = {"components": 0, "anomalies": 0, "be": 0, "kcl": 0, "opamps": 0}
        self.degrees = {}
        self.anomnum = {}
        self.components = {}
        self.component_keys = []
        self.ground = None
        self.nodenum = {}

    def process_netlist(self, netlist):
        # Iterate over components in the netlist file
        for data in netlist:
            # Skip comments and empty lines
            if data == [] or data[0][0] == "#":
                continue

            # Build the new component
            try:
                newcomp = Component(data)
            except ValueError:
                raise
            key = data[NCOL]
            # We will need to iterate over components twice
            # in the same order, so we save keys
            # TODO make this more memory efficient
            self.component_keys.append(key)
            self.components[key] = newcomp

            # Update the different component counts
            self.nums["components"] += 1
            curnodes = [data[ACOL], data[BCOL]]
            newnodes = [key for key in curnodes if key not in self.degrees]
            if data[TCOL] in NODE_TYPES_ANOM:
                self.anomnum[data[NCOL]] = self.nums["anomalies"]
                self.nums["anomalies"] += 1
            for node in newnodes:
                self.degrees[node] = 0
            for node in curnodes:
                self.degrees[node] += 1

        # Set ground node
        self.ground = find_ground_node(self.degrees)

        # Update node counts
        i = 0
        self.nodenum = {}
        for node in [k for k in self.degrees.keys() if k != self.ground]:
            self.nodenum[node] = i
            i += 1
        assert len(self.nodenum) == len(self.degrees) - 1

        # Update equations count
        logging.debug("nodenum={}".format(self.nodenum))
        self.nums["kcl"] = len(self.nodenum)
        self.nums["be"] = self.nums["anomalies"]
        logging.debug("nums={}".format(self.nums))
        logging.debug("anomnum={}".format(self.anomnum))
        # From now on nums shall become immutable


class Solution:
    def __init__(self, result, state, currents):
        self.result = result
        self.nodenum = state.nodenum
        self.nums = state.nums
        self.currents = currents
        self.ground = state.ground
        self.anomnum = state.anomnum

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
