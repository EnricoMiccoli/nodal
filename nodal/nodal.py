"""Central module of the nodal package.

Provides, among others, the classes meant for external usage:
    * Netlist: reads .csv files
    * Circuit: provides solve() method to compute electrical variables
    * Solution: printable object storing the computation results

Example use case:
    from nodal import Circuit, Netlist
    my_netlist = Netlist("path/to/netlist.csv")
    my_circuit = Circuit(my_netlist, sparse=True)
    my_solution = my_circuit.solve()
    print(my_solution)
"""

import csv
import logging

import numpy as np
import scipy as sp
import scipy.sparse as spsp
import scipy.sparse.linalg as spspla

import nodal.models as models

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


def build_opmodel(data):
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

    result = [input_resistor, output_resistor, vcvs]

    # rf is set to 0 when there is direct feedback without a resistor
    if rf != "0":
        result.append(feedback_resistor)
    else:
        assert neg == out

    return result


class Component:
    def __init__(self, data):
        self.check_input(data)

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

    def check_input(self, data):
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


class Netlist:
    def __init__(self, path):
        self.nums = {"components": 0, "anomalies": 0, "be": 0, "kcl": 0, "opamps": 0}
        self.degrees = {}
        self.anomnum = {}
        self.components = {}
        self.component_keys = []
        self.ground = None
        self.nodenum = {}
        self.opmodel_equivalents = []
        self.read_netlist(path)

    def process_component(self, data):
        # Skip comments and empty lines
        if data == [] or data[0][0] == "#":
            return

        # If the current component is an OPMODEL,
        # replace it with an equivalent circuit.
        if data[TCOL] == "OPMODEL":
            eq = build_opmodel(data)
            self.opmodel_equivalents.extend(eq)
            return

        # Otherwise, just build the component
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

    def read_netlist(self, path):
        try:
            infile = open(path, "r")
        except FileNotFoundError:
            logging.error(f"File '{path}' not found.")
            raise
        infile.close()

        with open(path, "r") as infile:
            reader = csv.reader(infile, skipinitialspace=True)

            # Iterate over components in the netlist file
            for data in reader:
                self.process_component(data)

            for data in self.opmodel_equivalents:
                self.process_component(data)

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


class Circuit:
    def __init__(self, netlist, sparse=False):
        if not isinstance(netlist, Netlist):
            raise TypeError("Input isn't a netlist")
        self.netlist = netlist
        self.sparse = sparse
        self.G, self.A, self.currents = self.build_model()

    def solve(self):
        try:
            if self.sparse:
                e = spspla.spsolve(self.G, self.A)
            else:
                e = np.linalg.solve(self.G, self.A)
        except (np.linalg.linalg.LinAlgError, sp.linalg.LinAlgError):
            logging.error("Model error: matrix is singular")
            logging.debug(self.G)
            raise
        return Solution(e, self.netlist, self.currents)

    def build_model(self):
        # Setup local variables
        nums = self.netlist.nums
        anomnum = self.netlist.anomnum
        components = self.netlist.components
        component_keys = self.netlist.component_keys
        ground = self.netlist.ground
        nodenum = self.netlist.nodenum

        # Initialize matrices
        n = nums["kcl"] + nums["be"]  # number of unknowns
        if self.sparse:
            G = spsp.dok_matrix((n, n), dtype=np.float64)
        else:
            G = np.zeros(shape=(n, n))
        A = np.zeros(n)
        currents = []

        # Iterate over components
        for key in component_keys:  # preserve order of iteration
            component = components[key]
            if component.anode != ground:
                i = nodenum[component.anode]
                assert 0 <= i <= nums["kcl"]
            else:
                i = None
            if component.bnode != ground:
                j = nodenum[component.bnode]
                assert 0 <= j <= nums["kcl"]
            else:
                j = None

            args = (component, i, j, ground, G, A, currents, anomnum, nums, nodenum)
            if component.type == "R":
                models.write_R(component, i, j, ground, G)
            elif component.type == "A":
                models.write_A(component, i, j, ground, A)
            elif component.type == "E":
                models.write_E(*args)
            elif component.type == "VCCS":
                models.write_VCVS(*args)
            elif component.type == "VCVS":
                models.write_VCVS(*args)
            elif component.type == "CCVS":
                models.write_CCVS(*args)
            elif component.type == "CCCS":
                models.write_CCCS(*args, components)
            elif component.type == "OPAMP":
                raise NotImplementedError
            else:
                raise ValueError(f"Unknown component type: {driver.type}")

        # Log and return
        logging.debug("currents={}".format(currents))
        logging.debug("G=\n{}".format(G))
        logging.debug("A=\n{}".format(A))
        if self.sparse:
            G = G.tocsr()
        return [G, A, currents]


class Solution:
    def __init__(self, result, netlist, currents):
        self.result = result
        self.nodenum = netlist.nodenum
        self.nums = netlist.nums
        self.currents = currents
        self.ground = netlist.ground
        self.anomnum = netlist.anomnum

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
