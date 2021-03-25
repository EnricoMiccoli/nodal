"""Constants shared between modules."""

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
