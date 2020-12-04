"""Core implementation of the nodal analysis procedure.

Provides write_COMPONENT() for all component types. These functions are
used exclusively by the Circuit.build_model() method when writing the
    G e = A
linear system.

The matrix G and the vector A are often referenced and modified by
these functions.
"""


def write_R(c, i, j, ground, G):
    try:
        conductance = 1 / c.value
    except ZeroDivisionError:
        raise ValueError("Model error: resistors can't have null resistance")
    if c.anode != ground:
        G[i, i] += conductance
    if c.bnode != ground:
        G[j, j] += conductance
    if c.bnode != ground and c.anode != ground:
        G[i, j] -= conductance
        G[j, i] -= conductance


def write_A(c, i, j, ground, A):
    current = c.value
    if c.anode != ground:
        A[i] += current
    if c.bnode != ground:
        A[j] -= current


def write_E(c, i, j, ground, G, A, currents, anomnum, nums, nodenum):
    tension = c.value
    k = anomnum[c.name]
    i = nums["kcl"] + k
    currents.append(c.name)
    A[i] += tension
    if c.anode != ground:
        j = nodenum[c.anode]
        assert G[i, j] == 0
        G[i, j] = 1
        G[j, i] = -1
    if c.bnode != ground:
        j = nodenum[c.bnode]
        assert G[i, j] == 0
        G[i, j] = -1
        G[j, i] = 1


def write_VCVS(c, i, j, ground, G, A, currents, anomnum, nums, nodenum):
    currents.append(c.name)
    r = c.value
    k = anomnum[c.name]
    i = nums["kcl"] + k
    c.cnode = c.pos_control
    c.dnode = c.neg_control
    # we need to write this BE:
    # ea - eb = r (ec - ed)
    # ea - eb - r ec + r ed = 0
    if c.anode != ground:
        j = nodenum[c.anode]
        assert G[i, j] == 0
        G[i, j] = 1
        G[j, i] = -1
    if c.bnode != ground:
        j = nodenum[c.bnode]
        assert G[i, j] == 0
        G[i, j] = -1
        G[j, i] = 1
    if c.cnode != ground:
        j = nodenum[c.cnode]
        G[i, j] += -r
    if c.dnode != ground:
        j = nodenum[c.dnode]
        G[i, j] += r


def write_VCCS(c, i, j, ground, G, currents, anomnum, nums, nodenum):
    currents.append(c.name)
    g = c.value
    k = anomnum[c.name]
    j = nums["kcl"] + k
    if c.anode != ground:
        i = nodenum[c.anode]
        assert G[i, j] == 0
        G[i, j] = -1
    if c.bnode != ground:
        i = nodenum[c.bnode]
        assert G[i, j] == 0
        G[i, j] = 1
    # we write the branch equation:
    # i_cccs = g (ec - ed)
    # i_cccs - g ec + g ed = 0
    i = j
    G[i, i] = +1
    c.cnode = c.pos_control
    c.dnode = c.neg_control
    if c.cnode != ground:
        j = nodenum[c.cnode]
        G[i, j] = -g
    if c.dnode != ground:
        j = nodenum[c.dnode]
        G[i, j] = +g


def write_CCVS(c, i, j, ground, G, A, currents, anomnum, nums, nodenum):
    r = c.value
    k = anomnum[c.name]
    i = nums["kcl"] + k
    currents.append(c.name)
    c.cnode = c.pos_control
    c.dnode = c.neg_control
    try:
        driver = cs[c.driver]
    except KeyError:
        raise KeyError(f"Driving component {c.driver} not found")
    assert c.cnode is not None
    assert c.dnode is not None
    assert driver is not None
    assert (c.cnode == driver.anode and c.dnode == driver.bnode) or (
        c.cnode == driver.bnode and c.dnode == driver.anode
    )
    if c.anode != ground:
        j = nodenum[c.anode]
        G[i, j] = 1
        G[j, i] = -1
    if c.bnode != ground:
        j = nodenum[c.bnode]
        G[i, j] = -1
        G[j, i] = 1

    # we write the branch equation:
    # v_cccv = r * i_driver
    # ea - eb - r * i_driver = 0
    if driver.type == "R":
        # i_driver = (ec - ed)/R_driver
        if c.cnode != ground:
            j = nodenum[c.cnode]
            G[i, j] = r / driver.value
        if c.dnode != ground:
            j = nodenum[c.dnode]
            G[i, j] = -r / driver.value
    elif driver.type in NODE_TYPES_ANOM:
        j = anomnum[driver.name]
        if driver.anode == c.pos_control:
            assert driver.bnode == c.neg_control
            G[i, j] = -r
        else:
            assert driver.anode == c.neg_control
            assert driver.bnode == c.pos_control
            G[i, j] = r
    elif driver.type == "A":
        A[i] = r * driver.value
    else:
        raise ValueError(f"Unknown component type: {driver.type}")


def write_CCCS(c, i, j, ground, G, A, currents, anomnum, nums, nodenum, components):
    currents.append(c.name)
    g = c.value
    k = anomnum[c.name]
    j = nums["kcl"] + k
    if c.anode != ground:
        i = nodenum[c.anode]
        assert G[i, j] == 0
        G[i, j] = -1
    if c.bnode != ground:
        i = nodenum[c.bnode]
        assert G[i, j] == 0
        G[i, j] = 1
    # we write the branch equation:
    # i_cccs = g * i_driver
    i = j
    assert G[i, i] == 0
    G[i, i] = 1
    try:
        driver = components[c.driver]
    except KeyError:
        raise KeyError(f"Driving component {c.driver} not found")
    # case 1: i_driver is unknown
    if driver.type == "R":
        c.cnode = c.pos_control
        c.dnode = c.neg_control
        assert (c.cnode == driver.anode and c.dnode == driver.bnode) or (
            c.cnode == driver.bnode and c.dnode == driver.anode
        )
        assert c.cnode is not None
        assert c.dnode is not None
        if c.cnode != ground:
            j = nodenum[c.cnode]
            assert G[i, j] == 0
            G[i, j] = +g / driver.value
        if c.dnode != ground:
            j = nodenum[c.dnode]
            assert G[i, j] == 0
            G[i, j] = -g / driver.value
    elif driver.type in NODE_TYPES_ANOM:
        j = anomnum[driver.name]
        if driver.anode == c.pos_control:
            assert driver.bnode == c.neg_control
            G[i, j] = -g
        else:
            assert driver.anode == c.neg_control
            assert driver.bnode == c.pos_control
            G[i, j] = g
    # case 2: i_driver is known
    elif driver.type == "A":
        assert A[i] == 0
        A[i] = g * driver.value
    else:
        raise ValueError(f"Unknown component type: {driver.type}")
