# Nodal.py
Nodal.py is a simple electrical circuit simulator that uses nodal analysis to solve linear networks made up of resistors and ideal current or voltage sources, both independent and controlled. The numerical work is done by the [numpy](https://www.numpy.org/) package.

## Usage example
Suppose we wanted to solve this circuit:

![Circuit diagram](doc/simple.png)

First we would prepare a file describing the circuit. Each line represents a component, using the format    
`name, type, value, first lead, second lead`.

Input file `netlist.csv`:
```
a1,A,1,1,3
r2,R,1,2,3
r3,R,1,1,2
```

Take notice of the orientation of the current generator: current flows toward the node connected to the first lead. Component type can either be
* `A` for current generators
* `R` for resistors

Default units are ampere and ohm.

We can then execute `$ python nodal.py netlist.csv` to get the list of node potentials (unit is volt).

Printed output:
```
Ground node: 2
e(1) = [ 1.]
e(3) = [-1.]
```

# Complete component specification

| | Type | Value | Lead + | Lead - | Control + | Control - | Driver |
|---|---|---|---|---|---|---|---|
Resistor | R | resistance (ohm) | first | second | NA | NA | NA |
Voltage source | E | voltage (volt) | + voltage | ground | NA | NA | NA |
Current source | A | current (ampere) | + current | ground | NA | NA | NA |
CCCS | CCCS | gain (adimensional) | + current | ground | + driving current | - driving current | component determining the current |
CCVS | CCVS | gain (ohm) | + tension | ground | + driving current | - driving current | component determining the current |
VCCS | VCCS | gain (ohm^-1) | + current | ground | + driving voltage | - driving voltage | NA |
VCVS | VCVS | gain (adimensional) | + voltage | ground | + driving voltage | - driving voltage | NA |
