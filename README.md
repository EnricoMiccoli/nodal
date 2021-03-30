# Nodal.py
Nodal.py is a simple electrical circuit simulator that uses nodal analysis to solve linear networks made up of resistors and ideal current or voltage sources, both independent and controlled.

The intendend use case of nodal is through the two provided commandline scripts. Even so, the package has been carefully written so that it might be used to programmatically simulate circuits from inside another application.

The numerical work is done by the [numpy](https://www.numpy.org/) package. For larger circuits, the user may choose to use sparce matrices, in which case computation is powered by [scipy](https://docs.scipy.org/doc/scipy/reference/sparse.linalg.html) instead.




### Table of contents
* [How to install](#installation)
* [Command line scripts](#command-line-scripts)
    * [Defining the input](#defining-the-input)
    * [Solving circuits: nodal-solver](#circuit-solver)
    * [Calculating resistance: nodal-resistance](#resistance-calculator)
    * [Component specification](#component-specification)
* [Usage examples](#usage-example)




## Installation
To get the stable release [download it](https://github.com/EnricoMiccoli/nodal/releases/latest) and run
```
$ pip install nodal-1.3.0.whl
```

To try out the latest developments use [`flit`](https://github.com/takluyver/flit) instead:
```
$ git clone https://github.com/EnricoMiccoli/nodal.git
$ cd nodal
$ flit install
```




## Command line scripts

### Defining the input
All scripts read their input from a [netlist](https://en.wikipedia.org/wiki/Netlist) in the csv format. Each line represents a component and contains informations about its electrical values (e.g. the resistance in ohm for a resistor) and all the connections it makes with other components.

Many annotated examples of netlists are available in the `doc/` directory.


### Circuit solver
Usage: `$ nodal-solver netlist.csv`. If the circuit is particularly large, one might use `$ nodal-solver -s netlist.csv` instead to save memory and runtime.

The output is printed in the following format:
```
e(node_name) = node potential (volt)
i(component_name) = current passing through the component (ampere)
```
The ground node is used as a potential reference. It is set at 0 volt.



### Resistance calculator
Usage: `$ nodal-resistance netlist.csv` or `$ nodal-resistance -s netlist.csv`.

This will print the equivalent resistance (unit is ohm) as seen through the nodes named `1` and `g` in the netlist.


### Component specification
| | Type | Value | Lead + | Lead - | Control + | Control - | Driver |
|---|---|---|---|---|---|---|---|
Resistor | R | resistance (ohm) | first | second | NA | NA | NA |
Voltage source | E | voltage (volt) | + voltage | ground | NA | NA | NA |
Current source | A | current (ampere) | + current | ground | NA | NA | NA |
CCCS | CCCS | gain (adimensional) | + current | ground | + driving current | - driving current | component determining the current |
CCVS | CCVS | gain (ohm) | + tension | ground | + driving current | - driving current | component determining the current |
VCCS | VCCS | gain (ohm^-1) | + current | ground | + driving voltage | - driving voltage | NA |
VCVS | VCVS | gain (adimensional) | + voltage | ground | + driving voltage | - driving voltage | NA |
OPAMP, simulated | OPMODEL | feedback resistance (ohm) | output | ground | + terminal | - terminal | NA |

Note that opamps are not supported directly, but they are instead internally simulated with an equivalent circuit. For example this is the circuit generated to model a voltage buffer:

![Voltage buffer circuit diagram](doc/buffer.png)




## Usage example
Suppose we wanted to solve this circuit:

![Circuit diagram](doc/simple.png)

First we would prepare a file describing the circuit. Each line represents a component, using the format    
`name, type, value, first lead, second lead`.

Input file `netlist.csv`:
```
a1,A,1,1,3
r1,R,1,1,2
r2,R,1,2,3
```

Take notice of the orientation of the current generator: current flows toward the node connected to the first lead. Component type can either be
* `A` for current generators
* `R` for resistors

Default units are ampere and ohm.

We can then execute `$ nodal-solver netlist.csv` to get the list of node potentials (unit is volt).

Printed output:
```
Ground node: 1
e(2)    = -1.0
e(3)    = -2.0
```

### A more advanced example

![Circuit diagram](doc/1.6.1.png)

_Example circuit from: [Esercizi da temi d'esame di Elettrotecnica](https://damore.faculty.polimi.it/download/temiDEsame.pdf)_

Input file `1.6.1.csv`:
```
r1,R,2,1,4
r2,R,2,1,g
r3,R,0.5,1,2
e1,E,8,4,g
a1,A,4,1,2
d1,CCCS,2,2,g,1,g,r2
```

Output:
```
Ground node: g
e(1)    = 2.0
e(2)    = -1.0
e(4)    = 8.0
i(d1)   = -2.0
i(e1)   = 3.0
```

### A more contrived example
[What is the resistance between two points that are a knight's move away on an infinite grid of 1 ohm resistors?](https://enricomiccoli.github.io/2019/03/20/xkcd-356-infinite-grid-resistors.html)

