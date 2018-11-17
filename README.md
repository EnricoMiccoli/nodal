# nodal
Solve linear circuits with nodal analysis.

Written in python, uses numpy for the numerical work.

Currently supports linear networks composed of ideal resistors and current generators.

## Example input and output
Input file `netlist.csv`:
```
a1,A,1,1,3
r2,R,1,1,2
r3,R,1,2,3
```
Execute with `$ python nodal.py netlist.csv`.

Printed output:
```
Ground node: 2
e(1) = [ 1.]
e(3) = [-1.]
```
