#!/usr/bin/env python
import sympy

E, p_x, p_y, p_z = sympy.symbols("E, p_x, p_y, p_z")
m = sympy.sqrt(E - p_x**2 - p_y**2 - p_z**2)

D = sympy.Matrix(
[[-m, 0, E-p_z, -p_x+sympy.I*p_y],
[0, -m, -p_x-sympy.I*p_y, E+p_z],
[E+p_z, p_x-sympy.I*p_y, -m, 0],
[p_x+sympy.I*p_y, E-p_z, 0, -m]])

print(D.eigenvals())
print(D.eigenvects())

# p*sigma
ps = sympy.Matrix(
[[E-p_z, -p_x+sympy.I*p_y],
[-p_x-sympy.I*p_y, E+p_z]])

print(ps.eigenvals())

# p*~sigma
pbs = sympy.Matrix(
[[E+p_z, +p_x-sympy.I*p_y],
[p_x+sympy.I*p_y, E-p_z]])

print(pbs.eigenvals())
