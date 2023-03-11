#!/usr/bin/env python
import sympy

# QFT, Peskin: Problem 17.4

z = sympy.Symbol("z")
p = sympy.Symbol("p")
x = sympy.Symbol("x") # ratio: p_perp/p

gamma0 = sympy.Matrix([
[0, 0, 1, 0],
[0, 0, 0, 1],
[1, 0, 0, 0],
[0, 1, 0, 0]
])

gamma1 = sympy.Matrix([
[0, 0, 0, 1],
[0, 0, 1, 0],
[0, -1, 0, 0],
[-1, 0, 0, 0]
])

gamma2 = sympy.Matrix([
[0, 0, 0, -1j],
[0, 0, 1j, 0],
[0, 1j, 0, 0],
[-1j, 0, 0, 0]
])

gamma3 = sympy.Matrix([
[0, 0, 1, 0],
[0, 0, 0, -1],
[-1, 0, 0, 0],
[0, 1, 0, 0]
])

g = sympy.Matrix([
[1, 0, 0, 0],
[0, -1, 0, 0],
[0, 0, -1, 0],
[0, 0, 0, -1]
])

d = sympy.Matrix([
[1, 0, 0, 0],
[0, 1, 0, 0],
[0, 0, 1, 0],
[0, 0, 0, 1]
])

def prod(vec1, vec2):
    ip = vec1[0] * vec2[0] - vec1[1] * vec2[1] - vec1[2] * vec2[2] - vec1[3] * vec2[3]
    return ip

def three_vertex(k, p, q):
    if k+p+q != sympy.Matrix([0, 0, 0 , 0]):
        print('momentums not conserved !')
    vertex = [
    [["", "", "", ""], ["", "", "", ""], ["", "", "", ""], ["", "", "", ""]],
    [["", "", "", ""], ["", "", "", ""], ["", "", "", ""], ["", "", "", ""]],
    [["", "", "", ""], ["", "", "", ""], ["", "", "", ""], ["", "", "", ""]],
    [["", "", "", ""], ["", "", "", ""], ["", "", "", ""], ["", "", "", ""]]]
    for m in range(4):
        for n in range(4):
            for r in range(4):
                vertex[m][n][r] = g[m, n]*(k-p)[r] + g[n, r]*(p-q)[m] + g[r, m]*(q-k)[n]
    return vertex

# define CONTRAVARIANT vectors
p = sympy.Matrix([1, 0, 0, 1])
q = sympy.Matrix([z, x, 0, z])
k = sympy.Matrix([1-z, -x, 0, 1-z])

# incoming gluons
e_pL = sympy.Matrix([0, 1, -1j, 0])
e_pR = sympy.Matrix([0, 1, 1j, 0])

# outgoing gluons
e_qLc = sympy.Matrix([0, 1, 1j, -x/z])
e_qRc = sympy.Matrix([0, 1, -1j, -x/z])
e_kLc = sympy.Matrix([0, 1, 1j, x/(1-z)])
e_kRc = sympy.Matrix([0, 1, -1j, x/(1-z)])

"""
g_R -> g_R g_R
"""
V = three_vertex(p, -q, -k)

factor = 0
for m in range(4):
    for n in range(4):
        for r in range(4):
            factor -= e_pR[m] * V[m][n][r] * e_qRc[n] * e_kRc[r]

print("f^{abc}:", sympy.simplify(sympy.nsimplify(factor)))

"""
g_R -> g_R g_L
"""
V = three_vertex(p, -q, -k)

factor = 0
for m in range(4):
    for n in range(4):
        for r in range(4):
            factor -= e_pR[m] * V[m][n][r] * e_qRc[n] * e_kLc[r]

print("f^{abc}:", sympy.simplify(sympy.nsimplify(factor)))

"""
g_R -> g_L g_R
"""
V = three_vertex(p, -q, -k)

factor = 0
for m in range(4):
    for n in range(4):
        for r in range(4):
            factor -= e_pR[m] * V[m][n][r] * e_qLc[n] * e_kRc[r]

print("f^{abc}:", sympy.simplify(sympy.nsimplify(factor)))

"""
g_R -> g_L g_L
"""
V = three_vertex(p, -q, -k)

factor = 0
for m in range(4):
    for n in range(4):
        for r in range(4):
            factor -= e_pR[m] * V[m][n][r] * e_qLc[n] * e_kLc[r]

print("f^{abc}:", sympy.simplify(sympy.nsimplify(factor)))
