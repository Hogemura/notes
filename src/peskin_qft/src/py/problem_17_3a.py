#!/usr/bin/env python
import sympy

# QFT, Peskin: Problem 17.3 (a)

t = sympy.Symbol("\\theta")

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

def slashed(vector):
    # vector^mu
    sl = vector[0] * gamma0 - vector[1] * gamma1 - vector[2] * gamma2 - vector[3] * gamma3
    return sl
def prod(vec1, vec2):
    ip = vec1[0] * vec2[0] - vec1[1] * vec2[1] - vec1[2] * vec2[2] - vec1[3] * vec2[3]
    return ip

# Dirac spinors
u_L = sympy.Matrix([0, 1, 0, 0])
v_L = sympy.Matrix([1, 0, 0, 0])
# define CONTRAVARIANT vectors
# momentums
p_1 = sympy.Matrix([1, 0, 0, 1])
p_2 = sympy.Matrix([1, 0, 0, -1])
k_1 = sympy.Matrix([1, sympy.sin(t), 0, sympy.cos(t)])
k_2 = sympy.Matrix([1, -sympy.sin(t), 0, -sympy.cos(t)])
# outgoing gluons
e_1Lc = sympy.Matrix([0, sympy.cos(t), 1j, -sympy.sin(t)])
e_1Rc = sympy.Matrix([0, sympy.cos(t), -1j, -sympy.sin(t)])
e_2Lc = sympy.Matrix([0, -sympy.cos(t), 1j, sympy.sin(t)])
e_2Rc = sympy.Matrix([0, -sympy.cos(t), -1j, sympy.sin(t)])

# q_L ~q_L -> g_R g_R
# t channel
iM = (v_L.T * gamma0 * slashed(e_2Rc) * (slashed(p_1) - slashed(k_1)) * slashed(e_1Rc) * u_L)[0, 0]
print(sympy.nsimplify(sympy.simplify(iM)))

# u channel
iM = (v_L.T * gamma0 * slashed(e_1Rc) * (slashed(p_1) - slashed(k_2)) * slashed(e_2Rc) * u_L)[0, 0]
print(sympy.nsimplify(sympy.simplify(iM)))

# 3 boson vertex
iM = (v_L.T * gamma0  * (prod(e_1Rc, e_2Rc) * slashed(k_2 - k_1) + slashed(e_2Rc) * prod(e_1Rc, - k_1 - 2*k_2) + slashed(e_1Rc) * prod(e_2Rc, 2*k_1 - k_2)) * u_L)[0, 0]
print(sympy.nsimplify(sympy.simplify(iM)))

# q_L ~q_L -> g_R g_L
# t channel
iM = (v_L.T * gamma0 * slashed(e_2Lc) * (slashed(p_1) - slashed(k_1)) * slashed(e_1Rc) * u_L)[0, 0]
print(sympy.nsimplify(sympy.simplify(iM)))

# u channel
iM = (v_L.T * gamma0 * slashed(e_1Rc) * (slashed(p_1) - slashed(k_2)) * slashed(e_2Lc) * u_L)[0, 0]
print(sympy.nsimplify(sympy.simplify(iM)))

# 3 boson vertex
iM = (v_L.T * gamma0  * (prod(e_1Rc, e_2Lc) * slashed(k_2 - k_1) + slashed(e_2Lc) * prod(e_1Rc, - k_1 - 2*k_2) + slashed(e_1Rc) * prod(e_2Lc, 2*k_1 - k_2)) * u_L)[0, 0]
print(sympy.nsimplify(sympy.simplify(iM)))

# q_L ~q_L -> g_L g_R
# t channel
iM = (v_L.T * gamma0 * slashed(e_2Rc) * (slashed(p_1) - slashed(k_1)) * slashed(e_1Lc) * u_L)[0, 0]
print(sympy.nsimplify(sympy.simplify(iM)))

# u channel
iM = (v_L.T * gamma0 * slashed(e_1Lc) * (slashed(p_1) - slashed(k_2)) * slashed(e_2Rc) * u_L)[0, 0]
print(sympy.nsimplify(sympy.simplify(iM)))

# 3 boson vertex
iM = (v_L.T * gamma0  * (prod(e_1Lc, e_2Rc) * slashed(k_2 - k_1) + slashed(e_2Rc) * prod(e_1Lc, - k_1 - 2*k_2) + slashed(e_1Lc) * prod(e_2Rc, 2*k_1 - k_2)) * u_L)[0, 0]
print(sympy.nsimplify(sympy.simplify(iM)))
