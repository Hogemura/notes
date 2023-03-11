#!/usr/bin/env python
import sympy
from sympy.simplify.fu import TR5

# QFT, Peskin: Problem 17.3 (b)

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
p_1 = sympy.Matrix([1, 0, 0, 1])
p_2 = sympy.Matrix([1, 0, 0, -1])
k_1 = sympy.Matrix([1, sympy.sin(t), 0, sympy.cos(t)])
k_2 = sympy.Matrix([1, -sympy.sin(t), 0, -sympy.cos(t)])
# incoming gluons
e_p1L = sympy.Matrix([0, 1, -1j, 0])
e_p1R = sympy.Matrix([0, 1, 1j, 0])
e_p2L = sympy.Matrix([0, -1, -1j, 0])
e_p2R = sympy.Matrix([0, -1, 1j, 0])
# outgoing gluons
e_k1Lc = sympy.Matrix([0, sympy.cos(t), 1j, -sympy.sin(t)])
e_k1Rc = sympy.Matrix([0, sympy.cos(t), -1j, -sympy.sin(t)])
e_k2Lc = sympy.Matrix([0, -sympy.cos(t), 1j, sympy.sin(t)])
e_k2Rc = sympy.Matrix([0, -sympy.cos(t), -1j, sympy.sin(t)])

"""
g_R g_R -> g_R g_R
"""
# t channel
V1 = three_vertex(-k_1, p_1, k_1-p_1)
V2 = three_vertex(-k_2, k_2-p_2, p_2)

factor = 0
for m in range(4):
    for n in range(4):
        for r in range(4):
            for s in range(4):
                for l in range(4):
                    for k in range(4):
                        factor += e_p1R[m] * e_p2R[n] * V1[r][m][k] * V2[s][l][n] * g[k, l] * e_k1Rc[r] * e_k2Rc[s]

print(sympy.expand_trig(sympy.simplify(sympy.nsimplify(factor))))

# u channel
V1 = three_vertex(-k_2, p_1, k_2-p_1)
V2 = three_vertex(-k_1, k_1-p_2, p_2)

factor = 0
for m in range(4):
    for n in range(4):
        for r in range(4):
            for s in range(4):
                for l in range(4):
                    for k in range(4):
                        factor += e_p1R[m] * e_p2R[n] * V1[s][m][k] * V2[r][l][n] * g[k, l] * e_k1Rc[r] * e_k2Rc[s]

print(sympy.expand_trig(sympy.simplify(sympy.nsimplify(factor))))

# s channel
V1 = three_vertex(-p_1-p_2, p_1, p_2)
V2 = three_vertex(k_1+k_2, -k_1, -k_2)

factor = 0
for m in range(4):
    for n in range(4):
        for r in range(4):
            for s in range(4):
                for l in range(4):
                    for k in range(4):
                        factor += e_p1R[m] * e_p2R[n] * V1[k][m][n] * V2[l][r][s] * g[l, k] * e_k1Rc[r] * e_k2Rc[s]

print(sympy.expand_trig(sympy.simplify(sympy.nsimplify(factor))))

# 4 bsoon
t1 = 0
t2 = 0
t3 = 0
for m in range(4):
    for n in range(4):
        for r in range(4):
            for s in range(4):
                t1 += e_p1R[m] * e_p2R[n] * (g[m, r] * g[n, s] - g[m, s] * g[n, r]) * e_k1Rc[r] * e_k2Rc[s]
                t2 += e_p1R[m] * e_p2R[n] * (g[m, n] * g[r, s] - g[m, s] * g[n, r]) * e_k1Rc[r] * e_k2Rc[s]
                t3 += e_p1R[m] * e_p2R[n] * (g[m, n] * g[r, s] - g[m, r] * g[n, s]) * e_k1Rc[r] * e_k2Rc[s]

print("f^{abe} f^{cde}:", sympy.expand_trig(sympy.simplify(sympy.nsimplify(t1))))
print("f^{ace} f^{bde}:", TR5(sympy.expand_trig(sympy.simplify(sympy.nsimplify(t2)))))
print("f^{ade} f^{bce}:", TR5(sympy.expand_trig(sympy.simplify(sympy.nsimplify(t3)))), "\n")

"""
g_R g_R -> g_R g_L
"""
# t channel
V1 = three_vertex(-k_1, p_1, k_1-p_1)
V2 = three_vertex(-k_2, k_2-p_2, p_2)

factor = 0
for m in range(4):
    for n in range(4):
        for r in range(4):
            for s in range(4):
                for l in range(4):
                    for k in range(4):
                        factor += e_p1R[m] * e_p2R[n] * V1[r][m][k] * V2[s][l][n] * g[k, l] * e_k1Rc[r] * e_k2Lc[s]

print(sympy.expand_trig(sympy.simplify(sympy.nsimplify(factor))))

# u channel
V1 = three_vertex(-k_2, p_1, k_2-p_1)
V2 = three_vertex(-k_1, k_1-p_2, p_2)

factor = 0
for m in range(4):
    for n in range(4):
        for r in range(4):
            for s in range(4):
                for l in range(4):
                    for k in range(4):
                        factor += e_p1R[m] * e_p2R[n] * V1[s][m][k] * V2[r][l][n] * g[k, l] * e_k1Rc[r] * e_k2Lc[s]

print(sympy.expand_trig(sympy.simplify(sympy.nsimplify(factor))))

# s channel
V1 = three_vertex(-p_1-p_2, p_1, p_2)
V2 = three_vertex(k_1+k_2, -k_1, -k_2)

factor = 0
for m in range(4):
    for n in range(4):
        for r in range(4):
            for s in range(4):
                for l in range(4):
                    for k in range(4):
                        factor += e_p1R[m] * e_p2R[n] * V1[k][m][n] * V2[l][r][s] * g[l, k] * e_k1Rc[r] * e_k2Lc[s]

print(sympy.expand_trig(sympy.simplify(sympy.nsimplify(factor))))

# 4 bsoon
t1 = 0
t2 = 0
t3 = 0
for m in range(4):
    for n in range(4):
        for r in range(4):
            for s in range(4):
                t1 += e_p1R[m] * e_p2R[n] * (g[m, r] * g[n, s] - g[m, s] * g[n, r]) * e_k1Rc[r] * e_k2Lc[s]
                t2 += e_p1R[m] * e_p2R[n] * (g[m, n] * g[r, s] - g[m, s] * g[n, r]) * e_k1Rc[r] * e_k2Lc[s]
                t3 += e_p1R[m] * e_p2R[n] * (g[m, n] * g[r, s] - g[m, r] * g[n, s]) * e_k1Rc[r] * e_k2Lc[s]

print("f^{abe} f^{cde}:", sympy.expand_trig(sympy.simplify(sympy.nsimplify(t1))))
print("f^{ace} f^{bde}:", TR5(sympy.expand_trig(sympy.simplify(sympy.nsimplify(t2)))))
print("f^{ade} f^{bce}:", TR5(sympy.expand_trig(sympy.simplify(sympy.nsimplify(t3)))), "\n")

"""
g_R g_R -> g_L g_L
"""
# t channel
V1 = three_vertex(-k_1, p_1, k_1-p_1)
V2 = three_vertex(-k_2, k_2-p_2, p_2)

factor = 0
for m in range(4):
    for n in range(4):
        for r in range(4):
            for s in range(4):
                for l in range(4):
                    for k in range(4):
                        factor += e_p1R[m] * e_p2R[n] * V1[r][m][k] * V2[s][l][n] * g[k, l] * e_k1Lc[r] * e_k2Lc[s]

print(sympy.expand_trig(sympy.simplify(sympy.nsimplify(factor))))

# u channel
V1 = three_vertex(-k_2, p_1, k_2-p_1)
V2 = three_vertex(-k_1, k_1-p_2, p_2)

factor = 0
for m in range(4):
    for n in range(4):
        for r in range(4):
            for s in range(4):
                for l in range(4):
                    for k in range(4):
                        factor += e_p1R[m] * e_p2R[n] * V1[s][m][k] * V2[r][l][n] * g[k, l] * e_k1Lc[r] * e_k2Lc[s]

print(sympy.expand_trig(sympy.simplify(sympy.nsimplify(factor))))

# s channel
V1 = three_vertex(-p_1-p_2, p_1, p_2)
V2 = three_vertex(k_1+k_2, -k_1, -k_2)

factor = 0
for m in range(4):
    for n in range(4):
        for r in range(4):
            for s in range(4):
                for l in range(4):
                    for k in range(4):
                        factor += e_p1R[m] * e_p2R[n] * V1[k][m][n] * V2[l][r][s] * g[l, k] * e_k1Lc[r] * e_k2Lc[s]

print(sympy.expand_trig(sympy.simplify(sympy.nsimplify(factor))))

# 4 bsoon
t1 = 0
t2 = 0
t3 = 0
for m in range(4):
    for n in range(4):
        for r in range(4):
            for s in range(4):
                t1 += e_p1R[m] * e_p2R[n] * (g[m, r] * g[n, s] - g[m, s] * g[n, r]) * e_k1Lc[r] * e_k2Lc[s]
                t2 += e_p1R[m] * e_p2R[n] * (g[m, n] * g[r, s] - g[m, s] * g[n, r]) * e_k1Lc[r] * e_k2Lc[s]
                t3 += e_p1R[m] * e_p2R[n] * (g[m, n] * g[r, s] - g[m, r] * g[n, s]) * e_k1Lc[r] * e_k2Lc[s]

print("f^{abe} f^{cde}:", sympy.expand_trig(sympy.simplify(sympy.nsimplify(t1))))
print("f^{ace} f^{bde}:", TR5(sympy.expand_trig(sympy.simplify(sympy.nsimplify(t2)))))
print("f^{ade} f^{bce}:", TR5(sympy.expand_trig(sympy.simplify(sympy.nsimplify(t3)))), "\n")
