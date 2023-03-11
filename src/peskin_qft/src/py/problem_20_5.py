#!/usr/bin/env python
import sympy

# QFT, Peskin: Problem 20.5 (b)

# upper component
p_c1 = sympy.Symbol("\\phi_1^+")
p_c2 = sympy.Symbol("\\phi_2^+")
# upper component (complex conj.)
p_c1c = sympy.Symbol("\\phi_1^-")
p_c2c = sympy.Symbol("\\phi_2^-")
# expectation
v_1 = sympy.Symbol("v_1")
v_2 = sympy.Symbol("v_2")
# lower Re
r_1 = sympy.Symbol("\\rho_1")
r_2 = sympy.Symbol("\\rho_2")
# lower Im
e_1 = sympy.Symbol("\\eta_1")
e_2 = sympy.Symbol("\\eta_2")
# coef
l_1, l_2, l_3, l_4, l_5 = sympy.symbols("\\lambda_1, \\lambda_2, \\lambda_3, \\lambda_4, \\lambda_5")


hd_1 = sympy.Matrix([p_c1, (v_1 + r_1 + sympy.I * e_1) / sympy.sqrt(2)])
hd_2 = sympy.Matrix([p_c2, (v_2 + r_2 + sympy.I * e_2) / sympy.sqrt(2)])
hd_1c = sympy.Matrix([p_c1c, (v_1 + r_1 - sympy.I * e_1) / sympy.sqrt(2)])
hd_2c = sympy.Matrix([p_c2c, (v_2 + r_2 - sympy.I * e_2) / sympy.sqrt(2)])

p1p1 = hd_1c.transpose() * hd_1
p2p2 = hd_2c.transpose() * hd_2
p1p2 = hd_1c.transpose() * hd_2
p2p1 = hd_2c.transpose() * hd_1

m_1 = sympy.Symbol("\\mu_1")
m_2 = sympy.Symbol("\\mu_2")
V = (- m_1**2 * p1p1 - m_2**2 * p2p2 + l_1 * p1p1**2 + l_2 * p2p2**2 + l_3 * p1p1 * p2p2 + l_4 * p1p2 * p2p1 + l_5 * (p1p2**2 + p2p1**2))[0]

dV_dv1 = sympy.diff(V, v_1)
print("\\frac{\partial V}{\partial v_1}(\phi_a = \langle\phi_a\rangle) =", sympy.latex(sympy.expand(dV_dv1.subs([(p_c1, 0), (p_c2, 0), (p_c1c, 0), (p_c2c, 0), (r_1, 0), (r_2, 0), (e_1, 0), (e_2, 0)]))))
dV_dv2 = sympy.diff(V, v_2)
print("\\frac{\partial V}{\partial v_2}(\phi_a = \langle\phi_a\rangle) =",sympy.latex(sympy.expand(dV_dv2.subs([(p_c1, 0), (p_c2, 0), (p_c1c, 0), (p_c2c, 0), (r_1, 0), (r_2, 0), (e_1, 0), (e_2, 0)]))))

m_1 = sympy.sqrt(l_1 * v_1**2 + (l_3/2 + l_4/2 + l_5) * v_2**2)
m_2 = sympy.sqrt(l_2 * v_2**2 + (l_3/2 + l_4/2 + l_5) * v_1**2)
V = (- m_1**2 * p1p1 - m_2**2 * p2p2 + l_1 * p1p1**2 + l_2 * p2p2**2 + l_3 * p1p1 * p2p2 + l_4 * p1p2 * p2p1 + l_5 * (p1p2**2 + p2p1**2))[0]

# phi
ddV11 = sympy.diff(sympy.diff(V, p_c1c), p_c1)
print("\\frac{\\partial^2V}{\\partial\\phi_1^-\\partial\\phi_1^+}(0) =", sympy.latex(sympy.expand(ddV11.subs([(p_c1, 0), (p_c2, 0), (p_c1c, 0), (p_c2c, 0), (r_1, 0), (r_2, 0), (e_1, 0), (e_2, 0)]))))

ddV12 = sympy.diff(sympy.diff(V, p_c1c), p_c2)
print("\\frac{\\partial^2V}{\\partial\\phi_1^-\\partial\\phi_2^+}(0) =", sympy.latex(sympy.expand(ddV12.subs([(p_c1, 0), (p_c2, 0), (p_c1c, 0), (p_c2c, 0), (r_1, 0), (r_2, 0), (e_1, 0), (e_2, 0)]))))

ddV21 = sympy.diff(sympy.diff(V, p_c2c), p_c1)
print("\\frac{\\partial^2V}{\\partial\\phi_2^-\\partial\\phi_1^+}(0) =", sympy.latex(sympy.expand(ddV21.subs([(p_c1, 0), (p_c2, 0), (p_c1c, 0), (p_c2c, 0), (r_1, 0), (r_2, 0), (e_1, 0), (e_2, 0)]))))

ddV22 = sympy.diff(sympy.diff(V, p_c2c), p_c2)
print("\\frac{\\partial^2V}{\\partial\\phi_2^-\\partial\\phi_2^+}(0) =", sympy.latex(sympy.expand(ddV22.subs([(p_c1, 0), (p_c2, 0), (p_c1c, 0), (p_c2c, 0), (r_1, 0), (r_2, 0), (e_1, 0), (e_2, 0)]))))

# rho
ddV11 = sympy.diff(sympy.diff(V, r_1), r_1)
print("\\frac{\\partial^2V}{\\partial\\rho_1\\partial\\rho_1}(0) =", sympy.latex(sympy.expand(ddV11.subs([(p_c1, 0), (p_c2, 0), (p_c1c, 0), (p_c2c, 0), (r_1, 0), (r_2, 0), (e_1, 0), (e_2, 0)]))))

ddV12 = sympy.diff(sympy.diff(V, r_1), r_2)
print("\\frac{\\partial^2V}{\\partial\\rho_1\\partial\\rho_2}(0) =", sympy.latex(sympy.expand(ddV12.subs([(p_c1, 0), (p_c2, 0), (p_c1c, 0), (p_c2c, 0), (r_1, 0), (r_2, 0), (e_1, 0), (e_2, 0)]))))

ddV22 = sympy.diff(sympy.diff(V, r_2), r_2)
print("\\frac{\\partial^2V}{\\partial\\rho_2\\partial\\rho_2}(0) =", sympy.latex(sympy.expand(ddV22.subs([(p_c1, 0), (p_c2, 0), (p_c1c, 0), (p_c2c, 0), (r_1, 0), (r_2, 0), (e_1, 0), (e_2, 0)]))))

# eta
ddV11 = sympy.diff(sympy.diff(V, e_1), e_1)
print("\\frac{\\partial^2V}{\\partial\\eta_1\\partial\\eta_1}(0) =", sympy.latex(sympy.expand(ddV11.subs([(p_c1, 0), (p_c2, 0), (p_c1c, 0), (p_c2c, 0), (r_1, 0), (r_2, 0), (e_1, 0), (e_2, 0)]))))

ddV12 = sympy.diff(sympy.diff(V, e_1), e_2)
print("\\frac{\\partial^2V}{\\partial\\eta_1\\partial\\eta_2}(0) =", sympy.latex(sympy.expand(ddV12.subs([(p_c1, 0), (p_c2, 0), (p_c1c, 0), (p_c2c, 0), (r_1, 0), (r_2, 0), (e_1, 0), (e_2, 0)]))))

ddV22 = sympy.diff(sympy.diff(V, e_2), e_2)
print("\\frac{\\partial^2V}{\\partial\\eta_2\\partial\\eta_2}(0) =", sympy.latex(sympy.expand(ddV22.subs([(p_c1, 0), (p_c2, 0), (p_c1c, 0), (p_c2c, 0), (r_1, 0), (r_2, 0), (e_1, 0), (e_2, 0)]))))
