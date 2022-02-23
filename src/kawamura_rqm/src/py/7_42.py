#!/usr/bin/env python
import numpy as np

zero = np.zeros([2, 2])
identity = np.identity(2)
metric = np.array([
[1, 0, 0, 0],
[0, -1, 0, 0],
[0, 0, -1, 0],
[0, 0, 0, -1]])

sigma_1 = np.array([[0, 1], [1, 0]])
sigma_2 = np.array([[0, -1j], [1j, 0]])
sigma_3 = np.array([[1, 0], [0, -1]])

sigma = [sigma_1, sigma_2, sigma_3]

gamma_0 = np.array([
[1, 0, 0, 0],
[0, 1, 0, 0],
[0, 0, -1, 0],
[0, 0, 0, -1]], dtype=np.complex128)
gamma_1 = np.concatenate([np.concatenate([zero, -sigma_1]), np.concatenate([sigma_1, zero])], axis=1)
gamma_2 = np.concatenate([np.concatenate([zero, -sigma_2]), np.concatenate([sigma_2, zero])], axis=1)
gamma_3 = np.concatenate([np.concatenate([zero, -sigma_3]), np.concatenate([sigma_3, zero])], axis=1)

gammaM_0 = np.concatenate([np.concatenate([zero, sigma_2]), np.concatenate([sigma_2, zero])], axis=1)
gammaM_1 = np.concatenate([np.concatenate([1j * sigma_3, zero]), np.concatenate([zero, 1j * sigma_3])], axis=1)
gammaM_2 = np.concatenate([np.concatenate([zero, sigma_2]), np.concatenate([-sigma_2, zero])], axis=1)
gammaM_3 = np.concatenate([np.concatenate([-1j * sigma_1, zero]), np.concatenate([zero, -1j * sigma_1])], axis=1)

gammaM_5 = 1j * gammaM_0 @ gammaM_1 @ gammaM_2 @ gammaM_3

check = 0
bool = gammaM_5 == np.concatenate([np.concatenate([sigma_2, zero]), np.concatenate([zero, -sigma_2])], axis=1)
if bool.all() == False:
    check += 1
print(check)

Unit = np.concatenate([np.concatenate([identity, sigma_2]), np.concatenate([sigma_2, -identity])], axis=1)

# Dirac
gammas_D = [gamma_0, gamma_1, gamma_2, gamma_3]

# Majorana
gammas_M = [gammaM_0, gammaM_1, gammaM_2, gammaM_3]

# RQM, Kawamura: eq(7.42) (p.106)

for i in range(4):
    bool = gammas_M[i] == Unit @ gammas_D[i] @ Unit * 1 / 2
    if bool.all() == False:
        check += 1
print(check)

# QFT, Sakamoto: check 6.12 (p.181)

for i in range(4):
    for k in range(4):
        bool = gammas_M[i] @ gammas_M[k] + gammas_M[k] @ gammas_M[i] == 2 * metric[i][k] * np.identity(4)
        if bool.all() == False:
            check += 1
print(check)

bool = np.conjugate(gammas_M[0].T) == gammas_M[0]
if bool.all() == False:
    check += 1

for i in range(1, 4):
    bool = np.conjugate(gammas_M[i].T) == -gammas_M[i]
    if bool.all() == False:
        check += 1
print(check)

for i in range(4):
    bool = np.conjugate(gammas_M[i]) == -gammas_M[i]
    if bool.all() == False:
        check += 1
print(check)

for i in range(4):
    bool = gammas_M[0] @ gammas_M[i].T @ gammas_M[0] == -gammas_M[i]
    if bool.all() == False:
        check += 1
print(check)

# QFT, Sakamoto: check 7.7 (p.199)

T_1 = np.array([
[0, 1, 0],
[1, 0, 1],
[0, 1, 0]
], dtype=np.complex128)

T_2 = np.array([
[0, -1j, 0],
[1j, 0, -1j],
[0, 1j, 0]
], dtype=np.complex128)

T_3 = np.array([
[2, 0, 0],
[0, 0, 0],
[0, 0, -2]
], dtype=np.complex128)

bool =  T_1 @ T_2 - T_2 @ T_1 == 1j * T_3
if bool.all() == False:
    check += 1
print(check)

bool =  T_2 @ T_3 - T_3 @ T_2 == 2j * T_1
if bool.all() == False:
    check += 1
print(check)

bool =  T_3 @ T_1 - T_1 @ T_3 == 2j * T_2
if bool.all() == False:
    check += 1
print(check)
