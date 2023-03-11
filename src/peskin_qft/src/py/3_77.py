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

# QFT, P & S: eq(3.77) (p.51)

for a in range(2):
    for b in range(2):
        for c in range(2):
            for d in range(2):
                LHS = 0
                if a == b and c == d:
                    LHS += 1 # mu = 0
                for mu in range(3):
                    LHS -= sigma[mu][a][b] * sigma[mu][c][d]
                if a > c:
                    RHS_1 = 1
                elif a < c:
                    RHS_1 = -1
                else:
                    RHS_1 = 0
                if b > d:
                    RHS_2 = 1
                elif b < d:
                    RHS_2 = -1
                else:
                    RHS_2 = 0
                if LHS != 2 * RHS_1 * RHS_2:
                    check += 1
print(check)

# QFT, P & S: eq(3.80) (p.52)

for a in range(2):
    for c in range(2):
        for mu in range(3):
            RHS = 0
            LHS = 0
            for b in range(2):
                if a > b:
                    LHS += sigma[mu][b][c]
                elif a < b:
                    LHS -= sigma[mu][b][c]
                if b > c:
                    RHS -= sigma[mu][b][a]
                elif b < c:
                    RHS += sigma[mu][b][a]
            if LHS != RHS:
                check += 1
print(check)
