#!/usr/bin/env python
import numpy as np

# Gell-Mann matrices
# QFT, Peskin. (20.40)(20.42)

t_1 = np.array([
[0, 1, 0],
[1, 0, 0],
[0, 0, 0]])
t_2 = np.array([
[0, -1j, 0],
[1j, 0, 0],
[0, 0, 0]])
t_3 = np.array([
[1, 0, 0],
[0, -1, 0],
[0, 0, 0]])
t_4 = np.array([
[0, 0, 1],
[0, 0, 0],
[1, 0, 0]])
t_5 = np.array([
[0, 0, -1j],
[0, 0, 0],
[1j, 0, 0]])
t_6 = np.array([
[0, 0, 0],
[0, 0, 1],
[0, 1, 0]])
t_7 = np.array([
[0, 0, 0],
[0, 0, -1j],
[0, 1j, 0]])
t_8 = np.array([
[1, 0, 0],
[0, 1, 0],
[0, 0, -2]]) * 3 ** (-0.5)

ts = [t_1/2, t_2/2, t_3/2, t_4/2, t_5/2, t_6/2, t_7/2, t_8/2]

Phi = np.array([
[1, 0, 0],
[0, 1, 0],
[0, 0, -2]])

for a in range(8):
    for b in range(8):
        M = (ts[a]@Phi - Phi@ts[a])@(ts[b]@Phi - Phi@ts[b])
        if np.trace(M) != 0:
            print('a=', a+1, 'b=', b+1, 'm^2=', np.trace(M)*(-2))

Phi = np.array([
[1, 0, 0],
[0, -1, 0],
[0, 0, 0]])

for a in range(8):
    for b in range(8):
        M = (ts[a]@Phi - Phi@ts[a])@(ts[b]@Phi - Phi@ts[b])
        if np.trace(M) != 0:
            print('a=', a+1, 'b=', b+1, 'm^2=', np.trace(M)*(-2))
