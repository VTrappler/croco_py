#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np


def get_angle(latu, lonu, spheroid='wgs84'):
    eps = 1e-22
    if spheroid == 'sphere':
        A = 6371000.0
        B = A
        E = np.sqrt(A * A - B * B)
    elif spheroid == 'clarke66':
        A = 6378206.4
        B = 6356583.8
        E = np.sqrt(A * A - B * B) / A
    elif spheroid == 'iau73':
        A = 6378160.
        B = 6356774.516
        E = np.sqrt(A * A - B * B) / A
    elif spheroid == 'wgs84':
        A = 6378137.
        E = 0.081819191
        B = np.sqrt(A**2 - (A * E)**2)

    EPS = E * E / (1 - E * E)

    latu = latu * np.pi / 180
    lonu = lonu * np.pi / 180

    latu[latu == 0] = eps
    M, L = latu.shape
    PHI1 = latu[:, :-1]
    XLAM1 = lonu[:, :-1]
    PHI2 = latu[:, 1:]
    XLAM2 = lonu[:, 1:]

    PHI2[PHI1 == PHI2] = 1e-14 + PHI2[PHI1 == PHI2]
    XLAM2[XLAM1 == XLAM2] = 1e-14 + XLAM2[XLAM1 == XLAM2]

    xnu1, xnu2 = [A / np.sqrt(1.0 - (E * np.sin(phi))**2)
                  for phi in [PHI1, PHI2]]
    TPSI2 = (1.0 - E**2) * np.tan(PHI2) + E**2 * xnu1 * np.sin(PHI1) / (xnu2 * np.cos(PHI2))

    DLAM = XLAM2 - XLAM1
    CTA12 = (np.cos(PHI1) * TPSI2 - np.sin(PHI1) * np.cos(DLAM)) / np.sin(DLAM)
    azim = np.arctan(1 / CTA12)

    DLAM2 = (np.abs(DLAM) < np.pi) * DLAM \
        + (DLAM >= np.pi) * (-2 * np.pi + DLAM) \
        + (DLAM <= -np.pi) * (2 * np.pi + DLAM)

    azim = azim + (azim < -np.pi) * 2 * np.pi - (azim >= np.pi) * 2 * np.pi
    azim = azim + np.pi * np.sign(-azim) * (np.sign(azim) != np.sign(DLAM2))
    angle = np.empty((azim.shape[0], azim.shape[1] + 2))
    angle[:, 1:-1] = (np.pi / 2.0) - azim
    angle[:, 0] = angle[:, 1]
    angle[:, -1] = angle[:, -2]

    return angle
