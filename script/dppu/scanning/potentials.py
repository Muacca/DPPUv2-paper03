"""
Potential Functions for Parameter Scanning
==========================================

Pre-computed potential functions for each topology × NY variant.
"""

import numpy as np

PI = np.pi
KAPPA = 1.0
L = 1.0


# S³×S¹ Potentials
def V_S3_FULL(r, V_param, eta, theta):
    if r <= 0:
        return 1e50
    term1 = V_param**2 * r**2
    term2 = 6 * V_param * KAPPA**2 * theta * (eta - 4) * r
    term3 = 9 * eta**2 + 72 * eta + 108
    return (2 * PI**2 * L * r * (term1 + term2 + term3)) / (3 * KAPPA**2)


def V_S3_TT(r, V_param, eta, theta):
    if r <= 0:
        return 1e50
    term1 = V_param**2 * r**2
    term2 = 12 * V_param * eta * KAPPA**2 * theta * r
    term3 = 9 * eta**2 + 72 * eta + 108
    return (2 * PI**2 * L * r * (term1 + term2 + term3)) / (3 * KAPPA**2)


def V_S3_REE(r, V_param, eta, theta):
    if r <= 0:
        return 1e50
    term1 = V_param**2 * r**2
    term2 = 6 * V_param * KAPPA**2 * theta * (eta + 4) * r
    term3 = 9 * eta**2 + 72 * eta + 108
    return (2 * PI**2 * L * r * (term1 + term2 + term3)) / (3 * KAPPA**2)


# T³×S¹ Potentials (Isotropic: R1=R2=R3=r)
def V_T3_FULL(r, V_param, eta, theta):
    if r <= 0:
        return 1e50
    prefactor = 16 * PI**4 * L / (3 * KAPPA**2)
    term1 = V_param**2 * r**3
    term2 = 6 * V_param * eta * KAPPA**2 * theta * r**2
    term3 = 9 * eta**2 * r
    return prefactor * (term1 + term2 + term3)


def V_T3_TT(r, V_param, eta, theta):
    if r <= 0:
        return 1e50
    prefactor = 16 * PI**4 * L / (3 * KAPPA**2)
    term1 = V_param**2 * r**3
    term2 = 12 * V_param * eta * KAPPA**2 * theta * r**2
    term3 = 9 * eta**2 * r
    return prefactor * (term1 + term2 + term3)


def V_T3_REE(r, V_param, eta, theta):
    return V_T3_FULL(r, V_param, eta, theta)


# Nil³×S¹ Potentials
def V_Nil3_FULL(r, V_param, eta, theta):
    if r <= 0:
        return 1e50
    term1 = 4 * V_param**2 * r**2
    term2 = 8 * V_param * KAPPA**2 * theta * (3 * eta + 1) * r
    term3 = 36 * eta**2 - 24 * eta - 9
    return (4 * PI**4 * L * r * (term1 + term2 + term3)) / (3 * KAPPA**2)


def V_Nil3_TT(r, V_param, eta, theta):
    if r <= 0:
        return 1e50
    term1 = 4 * V_param**2 * r**2
    term2 = 48 * V_param * eta * KAPPA**2 * theta * r
    term3 = 36 * eta**2 - 24 * eta - 9
    return (4 * PI**4 * L * r * (term1 + term2 + term3)) / (3 * KAPPA**2)


def V_Nil3_REE(r, V_param, eta, theta):
    if r <= 0:
        return 1e50
    term1 = 4 * V_param**2 * r**2
    term2 = 8 * V_param * KAPPA**2 * theta * (3 * eta - 1) * r
    term3 = 36 * eta**2 - 24 * eta - 9
    return (4 * PI**4 * L * r * (term1 + term2 + term3)) / (3 * KAPPA**2)


# Registry
POTENTIAL_FUNCTIONS = {
    ('S3', 'FULL'): V_S3_FULL,
    ('S3', 'TT'): V_S3_TT,
    ('S3', 'REE'): V_S3_REE,
    ('T3', 'FULL'): V_T3_FULL,
    ('T3', 'TT'): V_T3_TT,
    ('T3', 'REE'): V_T3_REE,
    ('Nil3', 'FULL'): V_Nil3_FULL,
    ('Nil3', 'TT'): V_Nil3_TT,
    ('Nil3', 'REE'): V_Nil3_REE,
}
