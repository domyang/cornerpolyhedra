import time
import itertools
import os

import numpy as np
from pyomo.opt import SolverFactory
from fractions import Fraction, gcd
from functools import reduce

from ReferenceModel import model
from enumerate_corner_polyhedra import gomory_system
from remove_redundancies import remove_redundancies

opt = SolverFactory('gurobi')

def shoot_n_times(group, g0, n, remove_redund=True):
    eqs, ineqs = gomory_system(group, g0)
    if remove_redund:
        eqs, ineqs = [], remove_redundancies(eqs, ineqs)
    faces = []
    order = group.order()
    for _ in range(n):
        v = uniform_nonnegative_nsphere(order-1)
        faces.append(solve_lp(v, eqs, ineqs))
    return faces


def uniform_nonnegative_nsphere(n):
    xs = np.random.randn(n)
    l = np.sqrt(sum(xs**2))
    return abs(xs) / l


def construct_params_file(filename, v, eqs, ineqs):
    n = len(v)
    m = len(ineqs)
    p = len(eqs)

    with open(filename, 'w') as f:
        f.write('param n := {} ;\n'.format(n))
        f.write('param m := {} ;\n'.format(m))
        f.write('param p := {} ;\n\n'.format(p))

        f.write('param v := \n')
        for i, vi in enumerate(v):
            f.write('{} {}\n'.format(i+1, vi))
        f.write(';\n\n')

        if len(eqs) > 0:
            f.write('param AEq :=\n')
            for i, eq in enumerate(eqs):
                for j, var in enumerate(eq[1:]):
                    f.write('{} {} {}\n'.format(i+1, j+1, var))
            f.write(';\n\n')

            f.write('param bEq := \n')
            for i, eq in enumerate(eqs):
                f.write('{} {}\n'.format(i+1, -eq[0]))
            f.write(';\n');

        if len(ineqs) > 0:
            f.write('param AIneq :=\n')
            for i, ineq in enumerate(ineqs):
                for j, var in enumerate(ineq[1:]):
                    f.write('{} {} {}\n'.format(i+1, j+1, var))
            f.write(';\n\n')

            f.write('param bIneq := \n')
            for i, ineq in enumerate(ineqs):
                f.write('{} {}\n'.format(i+1, -ineq[0]))
            f.write(';\n');


def lcm(a, b):
    return a * b // gcd(a, b)


def common_integer(numbers):
    """
    This will accept a list of floating points and return the list of integers
    that is a multiple of the floating points 
    Taken from https://stackoverflow.com/questions/36551393/how-to-get-a-common-integer-in-python-from-a-float-list
    """
    fractions = [Fraction(str(n)).denominator for n in numbers]
    multiple = reduce(lcm, [f.denominator for f in fractions])
    ints = [f * multiple for f in fractions]
    divisor = reduce(gcd, ints)
    return [int(n / divisor) for n in ints]


count = 0
def solve_lp(v, eqs, ineqs):
    filename = 'model_data_{}.dat'.format(count)
    construct_params_file(filename, v, eqs, ineqs)
    instance = model.create_instance(filename)
    results = opt.solve(instance)
    face = [-1] + [instance.pi[i+1]() for i in range(len(v))]
    mult = min(face[1:])
    face = [x*(1/mult) for x in face]
    os.remove(filename)
    return face


def count_faces(faces, uniqs):
    face_counts = {tuple(uniq): 0 for uniq in uniqs}
    for face in faces:
        for other in face_counts:
            if np.allclose(face, other):
                face_counts[other] += 1
    return face_counts

def remove_duplicates(list_of_lists):
    list_of_lists.sort()
    uniqs = [list_of_lists[0]]
    for face in list_of_lists[1:]:
        if np.allclose(face, uniqs[-1]):
            continue
        else:
            uniqs.append(face)
    return uniqs
