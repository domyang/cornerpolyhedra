import time
import os
import itertools
import csv
import signal
from multiprocessing import Pool

from sage.all import *

from remove_redundancies import remove_redundancies

backends = ['ppl']

def gomory_system(group, g0):
    """
    Creates the system of equations and inequalities which define the
    gomory system for a given group. Equations are stored in the form
    b + a1x1 + ... + anxn = 0 as a list of the n+1 coefficients.
    Inequalities are stored as
    b + a1x1 + ... + anxn >= 0 in a list of the n+1 coefficients.

    The indices of the arrays will correspond to the indices of the group
    passed in, e.g., if the group is ordered [e, x1, x2, ..., xn] where e is
    the identity, the coefficient of x1 will have the corresponding coef.

    Pulled from Theorem 18 of classic paper:
        Some Combinatorial Problems in Optimization

    Returns: Equations, Inequalities
    """
    n = group.order()
    g_plus = list(group)[1:]
    g_indices = {g:i+1 for i, g in enumerate(g_plus)}

    eqs = []

    if not g0.is_one():
        eq0 = [0]*n
        eq0[g_indices[g0]] += 1
        eq0[0] += -1
        eqs.append(eq0)

    for g in g_plus:
        if g == g0:
            continue
        eq = [0]*n
        eq[g_indices[g]] += 1
        eq[g_indices[g0 * (g.inverse())]] += 1
        eq[0] += -1
        eqs.append(eq)

    ieqs = []

    for g1 in g_plus:
        for g2 in g_plus:
            ieq = [0]*n
            ieq[g_indices[g1]] += 1
            ieq[g_indices[g2]] += 1
            if (g1*g2).is_one():
                continue
            ieq[g_indices[g1*g2]] += -1
            ieqs.append(ieq)

    for g in g_plus:
        ieq = [0]*n
        ieq[g_indices[g]] += 1
        ieqs.append(ieq)

    return eqs, ieqs

def poly_faces(poly):
    """
    This will return a list of a lists, each with n+1 components.
    This will return the matrix A and vector b in this form.
    Ax + b >= 0 .
    This will be stored as [b A]
    """
    facets = []
    for vertex in poly.vertices():
        # We first make all points integers
        point = vertex.vector()
        denom = point.denominator()
        point *= denom
        facets.append([-denom] + list(point))
    n = len(facets[0]) - 1
    for i in range(n):
        row = [0]*(n+1)
        row[i+1] = 1
        facets.append(row)

    return facets


def master_polyhedron(group, g0, backend, remove_redund=True):
    eqs, ieqs = gomory_system(group, g0)

    if remove_redund:
        eqs, ieqs = remove_redundancies(eqs, ieqs)

    polyhedron = Polyhedron(eqns=eqs, ieqs=ieqs, backend=backend)
    facets = poly_faces(polyhedron)

    corner_polyhedron = Polyhedron(ieqs=facets, backend=backend)
    return corner_polyhedron


def faces_verts_ind(poly):
    poly_faces = [face.vector() for face in poly.Hrepresentation()]

    verts = poly.Vrepresentation()
    vert_indices = [i for i, v in enumerate(verts) if v.is_vertex()]
    poly_verts = [verts[i].vector() for i in vert_indices]

    incidence_mat = poly.incidence_matrix()[vert_indices]

    return poly_faces, poly_verts, incidence_mat

def write_vectors(filename, vectors, faces):
    with open(filename, 'w') as f:
        if faces:
            d = len(vectors[0]) - 1
            eqn = ' + '.join('a_{}x_{}'.format(i, i) for i in range(d))
            f.write('# In order b + ' + eqn + ' >= 0\n')
            f.write('index,b,' + ','.join('a_{}'.format(i) for i in range(d)))
            f.write('\n')
        if not faces:
            d = len(vectors[0])
            f.write("Index," + ','.join('t_{}'.format(i) for i in range(d)))
            f.write('\n')
        for i, v in enumerate(vectors):
            f.write(str(i) + ',')
            f.write(','.join(map(str, v)))
            f.write('\n')

def write_matrix(filename, matrix):
    with open(filename, 'w') as f:
        f.write('# This encodes the incidence matrix, the first row')
        f.write(' lists the indices of the faces. The first column')
        f.write(' lists the indices of the vertices.\n')
        f.write('# Each element of the matrix is either 1 if a given vertex')
        f.write(' is adjacent to a face or 0 if not.\n')

        num_verts, num_faces = matrix.dimensions()
        f.write(',' + ','.join(map(str, range(num_faces))) + '\n')
        for i, row in enumerate(matrix):
            f.write(str(i) + ',')
            f.write(','.join(map(str, row)))
            f.write('\n')

def all_abelian_groups(n):
    """
    Return all abelian groups of order n
    """
    factorization = dict(factor(n))
    partitions = {}
    cyclic_groups = []
    for num, exponent in factorization.items():
        partitions = Partitions(exponent).list()
        prime_power_groups = []
        for part in partitions:
            cyclics = []
            for i in part:
                cyclics.append(CyclicPermutationGroup(num**i))
            prime_power_groups.append(direct_product_permgroups(cyclics))
        cyclic_groups.append(prime_power_groups)

    if len(cyclic_groups) == 1:
        return cyclic_groups[0]

    all_groups = []
    for groups in itertools.product(*cyclic_groups):
        all_groups.append(direct_product_permgroups(groups))

    return all_groups


def name_of_group(group):
    gens = group.group_generators()
    string = ''.join('c' + str(gen.order()) for gen in gens)
    return string


def group_iter(group):
    gens = group.gens()
    orders = [gen.order() for gen in gens]
    elems = {}
    for prod in itertools.product(*[range(order) for order in orders]):
        comps = []
        for gen, x in zip(gens, prod):
            comps.append(gen ** x)
        g = group[0]
        for comp in comps:
            g *= comp
        elems[prod] = g
    return elems

class Timeout:
    def __init__(self, seconds=1, error_message='Timeout'):
        self.seconds = seconds
        self.error_message = error_message
    def handle_timeout(self, signum, frame):
        raise RuntimeError(self.error_message)
    def __enter__(self):
        signal.signal(signal.SIGALRM, self.handle_timeout)
        signal.alarm(self.seconds)
    def __exit__(self, type, value, traceback):
        signal.alarm(0)

def run_test(group, g0, backend, group_name, g0_name, n):
    try:
        with Timeout(60*60):
            start_time = time.time()
            poly = master_polyhedron(group, g0, backend)
            tot_time = str(time.time() - start_time)
            faces, verts, matrix = faces_verts_ind(poly)
            write_vectors('polyhedra' + os.sep + 'group_{}_{}_{}_faces.csv'.format(n, group_name, g0_name), faces, True)
            write_vectors('polyhedra' + os.sep + 'group_{}_{}_{}_verts.csv'.format(n, group_name, g0_name), verts, False)
            write_matrix('polyhedra' + os.sep + 'group_{}_{}_{}_matrix.csv'.format(n, group_name, g0_name), matrix)
            fp = open('ppl_times.csv', 'a')
            writer = csv.writer(fp)
            writer.writerow([group_name, g0_name, backend, tot_time, len(faces), len(verts)])
    except RuntimeError:
        tot_time = 'N\A'




if __name__ == '__main__':
    poly_dir = 'polyhedra'

    fp = open('ppl_times.csv', 'a')
    writer = csv.writer(fp)
    writer.writerow(['Group', '$g_0$', 'backend', 'Computation Time', '\# Faces', '\# Vertices'])

    arg_list = []
    for n in range(2, 31):
        for group in all_abelian_groups(n):
            group_name = name_of_group(group)
            for decomp, g0 in group_iter(group).items():
                if not g0.is_one(): continue
                g0_name = ','.join(map(str, decomp))
                for backend in backends:
                    run_test(group, g0, backend, group_name, g0_name, n)
