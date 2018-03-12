import os

import numpy as np

def remove_redundancies(eqs, ieqs):
    file_num = os.getpid()
    file1 = 'tmp_{}.txt'.format(file_num)
    file2 = 'tmp_{}_.txt'.format(file_num)
    write_lrslib_format_hrep(eqs, ieqs, str(file_num), file1)
    run_redund(file1, file2)
    new_ieqs = read_lrslib_format_hrep(file2)

    new_eqs = []
    to_delete = []

    for i, ieq in enumerate(new_ieqs):
        for j, other in enumerate(new_ieqs[i+1:]):
            if np.allclose(np.array(ieq), -np.array(other)):
                new_eqs.append(ieq)
                to_delete.extend([i, i+j+1])

    for i in to_delete[::-1]:
        del new_ieqs[i]

    os.remove(file1)
    os.remove(file2)
    return new_eqs, new_ieqs

def write_lrslib_format_hrep(eqs, ieqs, name, filename):
    with open(filename, 'w') as f:
        f.write(name + '\n')
        f.write('H-representation\n')
        f.write('begin\n')
        m = len(ieqs) + len(eqs) * 2
        n = len(eqs[0])
        f.write('{} {} integer\n'.format(m, n))
        for ineq in ieqs:
             f.write(' '.join(map(str, ineq)) + '\n')

        for eq in eqs:
            f.write(' '.join(map(str, eq)) + '\n')
            f.write(' '.join(map(str, [-1*var for var in eq])) + '\n')
        f.write('end\n')


def run_redund(file1, file2):
    os.system('redund {} > {}'.format(file1, file2))


def read_lrslib_format_hrep(filename):
    with open(filename) as f:
        for line in f:
            if line.startswith('begin'):
                break

        ieqs = []

        for line in f:
            if line.startswith('*'):
                continue
            if line.startswith('end'):
                break

            try:
                ieqs.append(map(int, line.rstrip().split()))
            except:
                continue

    return ieqs
