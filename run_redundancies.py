import time
import csv

from sage.all import *

from enumerate_corner_polyhedra import gomory_system, all_abelian_groups, poly_faces, name_of_group
from remove_redundancies import remove_redundancies

def main():
    with open('redundancy_results.csv', 'w') as f:
        writer = csv.writer(f)
        writer.writerow(['Group', '$g_0$', '\# Constraints', '\# Reduced Constraints', 'Computation Time', 'Reduced Time'])
        for i in range(2, 20):
            for group in all_abelian_groups(i):
                eqs, ieqs = gomory_system(group, group[0])

                red_eqs, red_ieqs = remove_redundancies(eqs, ieqs)

                start_time = time.time()
                polyhedron = Polyhedron(eqns=eqs, ieqs=ieqs, backend='normaliz')
                facets = poly_faces(polyhedron)
                corner_polyhedron = Polyhedron(ieqs=facets, backend='normaliz')
                time1 = time.time() - start_time

                start_time = time.time()
                polyhedron = Polyhedron(eqns=red_eqs, ieqs=red_ieqs, backend='normaliz')
                facets = poly_faces(polyhedron)
                corner_polyhedron = Polyhedron(ieqs=facets, backend='normaliz')
                time2 = time.time() - start_time

                writer.writerow([name_of_group(group), 0, len(eqs) + len(ieqs), len(red_eqs) + len(red_ieqs), time1, time2])

if __name__ == '__main__':
    main()
