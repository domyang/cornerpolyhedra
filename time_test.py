from sage.all import *
from enumerate_corner_polyhedra import gomory_system, poly_faces, all_abelian_groups, name_of_group
import time
import signal


if __name__ == '__main__':
    with open('more.csv', 'w') as f:
        for i in range(2, 30):
            for group in all_abelian_groups(i):
                name = name_of_group(group)
                eqs, ieqs = gomory_system(group, group[0])

                corner_poly = Polyhedron(eqns=eqs, ieqs=ieqs, backend='normaliz')
                num_verts = corner_poly.n_vertices()
                num_faces = corner_poly.n_facets()

                string = ','.join(map(str, [name, num_faces, num_verts]))
                f.write(string)
                f.write('\n')
