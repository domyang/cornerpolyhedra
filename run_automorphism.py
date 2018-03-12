import csv
import time

from sage.all import *

import automorphisms as morph
from enumerate_corner_polyhedra import master_polyhedron, name_of_group

def main():
    timing_results = open('automorphism_results.csv', 'w')
    writer = csv.writer(timing_results)
    writer.writerow(['Group', '\# Automorphisms', 'Time Using Automorphisms', 'Time Without'])
    for i in range(12, 20):
        group = CyclicPermutationGroup(i)
        classes, mappings = morph.automorphism_orbits(i)
        start_time = time.time()
        for g in group:
            poly = master_polyhedron(group, g, 'normaliz')
        tot_time = time.time() - start_time

        start_time = time.time()
        for mapping in mappings:
            poly = master_polyhedron(group, list(group)[mapping], 'normaliz')
        tot_time2 = time.time() - start_time

        writer.writerow([name_of_group(group), len(classes), tot_time2, tot_time])


if __name__ == '__main__':
    main()
