import os
import csv
import time

from sage.all import *

from enumerate_corner_polyhedra import all_abelian_groups, name_of_group, group_iter
from shooting_theorem import shoot_n_times, remove_duplicates, count_faces


N = 1000
faces_dir = 'shooting_faces'
if not os.path.isdir(faces_dir):
    os.mkdir(faces_dir)

def main():
    timing_results = open(faces_dir + os.sep + 'timing_results.csv', 'w')
    time_writer = csv.writer(timing_results)
    time_writer.writerow(['Group', '$g_0$', 'Computation Time', '\# Unique Faces', 'Proportion'])
    for i in range(3, 31):
        for group in all_abelian_groups(i):
            group_name = name_of_group(group)
            g0 = group[0]
            filename = 'face_counts_{}_0.txt'.format(group_name)

            gens = group.gens()
            orders = [gen.order() for gen in gens]

            start_time = time.time()
            faces = shoot_n_times(group, g0, N, True)
            print(group_name, 'done_shooting')
            tot_time = time.time() - start_time
            uniqs = remove_duplicates(faces)
            num_uniqs = len(uniqs)
            time_writer.writerow([group_name, '0', tot_time, num_uniqs, num_uniqs / N])
            timing_results.flush()
            face_counts = count_faces(faces, uniqs)
            with open(faces_dir + os.sep + filename, 'w') as f:
                group_writer = csv.writer(f)
                f.write('count,b,' + ','.join('a_{}'.format(i) for i in range(i-1)))
                f.write('\n')
                for face, count in face_counts.items():
                    group_writer.writerow([count] + ['{:.2f}'.format(x) for x in face])

if __name__ == '__main__':
    main()
