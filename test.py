import os
import json
import numpy as np

path = '/Users/ymdt/src/highthroughputcellularbiology/data'

# for each sample
    # for each construct
        # for each section
            # for "pop_avg"
                # turn "num_of_mutations" into np.histogram("num_of_mutations").tolist()
            


print([f for f in os.listdir(path) if 'degree' in f])