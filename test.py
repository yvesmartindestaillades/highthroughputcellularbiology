import os
import json
import numpy as np

path = '/Users/ymdt/src/highthroughputcellularbiology/data'

# for each sample
    # for each construct
        # for each section
            # for "pop_avg"
                # turn "num_of_mutations" into np.histogram("num_of_mutations").tolist()
            

for sample in os.listdir(path):
    print(sample)
    if not sample.endswith('.json'):
        continue
    file = os.path.join(path, sample)
    with open(file, 'r') as f:
        data = json.load(f)
    for construct in data:
        if type(data[construct]) != dict:
            continue
        for section in data[construct]:
            if type(data[construct][section]) != dict:
                continue
            data[construct][section]["pop_avg"]["num_of_mutations"] = np.histogram(
                                                                            data[construct][section]["pop_avg"]["num_of_mutations"], 
                                                                            bins=np.arange(0,170,1),
                                                                            )[0].tolist()
                            
    with open(file, 'w') as f:
        json.dump(data, f)