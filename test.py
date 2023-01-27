import os
import json
import dreem 

# List of JSON files path
json_files = [os.path.join('data', f) for f in os.listdir('data') if f.endswith('S1.json')]

# Read JSON files
data = [json.load(open(json_file, 'r')) for json_file in json_files]

# Create study
study = dreem.draw.study.Study(
    data = data
)

# Sanity check
sample, construct, section = study.df.iloc[0][['sample', 'construct', 'section']]
plot = study.mutation_fraction(
    sample = sample,
    construct = construct,
    section = section
)

