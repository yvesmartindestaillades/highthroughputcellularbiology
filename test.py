from dreem.draw.study import Study
import json 

# load json
data = [
    json.load(open('/Users/ymdt/src/highthroughputcellularbiology/data/37degrees_1percent_1_S14_L001.json','r')),
    json.load(open('/Users/ymdt/src/highthroughputcellularbiology/data/45degrees_1_S18_L001.json','r')),
]
study = Study(
    data = data,
)