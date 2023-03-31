from config import *
from util import *
import generate_dataset
import dreem

if not os.path.exists(saved_feather):
    print('Loading JSONs...')
    data = []
    for f in os.listdir(path_data):
        if f.endswith('.json'):
            print(' ' + f, end=' ')
            data.append(json.load(open(path_data + f, 'r')))
    print('Done reading data.')
    
    print('Creating study...')
    study = dreem.draw.study.Study(
        data = data
    )
    print('Done creating study.')
    
    print('Filtering study...')
    study.df = study.df[study.df['min_cov'] > min_base_coverage].reset_index(drop=True)
    # only keep the references that have 8 sections 
    study.df = study.df[study.df['reference'].isin(study.df.groupby(['sample','reference']).filter(lambda x: len(x) == 8)['reference'])].reset_index(drop=True)
    print('Finding frame shift ROI...')
    study.df= generate_dataset.find_frame_shift_ROI(study)
    print('Done finding frame shift ROI.')
    
    print('Saving study to df.feather...')
    study.df.to_feather(saved_feather)
    print('Done saving study to df.feather.')
else:
    study = dreem.draw.study.Study()
    print('Reading study from df.feather...')
    study.df = pd.read_feather(saved_feather)
    print('Done reading study from df.feather.')

study.df['family'] = study.df['reference'].apply(lambda x: x.split('=')[1].split('-')[0]) 
