from config import *
from util import *
import generate_dataset

if not os.path.exists(saved_feather):
    study = dreem.draw.study.Study(
        data = [json.load(open(path_data + f, 'r')) for f in os.listdir(path_data) if f.endswith('.json')]
    )
    study.df['deltaG'] = study.df['deltaG'].apply(lambda x: 0 if x == 'void' else x)
    study.df['frame_shift_ROI'] = generate_dataset.find_frame_shift_ROI(study)
    study.df.to_feather(saved_feather)
else:
    study = dreem.draw.study.Study()
    study.df = pd.read_feather(saved_feather)

study.df = study.df[study.df['worst_cov_bases'] > min_base_coverage].reset_index(drop=True)

