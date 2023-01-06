import pandas as pd
import numpy as np
from tqdm import tqdm
from util import *
from config import *

# Load df and add family and flank
df = fasta_to_df('../data/reference.fasta')
df['family'] = df['construct'].apply(lambda s: s.split('_')[1].split('-')[0].split('=')[1])
df['flank'] = df['construct'].apply(lambda s: s.split('-')[2].split('=')[0])
        
# add section and name_section
df_out = pd.DataFrame()
for f in boundary.keys():
    df_temp = df.copy()
    df_temp['section_start'],  df_temp['section_end'] = df_temp['sequence'].apply(lambda x: boundary[f](x)[0]+1), df_temp['sequence'].apply(lambda x: boundary[f](x)[1])
    df_temp['section'] = df_temp.apply(lambda x: str(x['section_start']) + '-' + str(x['section_end']), axis=1)
    df_temp['name_section'] = f
    df_temp['barcode'] = df_temp['sequence'].apply(lambda x: x[boundary['barcode'](x)[0]:boundary['barcode'](x)[1]])
    df_temp['barcode_start'] = df_temp['sequence'].apply(lambda x: boundary['barcode'](x)[0]+1)
    df_temp['barcode_end'] = df_temp['sequence'].apply(lambda x: boundary['barcode'](x)[1])
    df_temp['secondary_signature'] = df_temp['sequence'].apply(lambda x: x[boundary['ROI'](x)[0]:boundary['ROI'](x)[1]])
    df_temp['secondary_signature_start'] = df_temp['sequence'].apply(lambda x: boundary['ROI'](x)[0]+1)
    df_temp['secondary_signature_end'] = df_temp['sequence'].apply(lambda x: boundary['ROI'](x)[1])
    df_out = pd.concat([df_out, df_temp], axis=0)

# drop sequence, reset index, save to csv
df = df_out.reset_index(drop=True)

df.to_csv('../data/library.csv', index=False)