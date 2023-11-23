'''
pip3 installed instrain from github (31/10/23) into 20230313_scripts_env
###! Successfully installed asteval-0.9.31 biopython-1.74 future-0.18.3 h5py-3.10.0 inStrain-1.8.0 iniconfig-2.0.0 lmfit-1.2.2 pluggy-1.3.0 pysam-0.22.0 pytest-7.4.3 uncertainties-3.1.7
'''

import inStrain
import inStrain.SNVprofile
import pandas as pd
from collections import defaultdict 

IS = inStrain.SNVprofile.SNVprofile('results/10_instrain/02_instrain_profile/A2-2')
genes_SNP_count_df = IS.get('genes_SNP_count')
genes_SNP_count_df.columns

genes_SNP_count_df = genes_SNP_count_df[genes_SNP_count_df['mm'] <= 5]
genes_SNP_count_df = genes_SNP_count_df.sort_values('mm')\
            .drop_duplicates(subset=['gene'], keep='last')\
            .sort_index().drop(columns=['mm'])
