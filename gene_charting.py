import pandas as pd
import numpy as np
import sort
import display
#import variant_charting
np.seterr(all="ignore")
pd.options.mode.chained_assignment = None  # default='warn'
cond = [
    ['breast', 'asian_census', 'census'],
    ['hepatocellular_carcinoma', 'asian_census', 'census'],
    ['lung', 'asian_census', 'census'],
    ['large_intestine', 'asian_census', 'census'],
    ['thyroid', 'asian_census', 'census']
]

def frequency_cal():
    num_asia = 0
    for j in range(0,5):
        ratio = dict()
        ratio_census = dict()
        for i in range (1,3):
            df = pd.read_csv(f"raw_output/census_mutation_output/cosmic_{cond[j][0]}_cancer_{cond[j][i]}_mutation.csv")
            df = df[df['Tier']==1]
            list = df['Gene name'].unique() 
            freq_list = [0]*len(list)
            ratio_list = [1.0]*len(list)
            gene_ratio = {
                        'list' : list, 
                        'freq': freq_list,  #how many time this gene is met over all patients
                        'ratio': ratio_list, #ratio of each gene
                        }
            sum_patient = df['id_individual'].nunique() 
            for k in range(0, len(list)):
                gene_ratio['freq'][k] = sort._count_gene_freq(list[k],df)
                gene_ratio['ratio'][k] = gene_ratio['freq'][k]/sum_patient
            if i==1:
                ratio = gene_ratio
            else:
                ratio_census = gene_ratio
        display._chart_create(ratio, ratio_census, cond[j][0], 'g', 'Châu Á', 'Thế giới')
        display.write_to_doc(ratio, ratio_census, cond[j][0])

frequency_cal()
            
            






