import pandas as pd
import numpy as np
import sort
import display
import docx
import math

np.seterr(all="ignore")
cond = [
    ['breast', 'asian_census', 'census'],
    ['hepatocellular_carcinoma', 'asian_census', 'census'],
    ['lung', 'asian_census', 'census'],
    ['large_intestine', 'asian_census', 'census'],
    ['thyroid', 'asian_census', 'census']
]

COSMIC_to_civicdb = {
    'molecular_profile' : 'MUTATION_REMARK',
    'disease' : 'DISEASE',
    'therapies' : 'DRUG_COMBINATION',
    'evidence_level' : 'ACTIONABILITY_RANK',
    'citation_id' : 'TRIAL_ID',
    'significance' : 'PRIMARY_OUTCOME_MEASURE',
}

COSMIC_source_conversion = pd.DataFrame(data = ['PubMed','DOI','Unindexed journals, Conference abstracts',
                                                'ClinicalTrials.gov', 'Corporate website', 'Submitted to clinicaltrials.gov'
                                                ,'FDA drug label','FDA drug label'], columns= ['source_name'])

def Civicfb_filter(evidence):
    evidence = evidence[evidence['evidence_type'] == "Predictive"]
    evidence = evidence[evidence['evidence_direction'] == "Supports"]
    evidence = evidence.query("evidence_level in ['A', 'B', 'C']")
    return evidence

def Civicdb_process():
    pd.read_table(f"raw_output/civicdb/civic_aevidence.tsv").to_csv("civic_aevidence.csv", index=False)
    evidence = Civicfb_filter(pd.read_csv(f"civic_aevidence.csv", index_col=False))
    for j in range(0,5):
        data_asian = display._load_data(pd.read_csv(f"output/{cond[j][0]}_cancer/drug/panel_gene_asia.csv"), display.listing_num)
        data_census = display._load_data(pd.read_csv(f"output/{cond[j][0]}_cancer/drug/panel_gene_world.csv"), display.listing_num)
        filter_a = pd.DataFrame().reindex_like(evidence).dropna()
        filter_w = pd.DataFrame().reindex_like(evidence).dropna()
        writer_a = pd.ExcelWriter(f"output/panels/final_{cond[j][0]}_panel_asia.xlsx", engine='xlsxwriter')
        writer_w = pd.ExcelWriter(f"output/panels/final_{cond[j][0]}_panel_world.xlsx", engine='xlsxwriter')
        relevant_entries = evidence[evidence['disease_type'] == cond[j][0]]
        relevant_entries.insert(6, "therapy_rank", 0)
        relevant_entries.insert(len(relevant_entries.columns), column="source", value="Civicdb")
        relevant_entries.drop(columns=['is_flagged'], axis=1, inplace=True)
        tier_1_drugs = pd.read_excel(f"raw_output/VN_approved.xlsx", sheet_name=cond[j][0])
        tier_2_drugs = pd.read_excel(f"raw_output/VN_unapproved.xlsx", sheet_name=cond[j][0])

        for i in range(0, len(data_asian)): 
            ########################
            gene_name = data_asian['list'][i]
            filter_a = relevant_entries[relevant_entries['molecular_profile'].str.contains(gene_name)].reset_index(drop=True)
            for k in range(0, len(tier_1_drugs)):
                filter_a.loc[(filter_a['therapies'].str.contains(tier_1_drugs['drug_name'][k], na=False)) 
                             & (filter_a['evidence_level'].isin(['A', 'B'])), 'therapy_rank'] = 1
            for k in range(0, len(tier_2_drugs)):
                filter_a.loc[(filter_a['therapies'].str.contains(tier_2_drugs['drug_name'][k], na=False)) 
                             & (filter_a['evidence_level'].isin(['A', 'B'])) & (filter_a['therapy_rank']==0), 'therapy_rank'] = 2
            filter_a.loc[(filter_a['evidence_level'].isin(['A', 'B'])) & (filter_a['therapy_rank']==0), 'therapy_rank'] = 3     
            filter_a.loc[(filter_a['evidence_level']=='C'), 'therapy_rank'] = 4
            ########################
            filter_a = filter_a.sort_values(['therapy_rank'], ascending=True, ignore_index=True).reset_index(drop=True)
            filter_a.to_excel(writer_a, sheet_name=gene_name)

        for i in range(0, len(data_census)): 
            ########################
            gene_name = data_census['list'][i]
            filter_w = relevant_entries[relevant_entries['molecular_profile'].str.contains(gene_name)].reset_index(drop=True)
            for k in range(0, len(tier_1_drugs)):
                filter_w.loc[(filter_w['therapies'].str.contains(tier_1_drugs['drug_name'][k], na=False)) 
                             & (filter_w['evidence_level'].isin(['A', 'B'])), 'therapy_rank'] = 1
            for k in range(0, len(tier_2_drugs)):
                filter_w.loc[(filter_w['therapies'].str.contains(tier_2_drugs['drug_name'][k], na=False)) 
                             & (filter_w['evidence_level'].isin(['A', 'B'])) & (filter_w['therapy_rank']==0), 'therapy_rank'] = 2
            filter_w.loc[(filter_w['evidence_level'].isin(['A', 'B'])) & (filter_w['therapy_rank']==0), 'therapy_rank'] = 3     
            filter_w.loc[(filter_w['evidence_level']=='C'), 'therapy_rank'] = 4
            ########################
            filter_w = filter_w.sort_values(['therapy_rank'], ascending=True, ignore_index=True).reset_index(drop=True)
            filter_w.to_excel(writer_w, sheet_name=gene_name)
        writer_a.close() 
        writer_w.close()

def COSMIC_filter(evidence):
    evidence = evidence.query("ACTIONABILITY_RANK in [1, 2, 3, 4]")
    #evidence = evidence[evidence['MUTATION_REMARK'].str.contains("unspecified") == False].reset_index(drop=True)
    return evidence


def _copy_data(civic_ent, cosmic_ent):
    inc = pd.DataFrame(index=cosmic_ent.index,columns=civic_ent.columns)
    civic_ent.drop(civic_ent.columns[0],axis=1,inplace=True)
    #print(cosmic_ent)
    for k in range(0, len(cosmic_ent)):
        inc['therapy_rank'][k] = 0
        inc['source'][k] = "COSMIC"
        inc['source_type'][k] = COSMIC_source_conversion['source_name'][cosmic_ent['SOURCE_TYPE'][k]-1]
        for ent in COSMIC_to_civicdb:
            inc[ent][k] = cosmic_ent[COSMIC_to_civicdb[ent]][k]
    civic_ent = pd.concat([civic_ent.reset_index(drop=True), inc.reset_index(drop=True)]).reset_index(drop=True)
    return civic_ent

def COSMIC_process():
    print(len(set(['asd', 'asd', 'aaa', 'aaa'])))
    gen_total = []
    evidence = COSMIC_filter(pd.read_table("raw_output/cosmic_actionability.tsv"))
    for j in range(0,5):
        relevant_entries = evidence[evidence['disease_type'] == cond[j][0]].reset_index(drop=True)
        filename_a = f"output/panels/final_{cond[j][0]}_panel_asia.xlsx"
        filename_w = f"output/panels/final_{cond[j][0]}_panel_world.xlsx"
        data_asian = display._load_data(pd.read_csv(f"output/{cond[j][0]}_cancer/drug/panel_gene_asia.csv"), display.listing_num)
        data_census = display._load_data(pd.read_csv(f"output/{cond[j][0]}_cancer/drug/panel_gene_world.csv"), display.listing_num)
        tier_1_drugs = pd.read_excel(f"raw_output/VN_approved.xlsx", sheet_name=cond[j][0])
        tier_2_drugs = pd.read_excel(f"raw_output/VN_unapproved.xlsx", sheet_name=cond[j][0])    
        for i in range(0, len(data_asian)): 
            gene_name = data_asian['list'][i]
            filter_a = pd.read_excel(f"output/panels/final_{cond[j][0]}_panel_asia.xlsx", sheet_name=gene_name)
            slice = relevant_entries[relevant_entries['GENE'] == gene_name].reset_index(drop=True)
            filter_a = _copy_data(filter_a, slice)
            ########################            
            for k in range(0, len(tier_1_drugs)):
                filter_a.loc[(filter_a['therapies'].str.contains(tier_1_drugs['drug_name'][k], na=False)) 
                             & (filter_a['evidence_level'].isin([1,2,3])), 'therapy_rank'] = 1
            for k in range(0, len(tier_2_drugs)):
                filter_a.loc[(filter_a['therapies'].str.contains(tier_2_drugs['drug_name'][k], na=False)) 
                             & (filter_a['evidence_level'].isin([1,2,3])) & (filter_a['therapy_rank']==0), 'therapy_rank'] = 2
            filter_a.loc[(filter_a['evidence_level'].isin([1,2,3])) & (filter_a['therapy_rank']==0), 'therapy_rank'] = 3     
            filter_a.loc[(filter_a['evidence_level']==4), 'therapy_rank'] = 4
            ########################
            filter_a = filter_a.sort_values(['therapy_rank'], ascending=True).reset_index(drop=True)            
            with pd.ExcelWriter(filename_a, mode="a", engine="openpyxl", if_sheet_exists='replace') as writer:
                filter_a.to_excel(writer, sheet_name=gene_name)  
        num = 0
        for i in range(0, len(data_census)): 
            gene_name = data_census['list'][i]
            filter_w = pd.read_excel(f"output/panels/final_{cond[j][0]}_panel_world.xlsx", sheet_name=gene_name)
            slice = relevant_entries[relevant_entries['GENE'] == gene_name].reset_index(drop=True)
            filter_w = _copy_data(filter_w, slice)
            ########################            
            for k in range(0, len(tier_1_drugs)):
                filter_w.loc[(filter_w['therapies'].str.contains(tier_1_drugs['drug_name'][k], na=False)) 
                             & (filter_w['evidence_level'].isin([1,2,3])), 'therapy_rank'] = 1
            for k in range(0, len(tier_2_drugs)):
                filter_w.loc[(filter_w['therapies'].str.contains(tier_2_drugs['drug_name'][k], na=False)) 
                             & (filter_w['evidence_level'].isin([1,2,3])) & (filter_w['therapy_rank']==0), 'therapy_rank'] = 2
            filter_w.loc[(filter_w['evidence_level'].isin([1,2,3])) & (filter_w['therapy_rank']==0), 'therapy_rank'] = 3     
            filter_w.loc[(filter_w['evidence_level']==4), 'therapy_rank'] = 4
            ########################
            filter_w = filter_w.sort_values(['therapy_rank'], ascending=True).reset_index(drop=True)      
            if len(filter_w) > 0 : 
                num += 1       
                gen_total.append(gene_name)
            with pd.ExcelWriter(filename_w, mode="a", engine="openpyxl", if_sheet_exists='replace') as writer:
                filter_w.to_excel(writer, sheet_name=gene_name)  
        print(f"{cond[j][0]} : {num} gen")
        
        #if cond[j][0] == 'hepatocellular_carcinoma': print(data_census)
    print(f"Number of unique gene = {len(set(gen_total))}, n = {display.listing_num}")

def COSMIC_res_process():

Civicdb_process()
COSMIC_process()
COSMIC_res_process()
            






