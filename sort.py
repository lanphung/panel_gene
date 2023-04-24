def _sort_by_gene(df): #sort by Gene name and ID_tumour to get the actual number of gene occurances
    df['gene_freq'] = df.groupby('Gene name', sort=False)['Gene name'].transform('count')
    df['id_tumour_freq'] = df.groupby('ID_tumour', sort=False)['ID_tumour'].transform('count')
    return df.sort_values(['gene_freq', 'id_tumour_freq'], ascending=[False, False])

def _sort_by_variant(df): #sort by Gene name and ID_tumour to get the actual number of gene occurances
    df['mutation_aa_freq'] = df.groupby('Mutation AA', sort=False)['Mutation AA'].transform('count')
    df['gene_freq'] = df.groupby('Gene name', sort=False)['Gene name'].transform('count')
    df['id_tumour_freq'] = df.groupby('ID_tumour', sort=False)['ID_tumour'].transform('count')
    return df.sort_values(['gene_freq', 'mutation_aa_freq', 'id_tumour_freq', ], ascending=[False, False, False])

def _count_gene_freq(gene_name, df):
    dg = df[df['Gene name'] == gene_name]
    return dg['id_individual'].nunique()
