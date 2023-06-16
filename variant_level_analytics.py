import csv
import pandas as pd
import matplotlib.pyplot as plt

year1 = '2019'
month1 = '01'
year2 = '2023'
month2 = '03'


def getDataFrame(filepath, date):
    df = pd.read_table(filepath)
    df = df[['Type', 'Name', 'GeneSymbol', 'ClinicalSignificance', 'LastEvaluated', 'Assembly', 'ReviewStatus', 'VariationID']]

    df = df[df['Assembly'] == 'GRCh37']
    df.drop(columns=['Assembly'])
    df_filtered = df[(df['ReviewStatus'] == 'reviewed by expert panel') | (df['ReviewStatus'] == 'practice guideline')]
    df_filtered.drop(columns = ['ReviewStatus'])

    df_filtered = df_filtered[(df_filtered['ClinicalSignificance'].str.contains('Pathogenic|Likely pathogenic|Uncertain significance|Likely benign|Benign', case=False))]

    df_filtered['LastEvaluated'] = df_filtered['LastEvaluated'].replace('-', date)

    df_filtered['LastEvaluated'] = pd.to_datetime(df_filtered['LastEvaluated'])

    # print(df_filtered)

    return df_filtered

def getGeneLevelSummary(old_df, new_df):
    old_gene_series = old_df.groupby(['GeneSymbol', 'ClinicalSignificance']).size()
    new_gene_series = new_df.groupby(['GeneSymbol', 'ClinicalSignificance']).size()

def get_transition_stats(current_summary, old_summary):
    # Get variationIDs that were pathogenic:
    curr_path_variants = set(current_summary.index[(current_summary['ClinicalSignificance'] == 'Pathogenic') | 
        (current_summary['ClinicalSignificance'] == 'Likely pathogenic')])

    curr_benign_variants = set(current_summary.index[(current_summary['ClinicalSignificance'] == 'Benign') | 
        (current_summary['ClinicalSignificance'] == 'Likely benign')])

    old_path_variants = set(old_summary.index[(old_summary['ClinicalSignificance'] == 'Pathogenic') | 
        (old_summary['ClinicalSignificance'] == 'Likely pathogenic')])

    old_benign_variants = set(old_summary.index[(old_summary['ClinicalSignificance'] == 'Benign') | 
        (old_summary['ClinicalSignificance'] == 'Likely benign')])

    if len(old_benign_variants.intersection(curr_path_variants)) > 0:
        print('benign to path: ')
        print(old_benign_variants.intersection(curr_path_variants))

    num_benign_to_path = len(old_benign_variants.intersection(curr_path_variants))
    num_path_to_benign = len(old_path_variants.intersection(curr_benign_variants))

    if len(curr_benign_variants - old_benign_variants) > 0:
        print('new benign: ')
        print(curr_benign_variants - old_benign_variants)

    if len(curr_path_variants - old_path_variants) > 0:
        print('new path: ')
        print(curr_path_variants - old_path_variants)

    path_difference = len(curr_path_variants) - len(old_path_variants)
    benign_difference = len(curr_benign_variants) - len(old_benign_variants)

    return num_benign_to_path, num_path_to_benign, path_difference, benign_difference

current_summary = getDataFrame('/Users/justin/Documents/repos/variant-reassessment/variant_summary_'+ year2 + '-'+ month2 +'.txt', '1-Jan-'+year2)
old_summary = getDataFrame('/Users/justin/Documents/repos/variant-reassessment/variant_summary_'+ year1 + '-'+ month1 + '.txt', '1-Mar-'+year1)

print('Number of old variants:')
print(len(old_summary.index))

print('Number of current variants:')
print(len(current_summary.index))

# # print("Uncertain significances:")
# # print(current_summary[current_summary['ClinicalSignificance'] == 'Uncertain significance'])

# # current_summary.set_index('VariationID', inplace=True)
# # old_summary.set_index('VariationID', inplace=True)

# all_variants = pd.concat([current_summary, old_summary])

# print('All variants: ')
# print(len(all_variants.index))

# new_variants = all_variants.drop_duplicates(subset='VariationID', keep=False, ignore_index=True)

# old_variants = all_variants.drop_duplicates(subset='VariationID', keep='first', ignore_index=True)
# old_variants = pd.concat([old_variants, new_variants])
# old_variants = old_variants.drop_duplicates(subset='VariationID', keep=False, ignore_index=True)

# print('Number of old variants:')
# print(len(old_variants.index))

# print('Number of new variants:')
# print(len(new_variants.index))

# #Create dataframe with old variants from current file.
# old_variants = old_variants.set_index('VariationID')
# current_summary = current_summary.set_index('VariationID')



# Remove all rows that are not present in the old summary
# current_summary = current_summary[current_summary.VariationID.isin(old_summary.VariationID)]

print(current_summary.head(5))

current_summary = current_summary.astype({'VariationID':'int'})

print(current_summary.head(5))

current_summary = current_summary.set_index('VariationID')
old_summary = old_summary.set_index('VariationID')

gene_list = old_summary['GeneSymbol'].unique()
# gene_list = [gene for line in gene_list for gene in line.split()]
curr_gene_list = current_summary['GeneSymbol'].unique()
old_grouped = old_summary.groupby(old_summary.GeneSymbol)
curr_grouped = current_summary.groupby(current_summary.GeneSymbol)

gene_transitions = {}
gene_path_differences = {}
gene_benign_differences = {}

print("Old Gene list length:")
print(len(gene_list))

print("Old Gene list:")
print(gene_list)

print("Current Gene list length:")
print(len(curr_gene_list))

print("Current Gene list:")
print(curr_gene_list)

for gene in gene_list:
    old_vars = old_grouped.get_group(gene)
    curr_vars = curr_grouped.get_group(gene)

    num_benign_to_path, num_path_to_benign, path_difference, benign_difference = get_transition_stats(curr_vars, old_vars)

    gene_transitions[gene] = num_benign_to_path, num_path_to_benign
    gene_path_differences[gene] = path_difference
    gene_benign_differences[gene] = benign_difference

print('Gene transitions:')
print(gene_transitions)
print('Gene pathogenic differences:')
print(gene_path_differences)
print('Gene benign differences:')
print(gene_benign_differences)

"""
plt.figure(1)
# Initialise the subplot function using number of rows and columns
# figure, axis = plt.subplots(2)
# plt.hist(gene_path_differences.values())
# plt[0].xlabel('Number of new pathogenic variants per gene')
# plt[0].ylabel('Frequency')
# plt[0].title('The frequency of new pathogenic 3 or 4 star variants per gene from 2019 to 2023')
# plt[1].hist(gene_benign_differences.values())
# plt[1].xlabel('Number of new benign variants per gene')
# plt[1].ylabel('Frequency')
# plt[1].title('The frequency of new genign 3 or 4 star variants per gene from 2019 to 2023')
# plt.show()

###
plt.hist(gene_path_differences.values())
plt.xlabel('Number of new pathogenic variants per gene')
plt.ylabel('Frequency')
plt.title('The frequency of new pathogenic 3 or 4 star variants per gene from 2019 to 2023')
plt.show()

plt.figure(2)
plt.hist(gene_benign_differences.values())
plt.xlabel('Number of new benign variants per gene')
plt.ylabel('Frequency')
plt.title('The frequency of new benign 3 or 4 star variants per gene from 2019 to 2023')
plt.show()

gene_path_differences.pop('PAH')

plt.figure(3)
plt.hist(gene_path_differences.values())
plt.xlabel('Number of new pathogenic variants per gene')
plt.ylabel('Frequency')
plt.title('The frequency of new pathogenic 3 or 4 star variants per gene from 2019 to 2023 excluding PAH')
plt.show()
"""



# Merge all variants into one dataframe, have column for appeared in 2019 but not in 2023

# Genes that Bob found:
# AKT3
# BRAF
# BRCA1
# BRCA2
# CCDC40;GAA
# CDH1
# CDH23
# CDKL5
# CFTR
# COCH
# ETHE1
# FANCI;POLG
# FOXG1
# GAA
# GJB2
# HRAS
# HTT
# ITGA2B
# ITGB3
# KCNQ4
# KLLN;PTEN
# KRAS
# LDLR
# MAP2K1
# MAP2K2
# MECP2
# MLH1
# MSH2
# MSH6
# MT-ATP6
# MTOR
# MYH7
# MYO6
# MYO7A
# PAH
# PDHA1
# PIK3CA
# PIK3R2
# PMS2
# POLG
# PTEN
# PTPN11
# RAF1
# RUNX1
# RYR1
# SHOC2
# SLC26A4
# SLC9A6
# SOS1
# TCF4
# TECTA
# TP53
# UBE3A
# USH2A