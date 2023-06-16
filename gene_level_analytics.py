import csv
import pandas as pd

# with open("gene_specific_summary.tsv") as recent_summary:
#     rec_sum_row = csv.reader(recent_summary, delimiter="\t")
#     for gene in rec_sum_row:

rec_sum = pd.read_table("/Users/justin/Documents/repos/variant-reassessment/gene_specific_summary.tsv")
old_sum = pd.read_table("/Users/justin/Documents/repos/variant-reassessment/gene_specific_summary_2018-01.tsv")
# Symbol, GeneID, Total_submissions, Total_alleles, Submissions_reporting_this_gene, Alleles_reported_Pathogenic_Likely_pathogenic, Gene_MIM_number, Number_uncertain, Number_with_conflicts
print(old_sum)

# Remove all rows that are not present in the old summary
rec_sum = rec_sum[rec_sum.GeneID.isin(old_sum.GeneID)]

rec_sum = rec_sum.set_index('GeneID')
old_sum = old_sum.set_index('GeneID')
# rec_sum = rec_sum.replace('-', 0)
# old_sum = old_sum.replace('-', 0)

rec_sum['Pathogenic_change'] = rec_sum['Alleles_reported_Pathogenic_Likely_pathogenic'] - old_sum['Alleles_reported_Pathogenic_Likely_pathogenic']

rec_sum.to_csv("/Users/justin/Documents/repos/variant-reassessment/Updated_summary.tsv", sep="\t")