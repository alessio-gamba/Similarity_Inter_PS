# Python 2.7
# The script requires two folders called 'in_data' and 'out_data'

from Omim_tools_1 import *


# 1) build new morbidmap txt from original morbidmap 
morbidmap = Read_morbidmap1("in_data/morbidmap.txt")
D_gene = Read_gene_info("in_data/Homo_sapiens.gene_info")
# print the output file
New_morbidmap(morbidmap, D_gene, "out_data/new_morbidmap.txt")


# 2) build new phenotypic series txt from original phenotypic series
D_morbidmap = Read_morbidmap2("out_data/new_morbidmap.txt")
valid = D_morbidmap.keys()
D_PS = Read_phen_ser("in_data/phenotypicSeries.txt", valid)
D_PS = Add_LH_to_PS(D_morbidmap, D_PS)
# print the output file
New_PS(D_PS, D_morbidmap, "out_data/new_PS.txt")


# 3) build:   MIM  -> HP annotation
#             gene -> GO annotation
D_annotation = Import_mim2hp('in_data/OMIM_ALL_FREQUENCIES_diseases_to_genes_to_phenotypes.txt')
Print_dict(D_annotation, 'out_data/new_mim2hp.txt')
D_annotation = Import_gene2go('in_data/gene2go')
Print_dict(D_annotation, 'out_data/new_gene2go.txt')
