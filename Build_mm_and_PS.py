from Omim_functions_1 import *


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


# Phenotypes :
# ------------
#
# Each Phenotype is followed by its MIM number, if different from that
# of the locus/gene, and then followed by its phenotype mapping
# key in parentheses (explanation below).
#
#
# Phenotype Mapping key - Appears in parentheses after a disorder :
# -----------------------------------------------------------------
#
# (1) The disorder is placed on the map based on its association with
# a gene, but the underlying defect is not known.
# (2) The disorder has been placed on the map by linkage or other
# statistical method; no mutation has been found.
# (3) The molecular basis for the disorder is known; a mutation has been
# found in the gene.
# (4) A contiguous gene deletion or duplication syndrome, multiple genes
# are deleted or duplicated causing the phenotype.
