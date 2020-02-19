# Similarity Inter Phenotypic Series

The datasets reported in the enclosed files have been generated and analyzed, as described in the manuscript "The Similarity of the Inherited Diseases (II): Clinical and Biological Similarity between the Phenotypic Series" by Gamba A, Salmona M, Cantù L and Bazzoni G (BMC Medical Genomics, 2020, submitted)

File 01_PS_id_names.txt is the list of the OMIM-derived Phenotypic Series (PS) that have been analyzed in this study. The list provides id and name of each PS.

File 02_similarity_coefficients_CSN.txt reports the clinical similarity coefficient (Sim_HPO) for each pair of PS (PSi_PSj), based on the shared disease phenotypes (i.e., the annotations of the diseases in Human Phenotype Ontology). The PS are the nodes of the Clinical Similarity Network (CSN) reported in Figure 2 of the manuscript, while the similarity coefficients represent the weight w of the edge connecting two PS in the CSN.

File 03_similarity_coefficients_BSNs.txt reports the clinical similarity coefficients (Sim_BP, Sim_CC, Sim_MF) for each pair of PS (PSi_PSj), based on the shared biological features (i.e., the annotations of the disease gene products in each of the three sub-ontologies of Gene Ontology, namely Biological Process, BP; Cellular Component, CC; Molecular Function, MF). The PS are the nodes of the Biological Similarity Networks (BSN) reported in Figures S2 (BSN-BP), S3 (BSN-CC) and S4 (BSN-MF) of the manuscript, while the similarity coefficients represent the weight w of the edge connecting the two PS in the relevant BSN. In the general BSN (reported in Figure 3), an edge (linking a given pair of PS) is the edge with the highest w among the three edges that link the same PS pair in the three sub-ontology BSN (Sim_GO_max).

'Build_mm_and_PS.py': this is the Python script for the building the correct structure of 'morbid map' and 'PS'.

'Similarity_go_hpo.py': this is the Python script for the calculation of similarity. It requires complete ontologies, as also annotations of genes and diseases, from both GO and HPO.

