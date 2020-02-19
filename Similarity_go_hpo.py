from Omim_tools_2 import *
	
t0 = clock()

# It is important to remember the [Typedef] at the end of the file "go-basic.obo"

inp = (
('in_data/hp.obo',       'HP:0000118', 'out_data/new_mim2hp.txt',  1, 'out_data/PS_Sim_HP.txt'),
('in_data/go-basic.obo', 'GO:0008150', 'out_data/new_gene2go.txt', 2, 'out_data/PS_Sim_BP.txt'),
('in_data/go-basic.obo', 'GO:0005575', 'out_data/new_gene2go.txt', 2, 'out_data/PS_Sim_CC.txt'),
('in_data/go-basic.obo', 'GO:0003674', 'out_data/new_gene2go.txt', 2, 'out_data/PS_Sim_MF.txt'))

inp = inp[3] # this is to change, depending by the input you want (0, 1, 2 or 3)

obo = inp[0]
root = inp[1]
annotation = inp[2]
p = inp[3]
out_file = inp[4]

ontology = Import_obo(obo) # file with [Typedef] removed

D_term = ontology[0] # this is the entire ontology
D_term_name = ontology[1] # names of the terms

D_term['HP:0000001'] = [] # add the root of hp that is not present
D_term['GO:0008150'] = [] # "biological_process"
D_term['GO:0005575'] = [] # "cellular_component"
D_term['GO:0003674'] = [] # "molecular_function"

print 'original ontology:', len(D_term)

D_term = Sub_ontology(D_term, root) # new ontology with only son of the root
print 'new ontology:', len(D_term)
valid_terms = set(D_term.keys()) # valid terms

############################################

#D_annotation = Import_db('OMIM_ALL_FREQUENCIES_diseases_to_genes_to_phenotypes.txt', 0, 3)
#D_annotation = Remove_prefix(D_annotation, 5) # remove the "OMIM:" prefix before the disease ID

D_annotation = Import_db(annotation, 0, 1)

print 'first len D_annotation', len(D_annotation)

D_annotation = Reduce_D_values(D_annotation, valid_terms) # keep only the valid hp terms
print 'second len D_annotation', len(D_annotation)

D_anc = Ancestors_of_all(D_term)
print 'len(D_anc)', len(D_anc)

D_ic = Calcul_Resnik(D_annotation, D_anc)
print 'len(D_ic)', len(D_ic)

D_ic_anc = Combine(D_ic, D_anc)
print 'len(D_ic_anc)', len(D_ic_anc)

#D_ic.clear()
#D_anc.clear()

D_ps = Import_db("out_data/new_PS.txt", 0, p) # PS -> gene
ps_values = Return_D_values(D_ps)
print 'len(ps_values)', len(ps_values)

D_annotation = Reduce_D_keys(D_annotation, ps_values)
print 'third len D_annotation', len(D_annotation)

Comb = list(combinations(D_annotation.keys(), 2))

D_ann1 = Rebuild_annotation1(D_annotation, D_term)

new_ann = Rebuild_annotation2(D_annotation, D_ic_anc)
D_ann2 = new_ann[0]
D_ann3 = new_ann[1]

all_sim = []
for e1, e2 in Comb:
	sim = Calcul_sim(e1, e2, D_ann2, D_ann3)
	all_sim.append(sim)

M = max(all_sim)

#### print out normalized similarities ####

out = open(out_file, 'w')
for e, s in zip(Comb, all_sim):
	sim = s / M
	out.write('%s\t%s\t%.4f\n' % (e[0], e[1], sim))
out.close()


#################
t1 = clock()
print "Elapsed time:", t1 - t0
#################

