from collections import Counter
from math import log
from itertools import combinations
from time import clock

##########################################
# It is important to remember the [Typedef] at the end of the file "go-basic.obo"


# import database with .obo format
def Import_obo(infile): # import one db as dictionary
	D_name = dict() # create a dictionary for the names of terms
	D = dict() # create a dictionary for all terms
	f = open(infile,'r') # open the file ".obo"
	inp = f.read() # read it
	f.close() # and close it
	inp = inp.split('[Term]\n') # split each term
	for term in inp[1:]: # skip all before the first "[Term]" (file header)
		term = term.split('\n')
		ID = term[0][4:] # the term ID, as "HP:0000001"
		parent = [] # list to fill with parents of the term
		for line in term:
			if line.startswith("name:"): # name of the term
				D_name[ID] = line[6:] # add name to the dictionary
			elif line.startswith("is_a:"): # "is_a" indicates one parent
				parent.append(line[6:16])
		if parent != []: # if the list of parents is no empty
			D[ID] = parent # the dictionary for each term has a list of its parents
	return (D, D_name)


# invert key-value of dictionary (child-parent) -> (parent-child)
def Invert(D1):
	D2 = dict()
	for key in D1:
		D2[key]=[]
	for key in D1:
		for val in D1[key]:
			D2[val].append(key)
	return D2


# calculates descendants of a given term
def Descendant(parent, D): # D is the disctionary to use
	desc = set(parent) # set with descendants
	while parent != set([]): # when a parent has no childs (end reached)
		child = [] # all children of the parent are put in this list
		for term in parent:
			child.extend(D[term])
		parent = set(child) # children become the new parents
		desc = desc | parent
	return desc


# return a subset of original ontology
def Sub_ontology(D1, root): # root is the origin of the new ontology
	D2 = Invert(D1) # 1st inversion
	desc = Descendant([root], D2)
	D3 = dict()
	for term in desc:
		D3[term] = D2[term]
	D4 = Invert(D3) # 2nd inversion
	return D4


# calculate ancestors of a single term
def Ancestors(child, D): # D is the disctionary to use
	anc = set(child) # set with descendants
	while child != set([]): # when a parent has no childs (end reached)
		parent = [] # all children of the parent are put in this list
		for term in child:
			parent.extend(D[term])
		child = set(parent) # children become the new parents
		anc = anc | child
	return anc


# calculates ancestors for all terms
def Ancestors_of_all(D1):
	D2 = dict() # dictionary of all ancestors for each terms
	for term in D1:
		anc = Ancestors([term], D1) # calculate the ancestors
		D2[term] = anc # the set of ancestors
	return D2


# open, read and close the given file
def Open_r(infile):
	f = open(infile,'r')
	inp = f.readlines()
	f.close()
	return inp


# import database (tab separated)
def Import_db(infile, p1, p2): # import one db as dictionary
	D = dict()
	inp = Open_r(infile)
	for line in inp[1:]: # skip the first line (file header)
		el = line[:-1].split('\t')
		try:
			D[el[p1]].add(el[p2])
		except KeyError:
			D[el[p1]] = set([el[p2]]) # one Key = one set of terms
	return D


# remove prefix from the keys
def Remove_prefix(D1, x): # the prefix in D keys is from x to end [x:] 
	D2 = dict()
	for key in D1:
		D2[key[x:]] = D1[key] # remove the prefix "OMIM:"
	return D2


# reduce a dictionary values using anly valid values
def Reduce_D_values(D1, valid):
	D2 = dict()
	for key, val in D1.items():
		newval = val & valid # new values
		if newval != set([]):
			D2[key] = newval
	return D2


# reduce a dictionary keys using only valid keys
def Reduce_D_keys(D1, valid):
	D2 = dict()
	for key in D1:
		if key in valid:
			D2[key] = D1[key]
	return D2
	
	
# return a set of all values of dictionary
def Return_D_values(D):
	values = set()
	for val in D.values():
		for v in val:
			values.add(v)
	return values


# count the frequency of each terms and calculate the IC
def Calcul_Resnik(D_ann, D_anc):# annotations, ancestors
	all_terms = [] # all terms of annotation
	for val in D_ann.values():
		ancestors = [] # all ancestors
		for term in val:
			ancestors.extend(D_anc[term])
		ancestors = list(set(ancestors))
		all_terms.extend(ancestors)
	D_count = Counter(all_terms) # count the presence of each term
	D_ic = dict()
	N = len(D_ann) * 1.0 # number of diseases
	#M = -log(1/N) # maximal possible value
	for term, p in D_count.items():
		D_ic[term] = -log(p/N)
	return D_ic


# combine two dictionaries of Ancestors and IC
def Combine(D_ic, D_anc):
	D = dict()
	for key in D_ic:
		ic_anc = []
		for term in D_anc[key]:
			ic_anc.append((D_ic[term], term))
		D[key] = set(ic_anc)
	return D


# 
def Rebuild_annotation1(D_ann, D_hpo):
	D = dict()
	for el in D_ann: # element in annotation as disease or gene
		terms = D_ann[el]
		parent = list(terms)
		for term in terms:
			parent.extend(list(D_hpo[term]))
		D[el] = set(parent)
	return D


# build two dictionaries used for calculate similarity
def Rebuild_annotation2(D_ann, D_anc):
	D1 = dict()
	D2 = dict()
	for el in D_ann: # element in annotation as disease or gene
		anc1 = [] # first type of ancestors
		anc2 = [] # second type of ancestors
		for term in D_ann[el]:
			anc1.append(D_anc[term])
			anc2.extend(list(D_anc[term]))
		D1[el] = anc1
		D2[el] = set(anc2)
	return (D1, D2)


# calculate similarity between two entities	
def Calcul_sim(e1, e2, D1, D2): # two element: e1, e2. two dictionary: D1, D2.
	mica = []
	for term in D1[e1]:
		M = max(term & D2[e2])
		mica.append(M[0])
	for term in D1[e2]:
		M = max(term & D2[e1])
		mica.append(M[0])
	sim = sum(mica)/len(mica)
	return sim

