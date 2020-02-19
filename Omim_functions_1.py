
# 1. open, read and close one given file
def Read_file(infile):
	f = open(infile,'r')
	inp = f.readlines()
	f.close()
	return inp


# 2. retrieve the file "morbidmap.txt"
def Read_morbidmap1(f): # read the file "f"
	inp = Read_file(f)
	OMIM = []
	for line in inp:
		if line.startswith('#'):
			continue
		line = line[:-1].split('\t')
		name = line[0] # disease name
		mimG = line[2] # MIM of Gene
		name = name.rstrip(', ')
		num = name[-3:]
		name = name[:-3] # remove (x) at the end
		name = name.rstrip(', ') # remove ',' or ' ' at the end
		mimD = name[-6:] # disease-MIM
		if mimD.isdigit():
			name = name[:-6] # remove the number
			name = name.rstrip(', ') # remove ',' or ' ' at the end
		else:
			mimD = mimG # if mim1 is not present use mim2
		OMIM.append([name, mimD, num, mimG])
	return OMIM


# 3. extract information from the file "Homo_sapiens.gene_info"
def Read_gene_info(f):
	D = dict()
	inp = Read_file(f)
	for line in inp:
		line = line[:-1].split('\t')
		if line[0] != "9606": # different from human
			continue
		dbXrefs = line[5].split('|')
		for el in dbXrefs:
			if el[:4] == 'MIM:':
				entrez = line[1] # the "GeneID" (entrez) is in column[1]
				symbol = line[2] # the "Symbol" is in column[2]
				prot = line[9] # protein-coding or not
				D[el[4:]] = [entrez, symbol, prot]
	return D


# 4. combine together the two files
def New_morbidmap(morbidmap, D, f):
	out = open(f, 'w')
	out.write('#Disease_name\tDisease_mim\tkey\tGene_mim\tentrez\tsymbol\ttype\n')
	for L1 in morbidmap:
		try:
			L2 = D[L1[3]]
		except KeyError:
			L2 = ['-','-','-']
		s = '\t'.join(L1 + L2)
		out.write(s + '\n')
	out.close()
	

# 5. retrieve the file "morbidmap.txt"
def Read_morbidmap2(f): # read the file "f"
	inp = Read_file(f)
	D = dict()
	for line in inp:
		line = line[:-1].split('\t')
		if line[0].startswith(('{', '?', '[')):
			continue
		if line[2] != '(3)':
			continue
		if line[6] != 'protein-coding':
			continue
		mim = line[1]
		entrez = line[4]
		try:
			D[mim].add(entrez)
		except KeyError:
			D[mim] = set([entrez])
	return D

		
# 6. open and extract information from the file "phenotypicSeries.txt"
def Read_phen_ser(f, valid): # read the file "f"
	D = dict()
	inp = Read_file(f)
	for line in inp:
		if line.startswith('#'):
			continue # ignore lines with '#'.
		line = line.split('\t')
		ps = line[0]
		disease = line[1]
		if disease not in valid:
			continue
		try:
			D[ps].add(disease)
		except KeyError:
			D[ps] = set([disease])
	return D


# 7. return a set of all values of a dictionary
def Return_values(D):
	val = []
	for el in D.values():
		val.extend(el)
	return set(val)


# 8. find cases of LH in morbid map dictionary 
def Find_LH(D):
	LH = [el for el in D if len(D[el]) > 1]
	return set(LH)


# 9. add LH to PS
def Add_LH_to_PS(D1, D2):
	mim_LH = Find_LH(D1)
	mim_PS = Return_values(D2)
	S = mim_LH - mim_PS
	for mim in S:
		D2['LH'+mim] = [mim]
	return D2


# 10. print out the new PS file 
def New_PS(D1, D2, f):
	out = open(f, 'w')
	out.write("#PS_ID\tMIM_ID\tentrez\n")
	for ps in D1:
		for mim in D1[ps]:
			try:
				gene = D2[mim]
			except KeyError:
				continue
			for g in gene:
				out.write("%s\t%s\t%s\n" % (ps, mim, g))
	out.close()
	

# 11. import database (tab separated)
def Import_db(f, p1, p2, n): # import one db as dictionary
	D = dict()
	inp = Read_file(f)
	for line in inp[n:]: # skip n lines (file header)
		el = line[:-1].split('\t')
		try:
			D[el[p1]].add(el[p2])
		except KeyError:
			D[el[p1]] = set([el[p2]]) # one Key = one set of terms
	return D


# 12. import Omim annotation (tab separated)
def Import_mim2hp(f): # import one db as dictionary
	D = dict()
	inp = Read_file(f)
	for line in inp[1:]: # skip 
		line = line[:-1].split('\t')
		if line[3] == 'HP:0000118': # 
			continue
		mim = line[0][5:]
		try:
			D[mim].add(line[3])
		except KeyError:
			D[mim] = set([line[3]]) # one Key = one set of terms
	return D


def Import_gene2go(f):
	D = dict()
	inp = Read_file(f)
	for line in inp:
		line = line[:-1].split('\t')
		if line[0] != "9606": # different from human
			continue
		if line[2] == 'GO:0008150': # "biological_process"
			continue
		if line[2] == 'GO:0005575': # "cellular_component"
			continue
		if line[2] == 'GO:0003674': # "molecular_function"
			continue
		try:
			D[line[1]].add(line[2])
		except KeyError:
			D[line[1]] = set([line[2]]) # one Key = one set of terms
	return D


def Print_dict(D, outfile):
	out = open(outfile, 'w')
	for key in D:
		for val in D[key]:
			out.write('%s\t%s\n' % (key,val))
	out.close()
	

# 13. reduce a dictionary values using only valid values
def Reduce_D_values(D1, valid):
	D2 = dict()
	for key, val in D1.items():
		newval = val & valid # new values
		if newval != set([]):
			D2[key] = newval
	return D2
