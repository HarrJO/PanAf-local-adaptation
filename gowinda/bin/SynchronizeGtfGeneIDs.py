#!/usr/bin/env python

import sys
import collections
from optparse import OptionParser, OptionGroup
from optparse import OptionParser,OptionGroup
import re

    #GO:0030518_steroid_hormone_receptor_signaling_pathway	CG2995	6	4	5.027203	-2.432819	5	5.0	1.000000
    #GO:0030518_steroid_hormone_receptor_signaling_pathway	CG1765	6	4	5.027203	-2.432819	5	5.0	1.000000
    #GO:0030518_steroid_hormone_receptor_signaling_pathway	CG32346	6	4	5.027203	-2.432819	5	5.0	1.000000
    #GO:0030518_steroid_hormone_receptor_signaling_pathway	CG1716	6	4	5.027203	-2.432819	5	5.0	1.000000
def get_geneset(inputfile):
	geneset=set([])
	for line in open(inputfile):
		line=line.rstrip()
		if(line.startswith("#")):
			continue
		a=line.split('\t')
		genestring=a[2]
		geneids=[]
		if(genestring.find(" ")):
			geneids.extend(genestring.split(" "))
		else:
			geneids.append(genestring)
		for gene in geneids:
			geneset.add(gene.lower())
	return geneset

def get_aliashash(inputfile):
	aliashash={}
	for line in open(inputfile):
		line=line.rstrip()
		genes=[]
		if(line.find("\t")):
			genes.extend(line.split("\t"))
		else:
			genes.add(line)
		for gene in genes:
			aliashash[gene.lower()]=genes
	return aliashash

def get_sync_geneid(geneid,geneset,aliaslist):
	novelgenid=geneid
	if(geneid.lower() not in aliaslist):
		return novelgenid
	for gene in aliaslist[geneid.lower()]:
		if gene.lower() in geneset:
			novelgenid=gene
	return novelgenid


parser = OptionParser()
parser.add_option("--gene-set", dest="geneset", help="gene set file")
parser.add_option("--gene-aliases", dest="aliaslist",help="List of gene aliases")
parser.add_option("--gtf",dest="gtf",help="the gtf file which should be synchronized with the alised gene IDs")
(options, args) = parser.parse_args()


geneset=get_geneset(options.geneset)
aliaslist=get_aliashash(options.aliaslist)

for line in open(options.gtf):
	line=line.rstrip()
	a=line.split('\t')
	if(len(a)!=9):
		continue
	attr=a[8]
	geneid=re.search(r'gene_id "([^"]+)";',attr).group(1)
	novelgeneid=get_sync_geneid(geneid,geneset,aliaslist)
	test=re.sub(r'gene_id "([^"]+)";','gene_id "'+novelgeneid+'";',attr)
	a[8]=test
	print '\t'.join(a)