#!/usr/bin/env python
import sys
import collections
from optparse import OptionParser, OptionGroup
from optparse import OptionParser,OptionGroup
import re

# 0             1        2          3       4       5           6       7       8
# chr4	dm3_flyBaseGene	exon	251356	251521	0.000000	+	.	gene_id "CG1674-RB"; transcript_id "CG1674-RB"; 
# chr4	dm3_flyBaseGene	start_codon	252580	252582	0.000000	+	.	gene_id "CG1674-RB"; transcript_id "CG1674-RB"; 
# chr4	dm3_flyBaseGene	CDS	252580	252603	0.000000	+	0	gene_id "CG1674-RB"; transcript_id "CG1674-RB"; 
#ID=FBgn0051973;Name=Cda5;fullname=Chitin deacetylase-like 5;Alias=FBgn0031210,FBgn0031211,FBgn0063656,CG2761,CG2776,BcDNA:RH43162,DmCDA5,CG31973;Ontology_term=SO:0000010,SO:0000087,GO:0006030,GO:0016810,GO:0005576,GO:0008061;Dbxref=FlyBase:FBan0031973,FlyBase_Annotation_IDs:CG31973,GB_protein:AAF51567,GB_protein:ABI31281,GB_protein:AAF51568,GB_protein:ABV53594,GB:AY129461,GB_protein:AAM76203,GB:BI613902,GB:CZ468962,GB:CZ486696,GB:CZ486697,UniProt/TrEMBL:Q8MQI4,UniProt/TrEMBL:Q9VPI3,UniProt/TrEMBL:Q9VPI4,INTERPRO:IPR002509,INTERPRO:IPR002557,INTERPRO:IPR011330,EntrezGene:33158,UniProt/TrEMBL:Q0E8V4,UniProt/TrEMBL:A8DYS5,InterologFinder:33158,BIOGRID:59423,FlyAtlas:CG31973-RA,GenomeRNAi:33158;gbunit=AE014134;derived_computed_cyto=21B1-21B1
#2L	FlyBase	gene	66721	71390	.	+	.	ID=FBgn0067779;Name=dbr;fullname=debra;Alias=FBgn0028472,FBgn0043386,CG11371,Dbr,EP456,BcDNA:LD26519,EP0456,Debra;Ontology_term=SO:0000010,SO:0000087,GO:0005771,GO:0005634,GO:0008270;Dbxref=FlyBase:FBan0011371,FlyBase_Annotation_IDs:CG11371,GB_protein:AAF51565,GB:AA942362,GB:AF184228,GB_protein:AAD55739,GB:AI109558,GB:AI402147,GB:AQ026007,GB:AQ073408,GB:AW942847,GB:AY047576,GB_protein:AAK77308,GB:CL527762,GB:CZ468654,GB:CZ470443,GB:CZ470758,GB:CZ471260,GB:CZ474676,GB:CZ474677,GB:CZ477845,GB:CZ482562,GB:CZ487746,UniProt/TrEMBL:Q961T5,UniProt/TrEMBL:Q9V474,INTERPRO:IPR012934,EntrezGene:33161,InterologFinder:33161,BIOGRID:59424,FlyAtlas:CG11371-RB,GenomeRNAi:33161;gbunit=AE014134;derived_computed_cyto=21B1-21B1
#2L	FlyBase	gene	71757	76211	.	+	.	ID=FBgn0031213;Name=galectin;fullname=galectin;Alias=CG11372,Dmgal;Ontology_term=SO:0000010,SO:0000087,GO:0016936,GO:0005829,GO:0015143;Dbxref=FlyBase:FBan0011372,FlyBase_Annotation_IDs:CG11372,GB_protein:AAF51564,GB_protein:ACZ94130,GB_protein:ACZ94129,GB:AF338142,GB_protein:AAL87743,GB:AY071294,GB_protein:AAL48916,GB:BI236456,GB:CZ467438,GB:CZ467439,GB:CZ468006,GB:CZ481733,GB:CZ481734,GB:CZ485874,GB:CZ485875,UniProt/TrEMBL:Q9VPI6,INTERPRO:IPR001079,INTERPRO:IPR008985,INTERPRO:IPR013320,EntrezGene:33162,UniProt/TrEMBL:E1JHP9,UniProt/TrEMBL:E1JHQ0,InterologFinder:33162,BIOGRID:59425,FlyAtlas:CG11372-RA,GenomeRNAi:33162;gbunit=AE014134;derived_computed_cyto=21B1-21B1
#2L	FlyBase	gene	76348	77783	.	+	.	ID=FBgn0031214;Name=CG11374;Ontology_term=SO:0000010,SO:0000087,GO:0016936,GO:0015143;Dbxref=FlyBase:FBan0011374,FlyBase_Annotation_IDs:CG11374,GB_protein:AAF51563,UniProt/TrEMBL:Q9VPI7,INTERPRO:IPR001079,INTERPRO:IPR008985,INTERPRO:IPR013320,EntrezGene:33163,InterologFinder:33163,FlyAtlas:CG11374-RA,GenomeRNAi:33163;gbunit=AE014134;derived_computed_cyto=21B1-21B1
#2L	FlyBase	gene	82456	87382	.	-	.	ID=FBgn0002931;Name=net;fullname=net;Alias=FBgn0031215,CG11450,Shout,Group IId,shout,bHLHc42;Ontology_term=SO:0000010,SO:0000087,GO:0005634,GO:0003700,GO:0008586,GO:0006355,GO:0007474,GO:0000122;Dbxref=FlyBase:FBan0011450,FlyBase_Annotation_IDs:CG11450,GB_protein:AAF51562,GB:AF303350,GB_protein:AAK14073,GB:CZ478919,GB:CZ478940,UniProt/TrEMBL:Q9VPI8,INTERPRO:IPR011598,UniProt/TrEMBL:Q7KH24,EntrezGene:45339,InterologFinder:45339,BIOGRID:69604,FlyAtlas:CG11450-RA,GenomeRNAi:45339,INTERACTIVEFLY:/gene/net1.htm;gbunit=AE014134;derived_computed_cyto=21B1-21B1
#2L	FlyBase	gene	94752	102086	.	+	.	ID=FBgn0031216;Name=CG11376;Alias=FBgn0025528,anon-EST:Posey54,Dm zir,FBgn 31216;Ontology_term=SO:0000010,SO:0000087,GO:0005525,GO:0005085,GO:0051020;Dbxref=FlyBase:FBan0011376,FlyBase_Annotation_IDs:CG11376,GB_protein:AAF51561,GB:AA391952,GB:AF083301,GB:AW942162,GB:AZ877629,GB:BT015273,GB_protein:AAT94502,GB:CZ468078,GB:CZ468665,GB:CZ469343,GB:CZ470525,GB:CZ470937,GB:CZ480810,GB:CZ480811,GB:CZ483732,UniProt/TrEMBL:Q9VPI9,EntrezGene:33165,INTERPRO:IPR021816,INTERPRO:IPR010703,InterologFinder:33165,BIOGRID:59427,FlyAtlas:CG11376-RA,GenomeRNAi:33165;gbunit=AE014134;derived_computed_cyto=21B1-21B2

def splitgenestring(genestring):
	genelist=[]
	if(genestring.find(',')==-1):
		genelist.append(genestring)
	else:
		genelist.extend(genestring.split(','))
	return genelist
	

#Author: Robert Kofler
parser = OptionParser()
parser.add_option("--input", dest="input", help="The input file; A .gff")
(options, args) = parser.parse_args()


for line in open(options.input):
	line=line.rstrip()

	a=line.split('\t');
	if(len(a)!=9):
		continue
	if(a[2]!="gene"):
		continue
	attr=a[8]
	
	genelist=[]
	idsearch=re.search(r"ID=([^;]+)",attr)
	namesearch=re.search(r"Name=([^;]+)",attr)
	aliassearch=re.search(r"Alias=([^;]+)",attr)
	fullnamesearch=re.search(r"fullname=([^;]+)",attr)
	
	if(idsearch):
		genelist.extend(splitgenestring(idsearch.group(1)))
	if(namesearch):
		genelist.extend(splitgenestring(namesearch.group(1)))
	if(aliassearch):
		genelist.extend(splitgenestring(aliassearch.group(1)))
	if(fullnamesearch):
		genelist.extend(splitgenestring(fullnamesearch.group(1)))
	
	print "\t".join(genelist)

    