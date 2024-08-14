import sys
import collections
from optparse import OptionParser, OptionGroup
from optparse import OptionParser,OptionGroup
import re

# 0             1        2          3       4       5           6       7       8
# chr4	dm3_flyBaseGene	exon	251356	251521	0.000000	+	.	gene_id "CG1674-RB"; transcript_id "CG1674-RB"; 
# chr4	dm3_flyBaseGene	start_codon	252580	252582	0.000000	+	.	gene_id "CG1674-RB"; transcript_id "CG1674-RB"; 
# chr4	dm3_flyBaseGene	CDS	252580	252603	0.000000	+	0	gene_id "CG1674-RB"; transcript_id "CG1674-RB"; 



def getMrnaDecoder(mrnalist):
	mdec={}
	for mrna in mrnalist:
		annotation=mrna[8]
		id=re.search(r"ID=([^;]+)",annotation).group(1)
		parent=re.search(r"Parent=([^;]+)",annotation).group(1)
		mdec[id]=parent
	return mdec


#Author: Robert Kofler
parser = OptionParser()
parser.add_option("--input", dest="input", help="The input file; A FlyBase .gff")
(options, args) = parser.parse_args()
exons=[]
cds=[]
mrnas=[]
for line in open(options.input):
	line=line.rstrip()

	a=line.split('\t');
	if(len(a)!=9):
		continue
	if(a[2]=="exon"):
		exons.append(a)
	if(a[2]=="CDS"):
		cds.append(a)
	if(a[2]=="mRNA"):
		mrnas.append(a)

mdec=getMrnaDecoder(mrnas)

for an in exons+cds:
	annotation=an[8]
	parent=re.search(r"Parent=([^;]+)",annotation).group(1)
	if(parent.find(',')!=-1):
		parent =parent.split(',')[0]
	if(parent not in mdec):
		continue
	geneparent=mdec[parent]
	an[8] = "gene_id \""+geneparent+"\";"
	print("\t".join(an))
		
	
	

#
#2L	FlyBase	exon	7529	8116	.	+	.	ID=FBgn0031208:1;Name=CG11023:1;Parent=FBtr0300689,FBtr0300690;parent_type=mRNA
#2L	FlyBase	exon	8193	9484	.	+	.	ID=FBgn0031208:3;Name=CG11023:3;Parent=FBtr0300689;parent_type=mRNA
#2L	FlyBase	exon	8193	8589	.	+	.	ID=FBgn0031208:2;Name=CG11023:2;Parent=FBtr0300690;parent_type=mRNA

#2L	FlyBase	CDS	8193	8610	.	+	1	ID=CDS_FBgn0031208:3_854;Name=CG11023-cds;Parent=FBtr0300689;parent_type=mRNA
#2L	FlyBase	CDS	8193	8589	.	+	1	ID=CDS_FBgn0031208:2_854;Name=CG11023-cds;Parent=FBtr0300690;parent_type=mRNA
#2L	FlyBase	CDS	8668	9276	.	+	0	ID=CDS_FBgn0031208:4_854;Name=CG11023-cds;Parent=FBtr0300690;parent_type=mRNA

#2L	FlyBase	mRNA	7529	9484	.	+	.	ID=FBtr0300690;Name=CG11023-RC;Parent=FBgn0031208;Dbxref=REFSEQ:NM_175941,FlyBase_Annotation_IDs:CG11023-RC;score_text=Strongly Supported;score=15
#2L	FlyBase	mRNA	9839	21372	.	-	.	ID=FBtr0078166;Name=l(2)gl-RB;Parent=FBgn0002121;Alias=CG2671-RB;Dbxref=FlyBase_Annotation_IDs:CG2671-RB,REFSEQ:NM_164348;score_text=Strongly Supported;score=11
#2L	FlyBase	mRNA	9839	21372	.	-	.	ID=FBtr0078167;Name=l(2)gl-RD;Parent=FBgn0002121;Alias=CG2671-RD;Dbxref=FlyBase_Annotation_IDs:CG2671-RD,REFSEQ:NM_164350;score_text=Weakly Supported;score=3
#2L	FlyBase	mRNA	9839	21372	.	-	.	ID=FBtr0078168;Name=l(2)gl-RE;Parent=FBgn0002121;Alias=CG2671-RE;Dbxref=FlyBase_Annotation_IDs:CG2671-RE,REFSEQ:NM_164351;score_text=Weakly Supported;score=3
#2L	FlyBase	mRNA	9839	21372	.	-	.	ID=FBtr0078169;Name=l(2)gl-RF;Parent=FBgn0002121;Alias=CG2671-RF;Dbxref=FlyBase_Annotation_IDs:CG2671-RF,REFSEQ:NM_164352;score_text=Weakly Supported;score=3

#2L	FlyBase	gene	66721	71390	.	+	.	ID=FBgn0067779;Name=dbr;fullname=debra;Alias=FBgn0028472,FBgn0043386,CG11371,Dbr,EP456,BcDNA:LD26519,EP0456,Debra;Ontology_term=SO:0000010,SO:0000087,GO:0005771,GO:0005634,GO:0008270;Dbxref=FlyBase:FBan0011371,FlyBase_Annotation_IDs:CG11371,GB_protein:AAF51565,GB:AA942362,GB:AF184228,GB_protein:AAD55739,GB:AI109558,GB:AI402147,GB:AQ026007,GB:AQ073408,GB:AW942847,GB:AY047576,GB_protein:AAK77308,GB:CL527762,GB:CZ468654,GB:CZ470443,GB:CZ470758,GB:CZ471260,GB:CZ474676,GB:CZ474677,GB:CZ477845,GB:CZ482562,GB:CZ487746,UniProt/TrEMBL:Q961T5,UniProt/TrEMBL:Q9V474,INTERPRO:IPR012934,EntrezGene:33161,InterologFinder:33161,BIOGRID:59424,FlyAtlas:CG11371-RB,GenomeRNAi:33161;gbunit=AE014134;derived_computed_cyto=21B1-21B1
#2L	FlyBase	gene	71757	76211	.	+	.	ID=FBgn0031213;Name=galectin;fullname=galectin;Alias=CG11372,Dmgal;Ontology_term=SO:0000010,SO:0000087,GO:0016936,GO:0005829,GO:0015143;Dbxref=FlyBase:FBan0011372,FlyBase_Annotation_IDs:CG11372,GB_protein:AAF51564,GB_protein:ACZ94130,GB_protein:ACZ94129,GB:AF338142,GB_protein:AAL87743,GB:AY071294,GB_protein:AAL48916,GB:BI236456,GB:CZ467438,GB:CZ467439,GB:CZ468006,GB:CZ481733,GB:CZ481734,GB:CZ485874,GB:CZ485875,UniProt/TrEMBL:Q9VPI6,INTERPRO:IPR001079,INTERPRO:IPR008985,INTERPRO:IPR013320,EntrezGene:33162,UniProt/TrEMBL:E1JHP9,UniProt/TrEMBL:E1JHQ0,InterologFinder:33162,BIOGRID:59425,FlyAtlas:CG11372-RA,GenomeRNAi:33162;gbunit=AE014134;derived_computed_cyto=21B1-21B1
#2L	FlyBase	gene	76348	77783	.	+	.	ID=FBgn0031214;Name=CG11374;Ontology_term=SO:0000010,SO:0000087,GO:0016936,GO:0015143;Dbxref=FlyBase:FBan0011374,FlyBase_Annotation_IDs:CG11374,GB_protein:AAF51563,UniProt/TrEMBL:Q9VPI7,INTERPRO:IPR001079,INTERPRO:IPR008985,INTERPRO:IPR013320,EntrezGene:33163,InterologFinder:33163,FlyAtlas:CG11374-RA,GenomeRNAi:33163;gbunit=AE014134;derived_computed_cyto=21B1-21B1
#2L	FlyBase	gene	82456	87382	.	-	.	ID=FBgn0002931;Name=net;fullname=net;Alias=FBgn0031215,CG11450,Shout,Group IId,shout,bHLHc42;Ontology_term=SO:0000010,SO:0000087,GO:0005634,GO:0003700,GO:0008586,GO:0006355,GO:0007474,GO:0000122;Dbxref=FlyBase:FBan0011450,FlyBase_Annotation_IDs:CG11450,GB_protein:AAF51562,GB:AF303350,GB_protein:AAK14073,GB:CZ478919,GB:CZ478940,UniProt/TrEMBL:Q9VPI8,INTERPRO:IPR011598,UniProt/TrEMBL:Q7KH24,EntrezGene:45339,InterologFinder:45339,BIOGRID:69604,FlyAtlas:CG11450-RA,GenomeRNAi:45339,INTERACTIVEFLY:/gene/net1.htm;gbunit=AE014134;derived_computed_cyto=21B1-21B1
