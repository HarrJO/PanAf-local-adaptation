import sys
import collections
from optparse import OptionParser, OptionGroup
from optparse import OptionParser,OptionGroup
import re

    #GO:0030518_steroid_hormone_receptor_signaling_pathway	CG2995	6	4	5.027203	-2.432819	5	5.0	1.000000
    #GO:0030518_steroid_hormone_receptor_signaling_pathway	CG1765	6	4	5.027203	-2.432819	5	5.0	1.000000
    #GO:0030518_steroid_hormone_receptor_signaling_pathway	CG32346	6	4	5.027203	-2.432819	5	5.0	1.000000
    #GO:0030518_steroid_hormone_receptor_signaling_pathway	CG1716	6	4	5.027203	-2.432819	5	5.0	1.000000

parser = OptionParser()
parser.add_option("--input", dest="input", help="input file")
(options, args) = parser.parse_args()

trans=collections.defaultdict(lambda:set([]))

for line in open(options.input):
	line=line.rstrip()
	a=line.split('\t')
	gocat=a[0]
	geneid=a[1]
	if('|' not in geneid):
		trans[gocat].add(geneid)


for gocat, geneids in trans.items():
	# GO:0030518_steroid_hormone_receptor_signaling_pathway
	m=re.search(r"([^_]+)_(.*)",gocat)
	goid=m.group(1)
	desc=m.group(2)
	genestring=" ".join(geneids)
	toprint=[]
	toprint.append(goid)
	toprint.append(desc)
	toprint.append(genestring)
	print("\t".join(toprint))
	
	
	#geneid=re.search(r'gene_id "([^"]+)";',a[8]).group(1)