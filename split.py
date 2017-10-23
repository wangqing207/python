#!usr/bin/python
import re
import sys

a = open(sys.argv[1],'r')
gindex=[4,5,6]
for line in a:
    cols=line.strip().split("\t")
    if line.startswith("\"Target_ID\""):
        #gindex=[i for i in range(len(cols)) if re.match("UCSC_REFGENE",cols[i])]
        print "\t".join(cols[:gindex[2]+1])+"\t"+"gene_names"+"\t"+"transcript_id"+"\t"+"body_type"+"\t"+"\t".join(cols[gindex[2]+1:])
    else:
#	print str(len(cols))+"\t"+str(gindex[1])
        if len(cols)<gindex[0]+1:
            print "\t".join(cols)
        else:
            if cols[gindex[0]]!="":
                genelist=cols[gindex[0]].strip().split(";")
                tranid=cols[gindex[1]].strip().split(";")
                body=cols[gindex[2]].strip().split(";")
                for i in range(len(genelist)):
                    print "\t".join(cols[:gindex[2]+1])+"\t"+genelist[i]+"\t"+tranid[i]+"\t"+body[i]+"\t"+"\t".join(cols[gindex[2]+1:])
            else:
                print "\t".join(cols[:gindex[2]+1])+"\t"+"\t".join(cols[gindex[2]+1:]) 
        




