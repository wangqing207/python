#!
import re, glob

myfiles = glob.glob("TCGA-*.txt")
print(myfiles)
exit
for myfile in myfiles:
    outfile = open("test" + myfile, "w")
    with open(myfile) as infile:
        for line in infile:
            cols =  line.strip().split("\t")
            outfile.write("\t".join(cols[0:2]) + "\n")
