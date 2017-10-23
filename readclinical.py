import json
import sys

data = json.loads(open('clinical.project-TCGA-KIRC.2017-08-23T05-38-34.681887.json').read())

for n, d in enumerate(data):
	diagnoses = d.get('diagnoses', [{}])[0]
	case_id = d.get('case_id', '')
	demographic = d.get('demographic', {})
	exposures = d.get('exposures', [{}])[0]
	
	#write title
	if n == 0:
		sys.stdout.write("\t".join(str(i) for i in diagnoses.keys()) + "\t")
		sys.stdout.write("case_id\t")
		sys.stdout.write("\t".join(str(i) for i in demographic.keys()) + "\t")
		sys.stdout.write("\t".join(str(i) for i in exposures.keys()) + "\n")
	
	if diagnoses != {}:
		sys.stdout.write("\t".join(str(i) for i in diagnoses.values()) + "\t")
	else:
		sys.stdout.write("\t".join([""]*22) + "\t")
	sys.stdout.write(case_id + "\t")
	if demographic != {}:
		sys.stdout.write("\t".join(str(i) for i in demographic.values()) + "\t")
	else:
		sys.stdout.write("\t".join([""]*10) + "\t")
	if exposures!= {}:
		sys.stdout.write("\t".join(str(i) for i in exposures.values()) + "\n")
	else:
		sys.stdout.write("\t".join([""]*12) + "\n")