import sys
import time
import random
from Bio import Entrez
ids=[]
infile=sys.argv[1]
for line in open(infile,'r'):
	line=line.strip()
	ids.append(line)
for i in range(1,len(ids)):
#  t = random.randrange(0,5)
	handle = ''
	Entrez.efetch(db="nucleotide", id=ids[i],rettype="fasta",email="jmzeng1314@163.com")
#  time.sleep(t)
	print(handle.read())