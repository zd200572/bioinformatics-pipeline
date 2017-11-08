fin = 
# Import the modules for interfacing with BLAST and parsing the output
from Bio.Blast import NCBIWWW, NCBIXML

# Blast the sequence of interest (in this case using the accession number
result_handle = NCBIWWW.qblast("blastn", "nr", filename)

# Parse the resulting output
blast_record = NCBIXML.read(result_handle)

# Loop over the alignments printing some output of interest
E_VALUE_THRESH = 0.004
for alignment in blast_record.alignments:
    result = alignment.title
    print('gi no.: '+str(result.split()))
    print('gi-desc: '+' '.join(str(result.split())))
    #print
##    for hsp in alignment.hsps:
##        if hsp.expect < E_VALUE_THRESH:
##            print
##            print '****Alignment****'
##            print 'sequence:', alignment.title
##            print 'length:', alignment.length
##            print 'e value:', hsp.expect
##            print hsp.query + '...'
##            print hsp.match + '...'
##            print hsp.sbjct + '...' 