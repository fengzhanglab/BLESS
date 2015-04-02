import sys
from Bio.Seq import Seq
from Bio import SeqIO

IONA = 'GAGGTAGTA' # BLESS handle
TRP = 'GAGGTAGTA'
LENOUT = 30

def truncate_seqs_fwd(filename):
	handle = open(filename,"rU")
	print filename + ' opened.'
	filename2 = filename.replace('.fastq','')
	handleout = open(filename2+'_genomic.fastq','w')

	rec_count = 0
	for record in SeqIO.parse(handle, "fastq"): #seqio is a package in biopython
		id = record.id #biopython object holds fastq
		seq = str(record.seq)
		#print seq
		qual = str(record.letter_annotations)
		#print qual
		iona_pos = str.find(seq,IONA)#is string present, rfind(start from back) changed to find(start from front) 

		record2 = []
		if iona_pos > -1:
			record2 = record[(iona_pos+len(IONA)):len(seq)]

		if len(record2) >= LENOUT:
			SeqIO.write(record2[0:LENOUT],handleout,"fastq")

	return

def truncate_seqs_rev(filename):
	handle = open(filename,"rU")
	print filename + ' opened.'
	filename2 = filename.replace('.fastq','')
	handleout = open(filename2+'_genomic.fastq','w')

	rec_count = 0
	for record in SeqIO.parse(handle, "fastq"):
		id = record.id
		seq = str(record.seq)
		#print seq
		qual = str(record.letter_annotations)
		#print qual
		trp_pos = str.find(seq,TRP)

		record2 = []
		if trp_pos > -1:
			record2 = record[(trp_pos+len(TRP)):len(seq)]

		if len(record2) >= LENOUT:
			SeqIO.write(record2[0:LENOUT],handleout,"fastq")

	return	

filename = sys.argv[1]
if str.find(sys.argv[2],'R1'):
	print filename
	truncate_seqs_fwd(filename)
elif str.find(sys.argv[2],'R2'):
	print filename
	truncate_seqs_rev(filename)	
