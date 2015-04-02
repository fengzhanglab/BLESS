import sys
import difflib

def find_seq_terminus(seq,seq0,idthresh,idsumthresh,mmthresh,mmsumthresh):
##################################################
#try using difflib for truncation mapping, or alternatively, try ssaha
#allow deletion or insertion lengths of no more than 1nt
#(normally illumina sequencing mismatches or indels will not exceed 1nt)
##################################################
	#idthresh = 1 
	#idsumthresh = 1 #indels are rare in illumina sequencing
	#mmthresh = 1
	#mmsumthresh = 2 #mismatches are more common in illumina sequencing

	s = difflib.SequenceMatcher(None, seq, seq0, autojunk=False)
	a = s.get_opcodes()
	
	print seq
	print seq0
	print a

	subterm = []
	idsum = 0
	mmsum = 0
	flag = 0
	ct = 0
	for sub in a:
		##print sub
		lensubseq = sub[2]-sub[1]
		lensubseq0 = sub[4]-sub[3]

		#move to sequence match start, consider bases removed from front of seq0 as insertion/deletion
		if (flag==0):
			#treat matches of len <= idthesh at start of seqeunce as deletions
			if (hash(sub[0])==hash('equal')) and (lensubseq0 > idthresh):
				##print 'equal and len greater than idthresh'				
				flag = 1
				idsum = idsum + sub[3]
			else:
				##print 'NOT equal and len greater than idthresh'
				continue

		if hash(sub[0])==hash('equal'):
			##print 'equal'
			subterm.append(sub[2])
		elif hash(sub[0])==hash('replace'):
			if max([lensubseq, lensubseq0])<=mmthresh:
				##print 'replace of len less than mmthresh'
				subterm.append(sub[2])
				mmsum = mmsum + max([lensubseq, lensubseq0])
			else:
				##print 'NOT replace of len less than mmthresh'
				break	
		elif hash(sub[0])==hash('insert'):
			if lensubseq0<=idthresh:
				##print 'insert of len less than idthresh'
				subterm.append(sub[2])
				idsum = idsum + lensubseq0
			else:
				##print 'NOT insert of len less than idthresh'
				break
		elif hash(sub[0])==hash('delete'):
			if lensubseq<=idthresh:
				##print 'delete of len less than idthresh'
				subterm.append(sub[2])
				idsum = idsum - lensubseq
			else:
				##print 'NOT delete of len less than idthresh'
				break

		if (abs(idsum) > idsumthresh) or (abs(mmsum) > mmsumthresh):
			##print 'idmmsumthresh exceeded'
			break
		
		ct = ct+1
		##print 'COUNT IS ' + str(ct)

	#return last subterminus for a series of subtermini containing:
	#1. individual insertions and deletions of length <= idthresh
	#2. total insertions and deletions compared to template string are <= idsumthresh 
	return subterm[ct-1]

def diff_find(seq,seq0,idthresh,idsumthresh,mmthresh,mmsumthresh):
##################################################
#NOTE:
#fundemental assumption of difflib violated. does not find longest contiguous matching sequence
#Therefore, need to write specific hashing algorithm for this purpose
#algorithm requires retrieval of longest sustring match with recursive end filling
#this problem leads to false negative matching of sequences with 5' errors below error thesholds
#does not lead to false positives and will perform at least to the level of explicit string matching
# 
#try using difflib for truncation mapping, or alternatively, try ssaha
#allow deletion or insertion lengths of no more than 1nt
#(normally illumina sequencing mismatches or indels will not exceed 1nt)
##################################################
	#idthresh = 1 
	#idsumthresh = 1 #indels are rare in illumina sequencing
	#mmthresh = 1
	#mmsumthresh = 1 #mismatches are more common in illumina sequencing

	s = difflib.SequenceMatcher(None, seq, seq0, autojunk=False)
	a = s.get_opcodes()
	##print a

	idsum = 0
	mmsum = 0
	flag = 0
	ind = 0
	ct = 0
	for sub in a:
		##print sub
		lensubseq = sub[2]-sub[1]
		lensubseq0 = sub[4]-sub[3]

		#move to sequence match start, consider bases removed from front of seq0 as insertion/deletion
		if (flag==0):
			#treat noncontiguous matches of len <= idthesh at start of seqeunce as deletions
			if (hash(sub[0])==hash('equal')) and (lensubseq0 > idthresh):
				##print 'equal and len greater than idthresh'
				flag = 1
				idsum = idsum + sub[3]
				ind = sub[1]
			else:
				##print 'NOT equal and len greater than idthresh'
				continue
				ct = ct+1	

		#treat noncontiguous matches of len <= idthesh at end of seqeunce as deletions
		if sub[4]==(len(seq0)):
			if hash(sub[0])==hash('equal'):
				break
			elif (len(seq0)-a[ct-1][4])>idthresh:
				return -1
			else:
				idsum = idsum+(len(seq0)-a[ct-1][4])				 

		if (abs(idsum) > idsumthresh) or (abs(mmsum) > mmsumthresh):
			##print 'idmmsumthresh exceeded'
			return -1

		if hash(sub[0])==hash('replace'):			
			if max([lensubseq, lensubseq0])<=mmthresh:
				##print 'replace of len less than mmthresh'
				mmsum = mmsum + max([lensubseq, lensubseq0])
			else:
				##print 'NOT replace of len less than mmthresh'
				return -1	
		elif hash(sub[0])==hash('insert'):
			if lensubseq0<=idthresh:
				##print 'insert of len less than idthresh'
				idsum = idsum + lensubseq0
			else:
				##print 'NOT insert of len less than idthresh'
				return -1
		elif hash(sub[0])==hash('delete'):
			if lensubseq<=idthresh:
				##print 'delete of len less than idthresh'
				idsum = idsum - lensubseq
			else:
				##print 'NOT delete of len less than idthresh'
				return -1			

		ct = ct+1
	#return starting index of seq0 within seq:
	#1. individual insertions and deletions of length <= idthresh
	#2. total insertions and deletions compared to template string are <= idsumthresh 
	return ind

