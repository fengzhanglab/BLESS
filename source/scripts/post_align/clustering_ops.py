import sys
import csv
import numpy
import pysam
import pickle
import difflib
import sys_ops
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.SeqRecord import SeqRecord
from Bio.SubsMat import MatrixInfo as matlist
from itertools import izip_longest

def fastq_get_seq_dict(fastqfile):
    try:
        fq_handle = open(fastqfile, "rU")
    except:
        sys_ops.throw_exception("Could not find file "+fastqfile)
        return  

    ct = 0
    READS = {}
    for record in SeqIO.parse(fq_handle, "fastq"):
        READS[hash(record.name)] = str(record.seq)
    
    fq_handle.close()

    return READS

def fastq_get_unique_bcs(bc_fastqfile,BC_READS_THRESHOLD):
    #DEPRECATED
    try:
        print "fastq_get_unique_bcs, opening "+bc_fastqfile
        bc_handle = open(bc_fastqfile, "rU")
    except:
        sys_ops.throw_exception("Could not find file "+bc_fastqfile)
        return      
    
    #read all unique bcs into list and record counts for each bc
    ct = 0
    BCSEQUNIQUE = {}
    for record in SeqIO.parse(bc_handle, "fastq"):
        bc = str(record.seq)
        if (hash(bc) not in BCSEQUNIQUE):
            BCSEQUNIQUE[hash(bc)] = [1, bc]
            ct = ct+1
        else:
            BCSEQUNIQUE[hash(bc)][0] = BCSEQUNIQUE[hash(bc)][0]+1

    bc_handle.close()

    BCSEQ = []
    for bchash in BCSEQUNIQUE:
        bc = BCSEQUNIQUE[bchash]
        if bc[0]>=BC_READS_THRESHOLD:
            BCSEQ.append(BCSEQUNIQUE[bchash])

    print "number of reads containing bc for clustering:"
    print ct

    return BCSEQ

def fastq_get_unique_bcs_inds(bc_fastqfile,BC_READS_THRESHOLD):
    try:
        print "fastq_get_unique_bcs, opening "+bc_fastqfile
        bc_handle = open(bc_fastqfile, "rU")
    except:
        sys_ops.throw_exception("Could not find file "+bc_fastqfile)
        return      
    
    #read all unique bcs into list and record counts for each bc
    ct = 0
    BCSEQUNIQUE = {}
    for record in SeqIO.parse(bc_handle, "fastq"):
        bc = str(record.seq)
        #sort bcs into dict by sequence
        if (hash(bc) not in BCSEQUNIQUE):
            BCSEQUNIQUE[hash(bc)] = [[], bc]
            BCSEQUNIQUE[hash(bc)][0].append(ct)
        else:
            BCSEQUNIQUE[hash(bc)][0].append(ct)
        ct+=1

    bc_handle.close()

    BCSEQ = []
    for bchash in BCSEQUNIQUE:
        bc = BCSEQUNIQUE[bchash]
        if bc[0]>=BC_READS_THRESHOLD:
            BCSEQ.append(BCSEQUNIQUE[bchash])

    print "number of reads containing bc for clustering:"
    print ct

    return BCSEQ

def list_get_unique_bcs(bc_seqs,BC_READS_THRESHOLD):   
    #DEPRECATED
    #read all unique bcs into list and record counts for each bc
    ct = 0
    BCSEQUNIQUE = {}
    for bc in bc_seqs:
        #sort bcs into dict by sequence
        if (hash(bc) not in BCSEQUNIQUE):
            BCSEQUNIQUE[hash(bc)] = [1, bc]
            ct = ct+1
        else:
            BCSEQUNIQUE[hash(bc)][0] = BCSEQUNIQUE[hash(bc)][0]+1

    BCSEQ = []
    for bchash in BCSEQUNIQUE:
        bc = BCSEQUNIQUE[bchash]
        if bc[0]>=BC_READS_THRESHOLD:
            BCSEQ.append(BCSEQUNIQUE[bchash])

    #print "number of reads containing bc for clustering:"
    #print ct

    return BCSEQ

def list_get_unique_bcs_inds(bc_seqs,BC_READS_THRESHOLD):   
    
    #read all unique bcs into list and record counts for each bc
    ct = 0
    BCSEQUNIQUE = {}
    for bc in bc_seqs:
        #sort bcs into dict by sequence
        if (hash(bc) not in BCSEQUNIQUE):
            BCSEQUNIQUE[hash(bc)] = [[], bc]
            BCSEQUNIQUE[hash(bc)][0].append(ct)
        else:
            BCSEQUNIQUE[hash(bc)][0].append(ct)
        ct+=1

    BCSEQ = []
    for bchash in BCSEQUNIQUE:
        bc = BCSEQUNIQUE[bchash]
        if len(bc[0])>=BC_READS_THRESHOLD:
            BCSEQ.append(BCSEQUNIQUE[bchash])

    #print "number of reads containing bc for clustering:"
    #print ct

    return BCSEQ

def cluster_bcs(BCSEQ,BC_THRESHOLD,BC_READS_THRESHOLD,BC_EXACT_MATCH,UMI_CLUSTER_METHOD):

    #aggregating all NSWMT sequences by sequence, just those that cluster in UMI space
    
    if str.find(UMI_CLUSTER_METHOD,'explicit')>=0:
        BCSEQC = threshold_cluster_uid_explicit(BCSEQ,BC_THRESHOLD)
    elif str.find(UMI_CLUSTER_METHOD,'prelinked')>=0:
        TEMP = threshold_cluster_uid_prelinked_setup(BCSEQ,BC_THRESHOLD)
        BCSEQC = threshold_cluster_uid_prelinked(TEMP,BC_THRESHOLD)
    else:
        sys_ops.throw_exception("Options for cluster_bcs must be either 'explicit' or 'prelinked'. Exiting...")
        return           

    #print 'number of unique BC clusters before filtering:'
    #print len(BCSEQC)

    #filter our all clusters with only n reads less than read threshold
    BCSEQC2 = []
    #iterate through all clusters and find clusters
    for bcs in BCSEQC:
        #if greater than 1 barcode in the cluser or greater than 1 read in a sincle bc cluster, cluster is not junk
        if BC_EXACT_MATCH==1:
            bcs2 = [bcs[0]]
        else:
            bcs2 = bcs

        ct = 0
        for bc in bcs2:
            ct = ct + bc[1]
        
            if (ct>=BC_READS_THRESHOLD):
                BCSEQC2.append(bcs2)
                break

    #print 'number of unique BC clusters after filtering:'
    #print len(BCSEQC2)

    return BCSEQC2 

def threshold_uid_match(uid0,uid,threshold):
    ct = 0
    mct = 0
    mmct = 0
    #for characters in uid compare string
    for c in uid0:
        #if character at position n matches template uid at position n, count match
        #if not, count as mismatch
        if c==uid[ct]:
            mct = mct+1
        else:        
            mmct = mmct+1

        #during each iteration, check mismatch ct
        #if mm ct is greater than threshold for match, return 0
        if mmct>threshold:
            return 0
        ct = ct+1

    #if iterate through characters in compare and template, and not return 0
    #string is within theshold, return 1
    return 1    

def threshold_cluster_uid_prelinked_setup(uid_list0,threshold):
    #build linkage_file for input into prelinked clustering -- is a python list-of-lists. 
    #uid_list0 is a dictionary of hash codes for all unique uids containing for each uid:
    #       1. a list of indicies of all occurances of the unique uid in originating list of reads
    #       2. the uid sequence
    #Each list-element of this list is a UID sequence that 1 or more reads IDENTICALLY CORRESPOND TO
    #It takes of form as follows uid_list = [[A0,B0,C0,[D0_0,D0_1,D0_2,...]],[A1,B1,C1,[D1_0,D1_1,D1_2,...]],...]
    #A0 is the SEQUENCE ITSELF (variable-type str)
    #B0 is number of reads that possess the sequence A0
    #C0 is the number of reads belonging to all UID sequences that are BOTH:
    #       1. Within a distance "threshold" of sequence A0 and
    #       2. Possessing read-abundances LESS THAN OR EQUAL TO B0
    #[D0_0,D0_1,D0_2,...] is a list of list-indices belonging to the sequences summed across for C0 (fitting the same criteria). 
    #THE FIRST OF THESE INDICES MUST CORRESPOND TO SEQUENCE A0 OUTPUTS FILE WITH EACH LINE: 
    #A_B_C where C is the cluster-index belonged to by unique UID sequence A, itself corresponding to B reads

    #prepare data structure for bc clustering
    ct = 0
    uid_list = []
    for el in uid_list0:
        uid = el[1]
        uidct = len(el[0])
        A0 = uid
        B0 = uidct

        #first index in D should correspond to matching sequence
        D0 = []
        D0.append(ct)

        ct2 = 0
        #iterate through all barcodes for barcodes to find barcodes within theshold mismatch 
        #could speed for low n thesholds by finding theshold mm space for a barcode and hashing over that space
        for el2 in uid_list0:

            uid2 = el2[1]
            uid2ct = len(el2[0])

            #diagonal has already been appended to beginning of D0
            #count index and continue
            if hash(uid)==hash(uid2):
                ct2 = ct2+1
                continue

            #C0 is the number of reads belonging to all UID sequences that are BOTH:
                #       1. Within a distance "threshold" of sequence A0 and
                #       2. Possessing read-abundances LESS THAN OR EQUAL TO B0   
            if uid2ct<=B0:
                #if within mm theshold, append to D and record barcode index in uid_list_rsorted
                if threshold_uid_match(uid,uid2,threshold)==1:
                    D0.append(ct2)

            ct2 = ct2+1

        C0 = 0
        for ind in D0:
            C0+=len(uid_list0[ind][0])


        uid_list.append([A0, B0, C0, D0])

        ct = ct+1

    return uid_list

def threshold_cluster_uid_prelinked(uid_list,threshold):
    #RECOMMENDED: threshold should be set = 1
    #P, subsample, prefix may be left as default
    #identical_uid0_file is a filename (for output purposes) that may be set as an empty string

    #load linkage_file into uid_list -- is a python list-of-lists.
    #Each list-element of this list is a UID sequence that 1 or more reads IDENTICALLY CORRESPOND TO
    #It takes of form as follows uid_list = [[A0,B0,C0,[D0_0,D0_1,D0_2,...]],[A1,B1,C1,[D1_0,D1_1,D1_2,...]],...]
    #A0 is the SEQUENCE ITSELF (variable-type str)
    #B0 is number of reads that possess the sequence A0
    #C0 is the number of reads belonging to all UID sequences that are BOTH:
    #       1. Within a distance "threshold" of sequence A0 and
    #       2. Possessing read-abundances LESS THAN OR EQUAL TO B0
    #[D0_0,D0_1,D0_2,...] is a list of list-indices belonging to the sequences summed across for C0 (fitting the same criteria). 
    #THE FIRST OF THESE INDICES MUST CORRESPOND TO SEQUENCE A0 OUTPUTS FILE WITH EACH LINE: 
    #A_B_C where C is the cluster-index belonged to by unique UID sequence A, itself corresponding to B reads

    #sort uid_list by decreasing read-number
    num_uid = len(uid_list)
    sorted_uid_list = sorted(uid_list, key=lambda row: -row[2]) #note: sorted_uid_list _REMAINS_ a pointer to uid_list
    index_vals = [-1 for i in range(num_uid)] 
    #print "assigning clusters ..."
        
    for sorted_uid_el in sorted_uid_list: 
        #index_vals, with indices corresponding to _original_ positions in pre-sorted uid_list, are initiated at -1 (stored in list at row[3])
        #UID's accepted into cluster with seed of index i, will be given value i in index_vals
        #UID's rejected from all classification are given index
        if index_vals[sorted_uid_el[3][0]] < 0: #if this seed has index -1 (has not been assigned to any seed itself)
            index_vals[sorted_uid_el[3][0]] = int(sorted_uid_el[3][0]) # set cluster seed to itself
            
        my_index_val = int(index_vals[sorted_uid_el[3][0]])
        
        for i in range(1,len(sorted_uid_el[3])):
            if index_vals[sorted_uid_el[3][i]] < 0: #connected read is unassigned -- assign to current cluster seed
                index_vals[sorted_uid_el[3][i]] = my_index_val        

    new_uid_list = []
    index_vals = [int(x) for x in index_vals]

    #consolidate clustered UID's   
    for i in range(num_uid):
        if i in index_vals:
            my_uid_list = [[uid_list[j][0], uid_list[j][1], uid_list[j][3][0]] for j in range(num_uid) if index_vals[j]==i]
            new_uid_list.append(sorted(my_uid_list, key=lambda row: -row[1]))       
    
    return new_uid_list    

def threshold_cluster_uid_explicit(uid_list,threshold):
    #prepare data structure for bc clustering
    
    uid_list_rsorted = sorted(uid_list,reverse=True)

    ct = 0
    uid_list_clusters = []
    ASSIGNED = [0]*len(uid_list_rsorted)
    # iterate through all barcode sequences
    for el in uid_list_rsorted:
        uid = el[1]

        #once all uids have been assigned, clustering is finished
        if (sum(ASSIGNED)==len(uid_list_rsorted)):
            break
        #if uid has been assigned, do not attempt to build cluster around uid
        elif ASSIGNED[ct]==1:
            ct = ct+1
            continue

        ct2 = 0
        uid_ind = []
        uid_list_cluster = []

        #first index in D should correspond to matching sequence
        uid_ind.append(ct)

        #iterate through all barcodes for barcodes to find barcodes within theshold mismatch 
        #could speed for low n thesholds by finding theshold mm space for a barcode and hashing over that space
        for el2 in uid_list_rsorted:
            uid2 = el2[1]
            #if matching, just count and continue
            #need to catch matches before they are incorporated into D
            if hash(uid)==hash(uid2):
                ct2 = ct2+1
                continue
            elif ASSIGNED[ct2]==1:
                ct2 = ct2+1
                continue                
            #if within mm theshold, append to D and record barcode index in uid_list_rsorted
            elif threshold_uid_match(uid,uid2,threshold)==1:
                uid_ind.append(ct2)
            ct2 = ct2+1

        for ind in uid_ind:
            uid_list_cluster.append([uid_list_rsorted[ind][1], uid_list_rsorted[ind][0], ind])
            ASSIGNED[ind] = 1 

        uid_list_clusters.append(uid_list_cluster)

        ct = ct+1

    return uid_list_clusters

def nncluster_chr_positions(alignments,input_format,READS_THRESHOLD,ALIGN_NN_THRESHOLD):
#input_format: 'cb'([chr, bound, id]) or 'cbb' (chr, lowerbound, upperbound, id)
    #sort by chromosome
    
    name_clusters = [] 
    alignment_clusters = [] 
    alignments_sort_chr = {}
    if hash(input_format)==hash('cb'):
        for alignment in alignments:
            if alignment[0] not in alignments_sort_chr:
                alignments_sort_chr[alignment[0]] = []
                alignments_sort_chr[alignment[0]].append([alignment[1], alignment[2]])              
            else:
                alignments_sort_chr[alignment[0]].append([alignment[1], alignment[2]])
   
        for chrid in alignments_sort_chr:
            positions = alignments_sort_chr[chrid]
            #sort by chromosomal position            
            ct = 0
            sorted_positions = sorted(positions) #note: sorted_alignments _REMAINS_ a pointer to alignments
            for position in sorted_positions:
                if ct==0:
                    name_clusters.append([position[1]])
                    alignment_clusters.append([[chrid, position[0]]])
                elif (position[0]-sorted_positions[ct-1][0])<=ALIGN_NN_THRESHOLD:
                    name_clusters[len(alignment_clusters)-1].append(position[1])
                    alignment_clusters[len(alignment_clusters)-1].append([chrid, position[0]])
                elif (position[0]-sorted_positions[ct-1][0])>ALIGN_NN_THRESHOLD:
                    name_clusters.append([position[1]])
                    alignment_clusters.append([[chrid, position[0]]])

                ct = ct+1

    elif hash(input_format)==hash('cbb'):
        for alignment in alignments:
            if alignment[0] not in alignments_sort_chr:
                alignments_sort_chr[alignment[0]] = []
                alignments_sort_chr[alignment[0]].append([alignment[1], alignment[2], alignment[3]])               
            else:
                alignments_sort_chr[alignment[0]].append([alignment[1], alignment[2], alignment[3]])
   
        for chrid in alignments_sort_chr:
            positions = alignments_sort_chr[chrid]
            #sort by chromosomal position            
            ct = 0
            sorted_positions = sorted(positions) #note: sorted_alignments _REMAINS_ a pointer to alignments

            #print '####################'
            #print 'sorted positions'
            #print sorted_positions
            #print '####################'
            
            for position in sorted_positions:
                if ct==0:
                    name_clusters.append([position[2]])
                    alignment_clusters.append([[chrid, position[0], position[1]]])
                else:
                    #modification 082914
                    #sorted defaults to sorting on the lower bound of position pair
                    #since sorting by ascending size, position[0] and [1] should be larger than sorted_positions[0]
                    #thf, in addition to proximity check for all bounds, check if positions are overlapping but neither LB or UB for the pairs are within threshold
                    #given conditions above, if (position[0]-sorted_positions[ct-1][1])<=0, then bounds are overlapping
                    temp = [abs(position[0]-sorted_positions[ct-1][0]), abs(position[0]-sorted_positions[ct-1][1]), abs(position[1]-sorted_positions[ct-1][0]), abs(position[1]-sorted_positions[ct-1][1])]
                    if (min(temp)<=ALIGN_NN_THRESHOLD) or ((position[0]-sorted_positions[ct-1][1])<=0):
                        name_clusters[len(alignment_clusters)-1].append(position[2])
                        alignment_clusters[len(alignment_clusters)-1].append([chrid, position[0], position[1]])
                    elif (min(temp)>ALIGN_NN_THRESHOLD):
                        name_clusters.append([position[2]])
                        alignment_clusters.append([[chrid, position[0], position[1]]])
                    
                ct = ct+1                

    #print '####################'
    #print 'name_clusters'
    #print name_clusters
    #print 'alignment_clusters'
    #print alignment_clusters
    #print '####################'

    #filter out clusters with reads less than ALIGN_CLUSTER_READS_THRESHOLD
    ct = 0
    name_clusters2 = []
    alignment_clusters2 = []
    for cluster in alignment_clusters:
        #first item in cluster list is ibc hashcode
        if len(cluster)>=READS_THRESHOLD:   
            alignment_clusters2.append(cluster)
            name_clusters2.append(name_clusters[ct])
        
        ct = ct + 1

    return(name_clusters2, alignment_clusters2) 

def nncluster_chr_positions_orientations(alignments,input_format,READS_THRESHOLD,ALIGN_NN_THRESHOLD):
#input_format: 'cb'([chr, bound, id]) or 'cbb' (chr, lowerbound, upperbound, id)
#all clustering output consists of pairs of two input loci with opposite orientation pairing 
#clusters do not have to have a consistent directionality

    name_clusters = [] 
    alignment_clusters = [] 
    alignments_sort_chr = {}
    if hash(input_format)==hash('cb'):
        for alignment in alignments:
            if alignment[0] not in alignments_sort_chr:
                alignments_sort_chr[alignment[0]] = []
                alignments_sort_chr[alignment[0]].append([alignment[1], alignment[2], alignment[3]])              
            else:
                alignments_sort_chr[alignment[0]].append([alignment[1], alignment[2], alignment[3]])
   
        for chrid in alignments_sort_chr:
            positions = alignments_sort_chr[chrid]
            #sort by chromosomal position            
            ct = 0
            orient0 = 0
            sorted_positions = sorted(positions)
            for position in sorted_positions:
                #find first element in cluster with the orientation of the first element in a pair
                if ct==0:
                    name_clusters.append([position[2]])
                    alignment_clusters.append([[chrid, position[0], position[1]]])
                    orient0 = position[1]
                elif ((not (position[1]==orient0)) and (abs(position[0]-sorted_positions[ct-1][0])<=ALIGN_NN_THRESHOLD)):
                    
                    #print '##########'
                    #print position[1]
                    #print orient0
                    #print '##########'
                    
                    #check for accompanying second element of pair
                    name_clusters[len(alignment_clusters)-1].append(position[2])
                    alignment_clusters[len(alignment_clusters)-1].append([chrid, position[0], position[1]])
                    orient0 = position[1]
                else:
                    name_clusters.append([position[2]])
                    alignment_clusters.append([[chrid, position[0], position[1]]])
                    orient0 = position[1]

                ct+=1     

    elif hash(input_format)==hash('cbb'):
        for alignment in alignments:
            if alignment[0] not in alignments_sort_chr:
                alignments_sort_chr[alignment[0]] = []
                alignments_sort_chr[alignment[0]].append([alignment[1], alignment[2], alignment[3], alignment[4]])               
            else:
                alignments_sort_chr[alignment[0]].append([alignment[1], alignment[2], alignment[3], alignment[4]])
   
        for chrid in alignments_sort_chr:
            positions = alignments_sort_chr[chrid]
            #sort by chromosomal position            
            ct = 0
            orient0 = 0
            sorted_positions = sorted(positions)
            for position in sorted_positions:

                if ct==0:
                    name_clusters.append([position[3]])
                    alignment_clusters.append([[chrid, position[0], position[1], position[2]]])
                    orient0 = position[2]
                else:
                    #modification 082914
                    #sorted defaults to sorting on the lower bound of position pair
                    #since sorting by ascending size, position[0] and [1] should be larger than sorted_positions[0]
                    #thf, in addition to proximity check for all bounds, check if positions are overlapping but neither LB or UB for the pairs are within threshold
                    #given conditions above, if (position[0]-sorted_positions[ct-1][1])<=0, then bounds are overlapping
                    temp = [abs(position[0]-sorted_positions[ct-1][0]), abs(position[0]-sorted_positions[ct-1][1]), abs(position[1]-sorted_positions[ct-1][0]), abs(position[1]-sorted_positions[ct-1][1])]
                    #DO NOT CHECK IF ARBITRARILY OVERLAPPING, WANT WITHIN THRESHOLD MATCHING OF ENDS OF CLUSTER
                    #if (not (position[1]==orient0)) and ((min(temp)<=ALIGN_NN_THRESHOLD) or ((position[0]-sorted_positions[ct-1][1])<=0)):
                    if (not (position[2]==orient0)) and (min(temp)<=ALIGN_NN_THRESHOLD):

                        #print '##########'
                        #print position[2]
                        #print orient0
                        #print '##########'

                        name_clusters[len(alignment_clusters)-1].append(position[3])
                        alignment_clusters[len(alignment_clusters)-1].append([chrid, position[0], position[1], position[2]])
                        orient0 = position[2]
                    else:
                        name_clusters.append([position[3]])
                        alignment_clusters.append([[chrid, position[0], position[1], position[2]]])
                        orient0 = position[2]

                ct+=1     

    ##############################
    for el in alignment_clusters:
        if len(el)>1:
            print el
    ##############################

    #filter out clusters with reads less than ALIGN_CLUSTER_READS_THRESHOLD
    ct = 0
    name_clusters2 = []
    alignment_clusters2 = []
    for cluster in alignment_clusters:
        #first item in cluster list is ibc hashcode
        if len(cluster)>=READS_THRESHOLD:   
            alignment_clusters2.append(cluster)
            name_clusters2.append(name_clusters[ct])
        
        ct = ct + 1

    return(name_clusters2, alignment_clusters2) 


