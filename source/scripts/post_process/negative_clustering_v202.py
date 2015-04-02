
# coding: utf-8

# In[1]:

import os
import csv
import sys
import numpy
from Bio.Seq import Seq


def main(file_name_1,file_name_2,file_out):

    # In[2]:

    print 'importing data'

    #Import the sample data and concatenate
    repids =[]
    raw_data = []
    with open(file_name_1, 'rb') as csvfile:
        csvreader = csv.reader(csvfile, delimiter=',')
        ct = 0
        for row in csvreader:
            raw_data+=[row]
            ct+=1
        repids+=[0]*ct #makes 1st list have repid = 0

    with open(file_name_2, 'rb') as csvfile:
        csvreader = csv.reader(csvfile, delimiter=',')
        ct = 0
        for row in csvreader:
            raw_data+=[row]
            ct+=1
        repids+=[1]*ct #makes 2nd list have repid = 1

    print len(raw_data)
    print len(repids)


    # In[15]:

    print 'negative clustering'

    ct = 0
    alignment_clusters = []
    for row in raw_data:
        alignment_clusters.append([int(row[0]),int(row[1]),int(row[2]),ct])
        ct+=1
    
    print 'clustering'

    #cluster all reads by position across an individual replicate
    clusters_loc = nncluster_chr_positions(alignment_clusters,'cbb',1,0)
    name_clusters2 = clusters_loc[0]
    alignment_clusters2 = clusters_loc[1]
    
    print 'compiling clustering data'

    OPP_DATA = compile_nncluster_rep2negct(raw_data,repids,name_clusters2,alignment_clusters2,[0,1])

    readct_neg = OPP_DATA[0]
    repreadct_neg = OPP_DATA[1]


    # In[16]:

    print 'writing sequences' 
    with open(file_out,'w') as csvfile:
        mywriter = csv.writer(csvfile, delimiter=',')
        ct = 0
        for row in raw_data:
            #do not write out negative control data
            if repids[ct]==1:
                break
            mywriter.writerow(row+[repids[ct]]+[readct_neg[ct]]+[repreadct_neg[ct]])
            ct+=1

            
# In[11]:

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
            upperbound = -1
            for position in sorted_positions:
                if ct==0:
                    name_clusters.append([position[2]])
                    alignment_clusters.append([[chrid, position[0], position[1]]])
                else:

                    #modification 141218
                    #if a preceding cluster is very large, check that the upperbound on the last cluster is actually greater
                    #otherwise, keep preceding upperbound
                    if sorted_positions[ct-1][1]>upperbound:
                        upperbound = sorted_positions[ct-1][1]

                    #modification 140829
                    #sorted defaults to sorting on the lower bound of position pair
                    #since sorting by ascending size, position[0] and [1] should be larger than sorted_positions[0]
                    #thf, in addition to proximity check for all bounds, check if positions are overlapping but neither LB or UB for the pairs are within threshold
                    #given conditions above, if (position[0]-sorted_positions[ct-1][1])<=0, then bounds are overlapping
                    temp = [abs(position[0]-sorted_positions[ct-1][0]), abs(position[0]-sorted_positions[ct-1][1]), abs(position[1]-sorted_positions[ct-1][0]), abs(position[1]-sorted_positions[ct-1][1])]
                    if (min(temp)<=ALIGN_NN_THRESHOLD) or ((position[0]-upperbound)<=0):
                        name_clusters[len(alignment_clusters)-1].append(position[2])
                        alignment_clusters[len(alignment_clusters)-1].append([chrid, position[0], position[1]])
                    elif (min(temp)>ALIGN_NN_THRESHOLD):
                        name_clusters.append([position[2]])
                        alignment_clusters.append([[chrid, position[0], position[1]]])
                    
                ct = ct+1                

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



# In[12]:

def compile_nncluster_rep2negct(data_in,repids,name_clusters,alignment_clusters,repid_pair):
    opp_repreadcts = [0]*len(data_in)
    opp_readcts = [0]*len(data_in)

    ct = 0
    for name_cluster in name_clusters:

        OUT = []

        readcts = {}
        repreadcts = {}
        #sampid1 uniquely specifies the sample
        for ind in name_cluster:
            repid = str(repids[ind])
            if repid not in repreadcts:
                readcts[repid] = int(data_in[ind][4])
                repreadcts[repid] = int(data_in[ind][9])
            else:
                readcts[repid]+=int(data_in[ind][4])
        
        for ind in name_cluster:
            repid = str(repids[ind])
            if repid==str(repid_pair[0]):
                if str(repid_pair[1]) in repreadcts:
                    opp_repreadcts[ind] = repreadcts[str(repid_pair[1])]
                    opp_readcts[ind] = readcts[str(repid_pair[1])]
            if repid==str(repid_pair[1]):
                if str(repid_pair[0]) in repreadcts:
                    opp_repreadcts[ind] = repreadcts[str(repid_pair[0])]
                    opp_readcts[ind] = readcts[str(repid_pair[0])]           

        ct+=1

    return [opp_readcts,opp_repreadcts]

file_name_1 = sys.argv[1]
file_name_2 = sys.argv[2]
file_out = sys.argv[3]
main(file_name_1,file_name_2,file_out)


