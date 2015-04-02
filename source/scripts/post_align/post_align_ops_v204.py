import sys
import csv
import numpy
import pysam
import sys_ops
import clustering_ops
from Bio.Seq import Seq
from operator import itemgetter 

ALIGN_NN_THRESHOLD = 30
ORIENTATION_NN_THRESHOLD = 4
OVERLAP_THRESHOLD = -5 #corresponds to an overlap of 6, negative indicates left strand first then right strand

def write_sams_location_orientation(sam_filename,sampid,outfile_clust_ts,outfile_clust_bs,CLUSTER_READS_THRESHOLD,write_output):

    print "reading input files"
    try:
        samfile = pysam.Samfile(sam_filename)
    except:
        sys_ops.throw_exception("Could not find file "+sam_filename)
        return          

    sam_alignments_ts = []
    sam_alignments_bs = []
    with open(outfile_clust_bs,'w') as csvfile_bs:
        mywriter_bs = csv.writer(csvfile_bs, delimiter=',')   
        with open(outfile_clust_ts,'w') as csvfile_ts:
            mywriter_ts = csv.writer(csvfile_ts, delimiter=',')         
            for read in samfile.fetch():

                #check that read is mapped (mapq>0) and then check read orientation
                if read.mapq>0:
                    if read.is_reverse:
                        #using len(read.seq) is not entirely accurate (should make function to parse cigar string and output mapping length)
                        mywriter_bs.writerow([read.rname,(read.aend-1),(read.aend-1),(read.aend-1),1,1,0,0,1,hash(sampid)])
                    else:
                        mywriter_ts.writerow([read.rname,read.pos,read.pos,read.pos,1,1,0,0,0,hash(sampid)])    
    samfile.close                           

    return 1

def precompile_clusters_sam_nobcs(rawdata_file,write_output):

    rawdata = []
    try:
        with open(rawdata_file, 'rb') as csvfile:
            print 'opened csv'
            csvreader = csv.reader(csvfile, delimiter=',')
            rawdata = [row for row in csvreader]
    except:
        sys_ops.throw_exception("Could not open file " + rawdata_file)
        sys.exit(1) 

    rawdata_sort = sorted(rawdata)
    chrids = map(itemgetter(0), rawdata_sort)

    ct = 0
    unique_chrids = {}
    for chrid in chrids:
        if chrid not in unique_chrids:
            unique_chrids[chrid] = [ct,0]
        else:
            unique_chrids[chrid][1] = ct

        ct+=1

    if write_output:
        print 'writing clusters' 
        for chrid in unique_chrids:
            out_file = rawdata_file.replace('.csv','_'+str(chrid)+'.csv')
            with open(out_file,'w') as csvfile:
                mywriter = csv.writer(csvfile, delimiter=',')
                for row in rawdata_sort[unique_chrids[chrid][0]:unique_chrids[chrid][1]]:
                    mywriter.writerow(row+[ct]) 

def compile_clusters_sam_nobcs(rawdata_file,outfile_clust,CLUSTER_READS_THRESHOLD,write_output):

    rawdata = []
    try:
        with open(rawdata_file, 'rb') as csvfile:
            print 'opened csv'
            csvreader = csv.reader(csvfile, delimiter=',')
            rawdata = [row for row in csvreader]
    except:
        sys_ops.throw_exception("Could not open file " + rawdata_file)
        sys.exit(1)         
  
    rawdata_chrid = map(itemgetter(0), rawdata)
    rawdata_chrlb = map(itemgetter(1), rawdata)
    rawdata_chrub = map(itemgetter(2), rawdata)
    rawdata_chrmed = map(itemgetter(3), rawdata)
    rawdata_readct = map(itemgetter(4), rawdata)
    rawdata_uniquect = map(itemgetter(5), rawdata)
    rawdata_empty1 = map(itemgetter(6), rawdata)
    rawdata_empty2 = map(itemgetter(7), rawdata)
    rawdata_orientid = map(itemgetter(8), rawdata)
    rawdata_sampid = map(itemgetter(9), rawdata)
    rawdata_repreadct = map(itemgetter(10), rawdata)

    ct = 0
    repreadcts = {}
    repreadct0 = 0   
    sam_alignments = []
    for el in rawdata_chrid:
        repreadct = rawdata_repreadct[ct]
        if str(repreadct) not in repreadcts:
            repreadcts[str(repreadct)] = 1
            repreadct0+=int(rawdata_repreadct[ct])
        sam_alignments.append([int(rawdata_chrid[ct]),int(rawdata_chrlb[ct]),int(rawdata_chrub[ct]),ct])   
        ct+=1

    #cluster all reads for an individual ibc by position 
    clusters_loc = clustering_ops.nncluster_chr_positions(sam_alignments,'cbb',CLUSTER_READS_THRESHOLD,ALIGN_NN_THRESHOLD)
    name_clusters = clusters_loc[0]
    alignment_clusters = clusters_loc[1]

    ct = 0
    clusters = []
    #find lbcs corresponding to names in name_clusters
    with open(outfile_clust,'w') as csvfile:
        mywriter = csv.writer(csvfile, delimiter=',')    
        for name_cluster in name_clusters:
            chrid = alignment_clusters[ct][0][0]
            chrlb = min([int(rawdata_chrlb[c]) for c in name_cluster])
            chrub = max([int(rawdata_chrub[c]) for c in name_cluster])
            chrmed = numpy.median(numpy.array([float(rawdata_chrmed[c]) for c in name_cluster]))
            chrlocs = [str([rawdata_chrlb[c],rawdata_chrub[c],rawdata_chrmed[c]]) for c in name_cluster]

            readct = len(name_cluster)             
            uniquect = {}
            for loc in chrlocs:
                if hash(str(loc)) not in uniquect:
                    uniquect[hash(str(loc))] = loc

            reps = {}
            repstr = ''
            locs_ts = []
            locs_bs = []
            locs_ts2 = []
            locs_bs2 = []
            locs_tsu = []
            locs_bsu = []
            locs_tsu2 = []
            locs_bsu2 = []        
            for c in name_cluster:
                sampid = int(rawdata_sampid[c])
                if sampid not in reps:
                    reps[sampid] = 1
                    repstr+=(':'+str(sampid)) 

                if int(rawdata_orientid[c])==0:
                    locs_ts.append(int(rawdata_chrlb[c]))
                elif int(rawdata_orientid[c])==1:
                    locs_bs.append(int(rawdata_chrlb[c]))      

            #remove bias due to PCR amplification bias, take only unique mappings
            locs_tsu = numpy.unique(numpy.array(locs_ts)).tolist()
            locs_bsu = numpy.unique(numpy.array(locs_bs)).tolist()

            if (len(locs_ts)>0 and len(locs_bs)>0):
                locs_ts2 = numpy.matrix([locs_ts]*len(locs_bs))
                locs_bs2 = numpy.matrix([locs_bs]*len(locs_ts))
                locs_diff = locs_ts2-numpy.transpose(locs_bs2)
                locs_g0 = float((locs_diff>=OVERLAP_THRESHOLD).sum())/float((locs_diff.shape[0]*locs_diff.shape[1]))
                locs_l0 = float((locs_diff<OVERLAP_THRESHOLD).sum())/float((locs_diff.shape[0]*locs_diff.shape[1]))

                locs_tsu2 = numpy.matrix([locs_tsu]*len(locs_bsu))
                locs_bsu2 = numpy.matrix([locs_bsu]*len(locs_tsu))
                locs_diffu = locs_tsu2-numpy.transpose(locs_bsu2)
                locs_g0u = float((locs_diffu>=OVERLAP_THRESHOLD).sum())/float((locs_diffu.shape[0]*locs_diffu.shape[1]))
                locs_l0u = float((locs_diffu<OVERLAP_THRESHOLD).sum())/float((locs_diffu.shape[0]*locs_diffu.shape[1]))

                with open(str.replace(outfile_clust,'.csv','')+'_'+str(chrlb)+'_'+str(chrub)+'.csv','w') as csvfile:
                    mywriterc = csv.writer(csvfile, delimiter=',')
                    for row in locs_diff:
                        mywriterc.writerow(row.tolist()[0])    

            else:
                locs_g0 = -1
                locs_l0 = -1
                locs_g0u = -1
                locs_l0u = -1                                    

            orientid = rawdata_orientid[name_cluster[0]]

            if locs_g0u != -1: #if one is zero, all the others are zero
                mywriter.writerow([chrid,chrlb,chrub,chrmed,readct,len(uniquect),0,0,len(reps),repreadct0,0,0,locs_g0u,locs_l0u])

            ct+=1                         

    return 1

