import os
import csv
import sys
import sys_ops
import job_ops
import clustering_ops
import post_align_ops_v204

SAMPLE_CLUSTER_READS_THRESHOLD = 1
REPLICATE_CLUSTER_READS_THRESHOLD = 1
JUNCTION_CLUSTER_READS_THRESHOLD = 1

IBC_THRESHOLD = 1
IBC_EXACT_MATCH = 0
IBC_READS_THRESHOLD = 1
LBC_THRESHOLD = 1
LBC_EXACT_MATCH = 0
LBC_READS_THRESHOLD = 1

suffix_sam = '_genomic.sam'
sam_subpath = 'Data/'

suffix_ibc = '_ibc.fastq'
suffix_lbc = '_lbc.fastq'
ibclbc_subpath = '../'

write_ibc_output = 0
write_sam_output = 0
write_cluster_output = 1

def main():
    print "starting main"

    op_request = sys.argv[1]
    filename = sys.argv[2]                                                        

    if op_request == "sample_write_sam_location_orientation":
        print "sample_write_sam_location_orientation"

        run_sample_write_sam_location_orientation(filename)

    if op_request == "precompile_clusters_sam_nobcs":
        print "precompile_clusters_sam_nobcs"

        precompile_clusters_sam_nobcs(filename) 

    if op_request == "compile_clusters_sam_nobcs":
        print "compile_clusters_sam_nobcs"

        compile_clusters_sam_nobcs(filename) 

def run_sample_write_sam_location_orientation(sam_filename):
    #cluster sams, then ibcs for sam clusts, then lbcs for ibc clusts

    sam_filename2 = sam_filename.replace(sam_subpath,ibclbc_subpath)

    outfile_ts = sam_filename.replace('.sam','_sam_clusters_ts.csv')
    outfile_bs = sam_filename.replace('.sam','_sam_clusters_bs.csv')
    
    ############################################################
    sampid = sam_filename.split('.')[0]
    print sampid
    ############################################################

    finished = post_align_ops_v204.write_sams_location_orientation(sam_filename,sampid,outfile_ts,outfile_bs,SAMPLE_CLUSTER_READS_THRESHOLD,write_cluster_output)
    if (os.path.isfile(outfile_ts)==False) or (os.path.isfile(outfile_bs)==False):
        sys_ops.throw_exception("Failed at cluster_sams on processing " + outfile)
        return

def precompile_clusters_sam_nobcs(rawdata_filename):
    #assemble output for parallel clustering jobs
    
    finished = post_align_ops_v204.precompile_clusters_sam_nobcs(rawdata_filename,write_cluster_output)
        
def compile_clusters_sam_nobcs(rawdata_filename):
    #assemble output for parallel clustering jobs

    outfile = rawdata_filename.replace('.csv','c.csv')
    
    finished = post_align_ops_v204.compile_clusters_sam_nobcs(rawdata_filename,outfile,SAMPLE_CLUSTER_READS_THRESHOLD,write_cluster_output)
    if os.path.isfile(outfile)==False:
        sys_ops.throw_exception("Failed at compile_clusters_sam on processing " + outfile)
        return

main() 
