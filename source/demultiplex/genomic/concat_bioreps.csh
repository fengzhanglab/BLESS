#! /bin/csh -f
#

#[[EMX-Sp101],[EMX-Sp102],[EMX-Sa101],[EMX-Sa102],[EMX1.3],[pUC]]
#repids = [['ATCACG','GTGGCC'],['CGATGT','GTTTCG'],['TTAGGC','CGTACG'],['TGACCA','GAGTGG'],['ACAGTG','ACTGAT'],['GCCAAT','ATTCCT','CAGATC','ACTTGA']]

# Use this script to concatenate together the bioreps into a single file
# The original files from which the bioreps came are already hash'd so they are identifiable

mkdir Data/_concat
cat Data/*ATCACG**ts.csv >> Data/_concat/1_2_HAALCADXX.1_2.0.unmapped.R1_2_001_genomic_sam_clusters_ts.csv
cat Data/*GTGGCC**ts.csv >> Data/_concat/1_2_HAALCADXX.1_2.0.unmapped.R1_2_001_genomic_sam_clusters_ts.csv
cat Data/*ATCACG**bs.csv >> Data/_concat/1_2_HAALCADXX.1_2.0.unmapped.R1_2_001_genomic_sam_clusters_bs.csv
cat Data/*GTGGCC**bs.csv >> Data/_concat/1_2_HAALCADXX.1_2.0.unmapped.R1_2_001_genomic_sam_clusters_bs.csv

cat Data/*CGATGT**ts.csv >> Data/_concat/1_2_HAALCADXX.1_2.1.unmapped.R1_2_001_genomic_sam_clusters_ts.csv
cat Data/*GTTTCG**ts.csv >> Data/_concat/1_2_HAALCADXX.1_2.1.unmapped.R1_2_001_genomic_sam_clusters_ts.csv
cat Data/*CGATGT**bs.csv >> Data/_concat/1_2_HAALCADXX.1_2.1.unmapped.R1_2_001_genomic_sam_clusters_bs.csv
cat Data/*GTTTCG**bs.csv >> Data/_concat/1_2_HAALCADXX.1_2.1.unmapped.R1_2_001_genomic_sam_clusters_bs.csv

cat Data/*TTAGGC**ts.csv >> Data/_concat/1_2_HAALCADXX.1_2.2.unmapped.R1_2_001_genomic_sam_clusters_ts.csv
cat Data/*CGTACG**ts.csv >> Data/_concat/1_2_HAALCADXX.1_2.2.unmapped.R1_2_001_genomic_sam_clusters_ts.csv
cat Data/*TTAGGC**bs.csv >> Data/_concat/1_2_HAALCADXX.1_2.2.unmapped.R1_2_001_genomic_sam_clusters_bs.csv
cat Data/*CGTACG**bs.csv >> Data/_concat/1_2_HAALCADXX.1_2.2.unmapped.R1_2_001_genomic_sam_clusters_bs.csv

cat Data/*TGACCA**ts.csv >> Data/_concat/1_2_HAALCADXX.1_2.3.unmapped.R1_2_001_genomic_sam_clusters_ts.csv
cat Data/*GAGTGG**ts.csv >> Data/_concat/1_2_HAALCADXX.1_2.3.unmapped.R1_2_001_genomic_sam_clusters_ts.csv
cat Data/*TGACCA**bs.csv >> Data/_concat/1_2_HAALCADXX.1_2.3.unmapped.R1_2_001_genomic_sam_clusters_bs.csv
cat Data/*GAGTGG**bs.csv >> Data/_concat/1_2_HAALCADXX.1_2.3.unmapped.R1_2_001_genomic_sam_clusters_bs.csv

cat Data/*ACAGTG**ts.csv >> Data/_concat/1_2_HAALCADXX.1_2.4.unmapped.R1_2_001_genomic_sam_clusters_ts.csv
cat Data/*ACTGAT**ts.csv >> Data/_concat/1_2_HAALCADXX.1_2.4.unmapped.R1_2_001_genomic_sam_clusters_ts.csv
cat Data/*ACAGTG**bs.csv >> Data/_concat/1_2_HAALCADXX.1_2.4.unmapped.R1_2_001_genomic_sam_clusters_bs.csv
cat Data/*ACTGAT**bs.csv >> Data/_concat/1_2_HAALCADXX.1_2.4.unmapped.R1_2_001_genomic_sam_clusters_bs.csv

cat Data/*GCCAAT**ts.csv >> Data/_concat/1_2_HAALCADXX.1_2.5.unmapped.R1_2_001_genomic_sam_clusters_ts.csv
cat Data/*ATTCCT**ts.csv >> Data/_concat/1_2_HAALCADXX.1_2.5.unmapped.R1_2_001_genomic_sam_clusters_ts.csv
cat Data/*CAGATC**ts.csv >> Data/_concat/1_2_HAALCADXX.1_2.5.unmapped.R1_2_001_genomic_sam_clusters_ts.csv
cat Data/*ACTTGA**ts.csv >> Data/_concat/1_2_HAALCADXX.1_2.5.unmapped.R1_2_001_genomic_sam_clusters_ts.csv
cat Data/*GCCAAT**bs.csv >> Data/_concat/1_2_HAALCADXX.1_2.5.unmapped.R1_2_001_genomic_sam_clusters_bs.csv
cat Data/*ATTCCT**bs.csv >> Data/_concat/1_2_HAALCADXX.1_2.5.unmapped.R1_2_001_genomic_sam_clusters_bs.csv
cat Data/*CAGATC**bs.csv >> Data/_concat/1_2_HAALCADXX.1_2.5.unmapped.R1_2_001_genomic_sam_clusters_bs.csv
cat Data/*ACTTGA**bs.csv >> Data/_concat/1_2_HAALCADXX.1_2.5.unmapped.R1_2_001_genomic_sam_clusters_bs.csv