#! /bin/csh -f
#

#Directory configuration
set pre_dir = ..
set dir = align
set c_dir = $dir"/concatenate"
set scriptsdir = ../../scripts
set index = /seq/regev_genome_portal/RESOURCES/human/Hg19/Hg19.fa
#set index = /seq/regev_genome_portal/RESOURCES/mouse/mm9/Mus_musculus_assembly9.fasta
set tdf_index = mm9

#bsub configuration
set queue = priority

#bowtie config
# --best (best alignment) -v (n mismatch) -S (sam ouput)
set bowtie = ../../scripts/bowtie-1.1.0/bowtie
set bowtie_params = "--best -v2 -S"

set NSPLIT = 2000000

##################################################
if($argv[1] == "run_samples_prealignment") then 
# RUN BOWTIE2 ALIGNMENTS
	foreach f ( `ls $pre_dir/*R1**_001.fastq` )
		set f2 = `echo $f | sed 's/R1/R2/'`
		set f_root = "$f"_""
		set f2_root = "$f2"_""
		split -l $NSPLIT $f $f_root
		split -l $NSPLIT $f2 $f2_root
		foreach ff ( `ls $f_root*` )
			set ff2 = `echo $ff | sed 's/R1/R2/'`
			python2.7 $scriptsdir"/pre_align/truncate_v200.py" $ff R1
			python2.7 $scriptsdir"/pre_align/truncate_v200.py" $ff2 R2
		end
	end	
endif

##################################################
#DEPRECATED
if($argv[1] == "run_samples_concat_prealignment") then 
# RUN BOWTIE2 ALIGNMENTS
	foreach f ( `ls $pre_dir/*R1**_001.fastq` )
		set f2 = `echo $f | sed 's/R1/R2/'`
		set f_root = "$f"_""
		set f2_root = "$f2"_""
		foreach ff ( `ls $f_root**` )
			set ff0 = `echo $ff | sed 's/.fastq//'`
			set ff0 = $ff0"_genomic.fastq"
			set ff2 = `echo $ff0 | sed 's/genomic.fastq/genomic0.fastq/'`
			set ff3 = `echo $ff0 | sed 's/R1/R2/'`
			mv $ff0 $ff2
			cat $ff2 >> $ff0
			cat $ff3 >> $ff0
		end
	end	
endif

##################################################
if($argv[1] == "run_samples_postprealignment") then 
# RUN BOWTIE2 ALIGNMENTS
	foreach f ( `ls $pre_dir/*_001.fastq` )
		set f_out = `echo $f | sed 's/.fastq/_/'`
		set f_out_gen = $f_out"genomic.fastq"
		foreach ff ( `ls $f_out**genomic.fastq` )
			set ffgen = $ff

			echo ""Writing: "$ffgen" to: "$f_out_gen"
			cat $ffgen >> $f_out_gen
		end							
	end	
endif

##################################################
if($argv[1] == "run_samples_cleanprealignment") then 
# RUN BOWTIE2 ALIGNMENTS
	foreach f ( `ls $pre_dir/*_001.fastq` )
		set f_root = $f"_"
		set f_out = `echo $f | sed 's/.fastq/_/'`
		set f_out_gen = $f_out"genomic.fastq"
		foreach ff ( `ls $f_out**_genomic.fastq` )
			set ffgen = $ff
			#rm $ffgen	
		end
		rm $f_root*							
	end	
endif

##################################################
if($argv[1] == "run_samples_alignment_fwd") then 
# RUN BOWTIE2 ALIGNMENTS
	foreach f ( `ls $pre_dir/*R1**001.fastq` )
		set f_out = `echo $f | sed 's/.fastq/_/'`
		foreach ff ( `ls $f_out**genomic.fastq` )
			set ffsam = `basename $ff`
			set ffsam2 = $dir"/"$ffsam
			set ffsam3 = `echo $ffsam2 | sed 's/.fastq/.sam/'`
			$bowtie $bowtie_params $index $ff $ffsam3
		end					
	end	
endif

##################################################
if($argv[1] == "run_samples_alignment_rev") then 
# RUN BOWTIE2 ALIGNMENTS
	foreach f ( `ls $pre_dir/*R2**001.fastq` )
		set f_out = `echo $f | sed 's/.fastq/_/'`
		foreach ff ( `ls $f_out**genomic.fastq` )
			set ffsam = `basename $ff`
			set ffsam2 = $dir"/"$ffsam
			set ffsam3 = `echo $ffsam2 | sed 's/.fastq/.sam/'`
			$bowtie $bowtie_params $index $ff $ffsam3		
		end					
	end	
endif

##################################################
if($argv[1] == "run_samples_alignment_paired_end") then 
# RUN BOWTIE2 ALIGNMENTS
	foreach f ( `ls $pre_dir/*R2**001.fastq` )
		set f_out = `echo $f | sed 's/.fastq/_/'`
		foreach ff ( `ls $f_out**genomic.fastq` )
			set ff2 = `echo $f | sed 's/R1/R2/'`
			set ffsam = `basename $ff`
			set ffsam2 = $dir"/"$ffsam
			set ffsam3 = `echo $ffsam2 | sed 's/.fastq/.sam/'`
			$bowtie $bowtie_params $index -1 $ff -2 $ff2 $ffsam3
		end					
	end		
endif

##################################################
if($argv[1] == "run_samples_postalignment") then 
# RUN BOWTIE2 ALIGNMENTS
	foreach f ( `ls $dir/*.sam` )
		python2.7 $scriptsdir"/post_align/main_v204.py" sample_write_sam_location_orientation $f		
	end	
endif

##################################################
if($argv[1] == "run_samples_precompilealignment") then 
# RUN BOWTIE2 ALIGNMENTS
	foreach f ( `ls $c_dir/**.csv` )
		python2.7 $scriptsdir"/post_align/main_v204.py" precompile_clusters_sam_nobcs $f		
	end	
endif

################################################## (from DAS)
if($argv[1] == "run_samples_compilealignment") then 
# RUN BOWTIE2 ALIGNMENTS
	foreach f ( `ls $c_dir/**_ts_**.csv` )
		set f2 = `echo $f | sed 's/_ts_/_bs_/'`
		set f3 = `echo $f | sed 's/_ts_/_/'`
		cat $f >> $f3
		cat $f2 >> $f3
		python2.7 $scriptsdir"/post_align/main_v204.py" compile_clusters_sam_nobcs $f3
	end	
endif

################################################################################
################################################################################

##################################################
if($argv[1] == "run_samples_sam_2_bam") then 
# RUN BOWTIE2 ALIGNMENTS
	foreach f ( `ls $dir/*.sam` )
		set out = `echo $f | sed 's/.sam/.bam/'`
		samtools view -bS -o $out $f
	end	
endif

##################################################
if($argv[1] == "run_samples_sort_bam") then 
# RUN BOWTIE2 ALIGNMENTS
	foreach f ( `ls $dir/*.bam` )
		set out = `echo $f | sed 's/.bam/.sorted/'`	
		samtools sort $f $out
	end	
endif

##################################################
if($argv[1] == "run_samples_index_sortedbam") then 
# RUN BOWTIE2 ALIGNMENTS
	foreach f ( `ls $dir/*sorted.bam` )
		set out = `echo $f | sed 's/.sorted.bam/.sorted.bam.bai/'`	
		samtools index $f $out
	end	
endif

##################################################
if($argv[1] == "run_samples_sortedbam_2_tdf") then 
# RUN BOWTIE2 ALIGNMENTS
	foreach f ( `ls $dir/*.sorted.bam` )
		set out = `echo $f | sed 's/.sorted.bam/.tdf/'`	
		"$scriptsdir"/IGVTools/igvtools"" count $f $out $tdf_index
	end	
endif
