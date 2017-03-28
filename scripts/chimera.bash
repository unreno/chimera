#!/usr/bin/env bash


script=`basename $0`

function usage(){
	echo
	echo "Usage: (NO EQUALS SIGNS)"
	echo
	echo "$script GZIPPED_BAM_SOURCE_FILE_LIST"
	echo
	echo "Example:"
	echo "$script ~/genepi2/ccls/data/raw/TCGA_Glioma_HERV52/TCGA-??-????-10*gz"
	echo
	exit
}

[ $# -eq 0 ] && usage

#	Print a lot more stuff
set -x

{
	echo "Starting at ..."
	date

#	/bin/rm -rf working_dir
	mkdir -p working_dir
	cd working_dir
	working_dir=$PWD

	while [ $# -ne 0 ] ; do
		sample=$1

		echo $sample
		cd $working_dir

		#	sample = /genepi2/ccls/data/raw/TCGA_Glioma_HERV52/TCGA-06-0155-10A-01D-0703-09.HERV52.bam.gz
		gz_link=${sample##*/}      # TCGA-06-0155-10A-01D-0703-09.HERV52.bam.gz
		bam=${gz_link%%.gz}        # TCGA-06-0155-10A-01D-0703-09.HERV52.bam
#		fastq_base=${bam%%.bam}    # TCGA-06-0155-10A-01D-0703-09.HERV52
		sample_base=${gz_link%%.*} # TCGA-06-0155-10A-01D-0703-09

		mkdir -p $sample_base
		cd $sample_base

		ln -s $sample
		gunzip -f $gz_link
		#	will also remove link

		bamToFastq -i $bam -fq $sample_base.1.fastq -fq2 $sample_base.2.fastq

		rm $bam

		chimera_paired_local.bash -v SVAs_and_HERVKs --distance 15 -1 $sample_base.1.fastq -2 $sample_base.2.fastq

		chimera_unpaired_local.bash -v SVAs_and_HERVKs --distance 15 $sample_base.1.fastq,$sample_base.2.fastq

		shift
	done

	echo
	echo "Finished at ..."
	date

} 1>>$script.out 2>&1 &
