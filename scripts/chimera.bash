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

		sample_base=${sample##*/}
		sample_base=${sample_base%%.*}

		mkdir -p $sample_base
		cd $sample_base

		ln -s $sample
		gunzip -f $sample	
		#	will also remove link

		bam=${sample%%.gz}
		echo $bam

		fastq_base=${bam%%.bam}
		echo $fastq_base

		bamToFastq -i $bam -fq $fastq_base.1.fastq -fq2 $fastq_base.2.fastq

		chimera_paired_local.bash -v SVAs_and_HERVKs -1 $fastq_base.1.fastq -2 $fastq_base.2.fastq

		chimera_unpaired_local.bash -v SVAs_and_HERVKs $fastq_base.1.fastq,$fastq_base.2.fastq

		shift
	done

	echo
	echo "Finished at ..."
	date

} 1>>$script.out 2>&1 &
