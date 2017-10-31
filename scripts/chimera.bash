#!/usr/bin/env bash


script=`basename $0`
viral='herv_k113'
human='hg19'
threads=2
distance=15

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

while [ $# -ne 0 ] ; do
	case $1 in
		-h|--h*)
			shift; human=$1; shift ;;
		-v|--v*)
			shift; viral=$1; shift ;;
		-t|--t*)
			shift; threads=$1; shift ;;
		-d|--d*)
			shift; distance=$1; shift ;;
		-*)
			echo ; echo "Unexpected args from: ${*}"; usage ;;
		*)
			break;;
	esac
done


[ $# -eq 0 ] && usage

#/bin/rm -rf working_dir
mkdir -p working_dir
cd working_dir
working_dir=$PWD


#	Print a lot more stuff
set -x

{
	echo "Starting at ..."
	date

	while [ $# -ne 0 ] ; do
		sample=$1

		echo $sample
		cd $working_dir

		#	sample = /genepi2/ccls/data/raw/TCGA_Glioma_HERV52/TCGA-06-0155-10A-01D-0703-09.HERV52.bam.gz
		gz_link=${sample##*/}      # TCGA-06-0155-10A-01D-0703-09.HERV52.bam.gz
		bam=${gz_link%%.gz}        # TCGA-06-0155-10A-01D-0703-09.HERV52.bam
		sample_base=${gz_link%%.*} # TCGA-06-0155-10A-01D-0703-09

		mkdir -p $sample_base
		cd $sample_base

		ln -s $sample
		gunzip -f $gz_link
		#	will also remove link

#		bamToFastq -i $bam -fq $sample_base.1.fastq -fq2 $sample_base.2.fastq
#		samtools fastq $bam -1 $sample_base.1.fastq -2 $sample_base.2.fastq
		samtools fasta $bam -1 $sample_base.1.fasta -2 $sample_base.2.fasta

		rm -f $bam

#		chimera_paired_local.bash --human $human --threads $threads --viral $viral --distance $distance \
#			-1 $sample_base.1.fastq -2 $sample_base.2.fastq
#
#		chimera_unpaired_local.bash --human $human --threads $threads --viral $viral --distance $distance \
#			$sample_base.1.fastq,$sample_base.2.fastq
#
#		rm $sample_base.1.fastq $sample_base.2.fastq

		chimera_paired_local.bash --human $human --threads $threads --viral $viral --distance $distance \
			-1 $sample_base.1.fasta -2 $sample_base.2.fasta

		chimera_unpaired_local.bash --human $human --threads $threads --viral $viral --distance $distance \
			$sample_base.1.fasta,$sample_base.2.fasta

		rm $sample_base.1.fasta $sample_base.2.fasta

		shift
	done

	#	This is kinda important!
	cd $working_dir

	for q in 20 10 00 ; do

		for p in paired unpaired ; do

			insertion_points_to_table.bash \*.${p}\*Q${q}\*points > ${p}_insertion_points_table.Q${q}.csv
			#	= tmpfile. + EXACTLY AS ABOVE + .* (for timestamp)
			mv tmpfile.\*.${p}\*Q${q}\*points.* ${p}_insertion_points.hg19.Q${q}

			chimera_csv_table_group_rows.bash ${p}_insertion_points_table.Q${q}.csv > ${p}_insertion_points_table.Q${q}.grouped.csv

			overlappers_to_table.bash \*.${p}\*Q${q}\*overlappers > ${p}_overlappers_table.Q${q}.csv
			#	= tmpfile. + EXACTLY AS ABOVE + .* (for timestamp)
			mv tmpfile.\*.${p}\*Q${q}\*overlappers.* ${p}_overlappers.hg19.Q${q}

		done	#	paired unpaired

	done	#	20 10 00

	echo
	echo "Finished at ..."
	date

} 1>>$working_dir/$script.out 2>&1 

#	Don't run in background by default
#	&
