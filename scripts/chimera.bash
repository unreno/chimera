#!/usr/bin/env bash


script=$( basename $0 )
viral='herv_k113'
human='hg19'
threads=2
distance=15
paired_and_or_unpaired="paired unpaired"

function usage(){
	echo
	echo "Usage: (NO EQUALS SIGNS)"
	echo
	echo "$script [--human STRING] [--viral STRING] [--threads INTEGER] [--distance INTEGER] <SOURCE FILE LIST>"
	echo
	echo "SAMs and BAMs MUST be laned and sorted by name as it is to be split"
	echo "They can be gzipped or not. Whatever."
	echo
	echo "fasta and fastq files, gzipped or not, can also be in the list"
	echo "They must be paired with lane 1 first and lane 2 second."
	echo
	echo "Defaults:"
	echo "  human ..... : $human"
	echo "  viral ..... : $viral"
	echo "  threads ... : $threads (for bowtie2)"
	echo "  distance .. : $distance (for overlappers)"
	echo
	echo
	echo "Example:"
	echo "$script ~/genepi2/ccls/data/raw/TCGA_Glioma_HERV52/TCGA-??-????-10*bam.gz"
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
		--paired_and_or_unpaired)
			shift; paired_and_or_unpaired=$1; shift ;;
		-*)
			echo ; echo "Unexpected args from: ${*}"; usage ;;
		*)
			break;;
	esac
done


[ $# -eq 0 ] && usage

base_dir=$PWD

#/bin/rm -rf working_dir
mkdir -p working_dir
cd working_dir
working_dir=$PWD


#	Print a lot more stuff
set -x

samtools --version
bowtie2 --version

{
	echo "Starting at ..."
	date

	while [ $# -ne 0 ] ; do
		#	on second loop, would be in working_dir which would generate bad path for that sample
		cd $base_dir

		#	realpath returns a full absolute path without links
		#sample=$( realpath $1 )
		#	Keep the link names rather than the targets
		sample=$( realpath --no-symlinks $1 )
		#	could also use "readlink -f $1"

		echo $sample
		if [ ${sample:(-7)} == ".bam.gz" \
				-o ${sample:(-7)} == ".sam.gz" \
				-o ${sample:(-4)} == ".bam" \
				-o ${sample:(-4)} == ".sam" \
		] ; then
			echo "Filetype: sam or bam"
			filetype="_am"
		elif [ ${sample:(-9)} == ".fasta.gz" \
				-o ${sample:(-9)} == ".fastq.gz" \
				-o ${sample:(-6)} == ".fasta" \
				-o ${sample:(-6)} == ".fastq" \
				-o ${sample:(-6)} == ".fa.gz" \
				-o ${sample:(-6)} == ".fq.gz" \
				-o ${sample:(-3)} == ".fa" \
				-o ${sample:(-3)} == ".fq" \
		] ; then
			echo "Filetype: fasta or fastq"
			filetype="fast_"
			fast_1_with_path=$( realpath $1 )
			shift
			fast_2_with_path=$( realpath $1 )
		else
			echo "Unknown filetype so skipping"
			break
		fi

		cd $working_dir

		sample_dirname=$( basename $sample )	#	just the file name, no path

		#	strip off any expected extensions
		sample_dirname=${sample_dirname%.gz}
		#	really can only be one of the following
		sample_dirname=${sample_dirname%.bam}
		sample_dirname=${sample_dirname%.sam}
		sample_dirname=${sample_dirname%.fasta}
		sample_dirname=${sample_dirname%.fastq}
		sample_dirname=${sample_dirname%.fa}
		sample_dirname=${sample_dirname%.fq}
		mkdir -p $sample_dirname
		cd $sample_dirname

		sample_basename=$( basename $sample )	#	just the file name, no path
		#	file list could be sample1.bam and sample1.bam.gz
		#	so remember any and all of the extensions
#	this may be moot now
		initial_sample_basename=$sample_basename

		#	for .sam.gz and .bam.gz
		if [ ${sample_basename:(-5)} == "am.gz" ] ; then
			ln -s $sample

			#	need the -f to gunzip a soft linked file.
			#	this will also remove the link
			gunzip -f $sample_basename

			sample_basename=${sample_basename%.gz}
		fi

		#	for .sam and .bam
		if [ ${filetype} == "_am" ] ; then
			echo "sam or bam"

			#	If not initially gzipped, need to link the sam/bam file.
			if [ ${sample:(-2)} == "am" ] ; then
				ln -s $sample
			fi

			fast_1=$initial_sample_basename.1.fasta
			fast_2=$initial_sample_basename.2.fasta
			samtools fasta $sample_basename -1 $fast_1 -2 $fast_2
			rm -f $sample_basename
		elif [ ${filetype} == "fast_" ] ; then
			echo "fasta or fastq"
			fast_1=$( basename $fast_1_with_path )
			fast_2=$( basename $fast_2_with_path )
			ln -s $fast_1_with_path
			ln -s $fast_2_with_path
		fi

		for p in $paired_and_or_unpaired ; do

			if [ $p == "paired" ] ; then
				chimera_paired_local.bash --human $human --threads $threads \
					--viral $viral --distance $distance \
					-1 $fast_1 -2 $fast_2
			fi

			if [ $p == "unpaired" ] ; then
				chimera_unpaired_local.bash --human $human --threads $threads \
					--viral $viral --distance $distance ${fast_1},${fast_2}
			fi

		done

		rm -f $fast_1 $fast_2

		#	write protect everything
		chmod -R -w .

		shift
	done

	#	This is kinda important!
	cd $working_dir

	for q in 20 10 00 ; do

		for p in $paired_and_or_unpaired ; do

			for h in ${human/,/ } ; do

				#	Given that some of our runs were very large, an actual file listing
				#	was too much for the shell to handle. Had to use a glob pattern that
				#	will be used by find. This has gotten rather awkward now.

				#	And includes files that shouldn't be included

				chimera_insertion_points_to_table.bash \*.${p}\*.${h}.\*Q${q}\*points \
					> ${p}_insertion_points_table.${h}.Q${q}.csv
				#	This script generates a tmpfile that is just a list of ALL the insertion points.
				#	It is unused, but kept for curiousity.
				#	= tmpfile. + EXACTLY AS ABOVE + .* (for timestamp)
				mv tmpfile.\*.${p}\*.${h}.\*Q${q}\*points.* ${p}_insertion_points.${h}.Q${q}

				chimera_csv_table_group_rows.bash ${p}_insertion_points_table.Q${q}.csv \
					> ${p}_insertion_points_table.${h}.Q${q}.grouped.csv

				#	this is a TINY bit different as it preserves full file names.
				chimera_overlappers_to_table.bash \*.${p}\*.${h}.\*Q${q}\*overlappers \
					> ${p}_overlappers_table.${h}.Q${q}.csv
				#	This script generates a tmpfile that is just a list of ALL the insertion points.
				#	It is unused, but kept for curiousity.
				#	= tmpfile. + EXACTLY AS ABOVE + .* (for timestamp)
				mv tmpfile.\*.${p}\*.${h}.\*Q${q}\*overlappers.* ${p}_overlappers.${h}.Q${q}

			done	#	human refs, possibly separated by comma "hg19,hg38"

		done	#	paired unpaired

	done	#	20 10 00

	echo
	echo "Finished at ..."
	date

} 1>>$working_dir/$script.out 2>&1 

#	Don't run in background by default
#	&
