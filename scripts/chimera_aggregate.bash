#!/usr/bin/env bash


script=$( basename $0 )
#viral='herv_k113'
human='hg38'
#threads=4
#distance=15
paired_and_or_unpaired="paired,unpaired"
#bowtie_viral_params="--bowtie_viral_params '--very-sensitive-local'"

function usage(){
	echo
	echo "Usage: (NO EQUALS SIGNS)"
	echo
#	echo "$script [--human STRING] [--viral STRING] [--threads INTEGER] [--distance INTEGER] <SOURCE FILE LIST>"
	echo "$script [--human STRING] [--paired_and_or_unpaired STRING]"
	echo
#	echo "SAMs and BAMs MUST be laned and sorted by name as it is to be split"
#	echo "They can be gzipped or not. Whatever."
#	echo
#	echo "fasta and fastq files, gzipped or not, can also be in the list"
#	echo "They must be paired with lane 1 first and lane 2 second."
	echo
	echo "human : human reference basename (comma separated if multiple)"
	echo
	echo "paired_and_or_unpaired : paired or unpaired (or both comma separated)"
	echo
	echo "Defaults:"
	echo "  human ..... : $human"
	echo "  paired_and_or_unpaired ..... : $paired_and_or_unpaired"
#	echo "  viral ..... : $viral"
#	echo "  threads ... : $threads (for bowtie2)"
#	echo "  distance .. : $distance (for overlappers)"
	echo
	echo
#	echo "Example:"
#	echo "$script ~/genepi2/ccls/data/raw/TCGA_Glioma_HERV52/TCGA-??-????-10*bam.gz"
	echo
	exit
}

while [ $# -ne 0 ] ; do
	case $1 in
		-h|--h*)
			shift; human=$1; shift ;;
#		-v|--v*)
#			shift; viral=$1; shift ;;
#		-t|--t*)
#			shift; threads=$1; shift ;;
#		-d|--d*)
#			shift; distance=$1; shift ;;
#		--paired_and_or_unpaired)
		-p|--p*)
			shift; paired_and_or_unpaired=$1; shift ;;
#		--bowtie_viral_params)
##			shift; bowtie_viral_params="--bowtie_viral_params \"${1}\""; shift ;;
##	simply can't pass them together.
#			shift; bvp1="--bowtie_viral_params"; bvp2="$1"; shift ;;
		-*)
			echo ; echo "Unexpected args from: ${*}"; usage ;;
		*)
			break;;
	esac
done


[ $# -ne 0 ] && usage

#base_dir=$PWD

working_dir=$PWD

#	Print a lot more stuff
set -x

{
	echo "Starting at ..."
	date

#	#	This is kinda important!
#	cd $working_dir

	for q in 20 10 00 ; do

		for p in ${paired_and_or_unpaired/,/ } ; do

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



				#	chimera_csv_table_group_rows.bash input NEEDS to be sorted
				head -1 ${p}_insertion_points_table.${h}.Q${q}.csv > ${p}_insertion_points_table.${h}.Q${q}.sorted.csv
				tail -n +2 ${p}_insertion_points_table.${h}.Q${q}.csv \
					| sort -t \| -k 1,1 -k 2n,2 >> ${p}_insertion_points_table.${h}.Q${q}.sorted.csv

				chimera_csv_table_group_rows.bash ${p}_insertion_points_table.${h}.Q${q}.sorted.csv \
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
