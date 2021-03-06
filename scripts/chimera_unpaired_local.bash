#!/usr/bin/env bash


#
#	Can't pass options when using env
#
#	with the -x the commands are sent to STDERR before execution
#
#	bowtie output goes to stderr for some reason
#	probably because the SAM file usually goes to stdout
#	so, wrap everything in curly braces and direct both
#	to files.
#
#	Explicit redirection within the block will override this
#


#	This really needs to be run in the data directory.
#
#	Eventually, may want to pass number of cpus or threads so
#	execs can use the same number.

script=$( basename $0 )
basedir=$( dirname $0 )
human='hg38'
viral='herv_k113'
threads=4
distance=10
bowtie_viral_params="--very-sensitive-local"
#bowtie_viral_params="--local --mp 4 --ma 4 -D 15 -R 5 -i S,1,0.25 -L 3 -N 1"

function usage(){
	echo
	echo "Usage: (NO EQUALS SIGNS)"
	echo
	echo "Searches for reads that are part viral, part human."
	echo "Completely ignores any laning."
	echo
	echo "$script [--human STRING] [--viral STRING] [--threads INTEGER] [--distance INTEGER] <fastq file(s)>"
	echo
	echo "Defaults:"
	echo "  human ..... : $human"
	echo "  viral ..... : $viral"
	echo "  threads ... : $threads (for bowtie2)"
	echo "  distance .. : $distance (for overlappers)"
	echo
	echo "Note: all output files will be based on the working directory's name"
	echo
	echo "Note: fastq files can be gzipped (.gz) or bzipped (.bz2),"
#	Not true, I now join them. This doesn't work in other scripts where they are laned.
#	I would have to add a special parameter to determine which is which.
#	echo " but MUST be separated by commas without white space."
#	echo "They will be passed directly to bowtie2 as such."
#	echo "For example:"
#	echo "  a.fastq.gz,b.fastq.gz"
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
		--bowtie_viral_params)
			shift; bowtie_viral_params=${1}; shift ;;
#			shift; bowtie_viral_params="--bowtie_viral_params '${1}'"; shift ;;
		-*)
			echo ; echo "Unexpected args from: ${*}"; usage ;;
		*)
			break;;
	esac
done


#       Basically, this is TRUE AND DO ...
[ $# -eq 0 ] && usage

base=$( basename $PWD )

#	Print a lot more stuff
set -x

samtools --version
bowtie2 --version

{
	echo "Starting at ..."
	date

	#indexes=/Volumes/cube/working/indexes
	#	leading with the ": " stops execution
	#	just ${BOWTIE2_INDEXES:"/Volumes/cube/working/indexes"}
	#	would try to execute the result.  I just want the OR/EQUALS feature
	#: ${BOWTIE2_INDEXES:="/Volumes/cube/working/indexes"}
	: ${BOWTIE2_INDEXES:="$HOME/BOWTIE2_INDEXES"}

	#	they MUST be exported, to be picked up by bowtie2
	export BOWTIE2_INDEXES

	#	One line "if-then-else" to determine filetype by last character of first file name.
#	[ "${1:(-1)}" == 'q' ] && filetype='-q' || filetype='-f'
	[ "${1:(-1)}" == "q" -o "${1:(-5)}" == "q.bz2" -o "${1:(-4)}" == "q.gz" ] \
		&& filetype='-q' || filetype='-f'

	#	Create a string of all files to be passed to bowtie2
	files=$(echo $* | awk '{printf "-U ";for(i=1;i<NF;i++){printf "%s,",$i};print $NF}')

	base="$base.bowtie2.$( basename ${viral} ).very_sensitive_local.unpaired"


	#	I think that I can do this in one step.
	aligned="$base.aligned"
	#	--very-sensitive-local \
	#	--local --mp 4 --ma 3 -D 25 -R 4 -i S,1,0.50 -L 3 -N 1 \
	#	--local --mp 4 --ma 4 -D 15 -R 5 -i S,1,0.25 -L 3 -N 1 \
echo ":${bowtie_viral_params}:"
	bowtie2 --threads $threads -x $viral \
		${bowtie_viral_params} \
		$filetype $files \
		| samtools view -h -F 4 - \
		| gawk -v base=$aligned -f $basedir/chimera_unpaired_trim_aligned_to_fastas.awk


	for human_ref in ${human/,/ } ; do
		human_ref_base=$( basename ${human_ref} )

		for pre_or_post in pre post ; do
	
			if [ -f $aligned.$pre_or_post.fasta ] ; then

				#	Align the chimeric reads to the human reference.
				bowtie2 -x $human_ref --threads $threads -f -U $aligned.$pre_or_post.fasta \
					-S $aligned.$pre_or_post.bowtie2.${human_ref_base}.sam
				status=$?
				if [ $status -ne 0 ] ; then
					date
					echo "bowtie failed with $status"
					exit $status
				fi
	
				#	SORT and convert to bam and remove the sam.
				samtools sort -o $aligned.$pre_or_post.bowtie2.${human_ref_base}.position.bam \
					$aligned.$pre_or_post.bowtie2.${human_ref_base}.sam
				rm $aligned.$pre_or_post.bowtie2.${human_ref_base}.sam
	
				#	Index now so don't have to before running IGV
				samtools index $aligned.$pre_or_post.bowtie2.${human_ref_base}.position.bam

			else
				echo "No $aligned.$pre_or_post.fasta found"
			fi
	
		done	#	for pre_or_post in pre post ; do
	
		if [ -f $aligned.pre.bowtie2.${human_ref_base}.position.bam -a \
				 -f $aligned.post.bowtie2.${human_ref_base}.position.bam ] ; then

			samtools merge $aligned.bowtie2.${human_ref_base}.position.bam \
				$aligned.pre.bowtie2.${human_ref_base}.position.bam \
				$aligned.post.bowtie2.${human_ref_base}.position.bam

		elif [ -f $aligned.pre.bowtie2.${human_ref_base}.position.bam ] ; then

			cp $aligned.pre.bowtie2.${human_ref_base}.position.bam $aligned.bowtie2.${human_ref_base}.position.bam

		elif [ -f $aligned.post.bowtie2.${human_ref_base}.position.bam ] ; then

			cp $aligned.post.bowtie2.${human_ref_base}.position.bam $aligned.bowtie2.${human_ref_base}.position.bam

		else

			echo "Likely gonna crash. No pre or post found."

		fi

		samtools index $aligned.bowtie2.${human_ref_base}.position.bam

		#	find insertion points
		#	then find those with the signature overlap
	
		#	 f = ALL/YES
		#	 F = NONE/NOT	(results in double negatives)
		#	 4 = not aligned
		#	 8 = mate not aligned
		#	16 = reverse complement
	
		echo "Seeking insertion points and overlaps"
	
		for q in 20 10 00 ; do
			echo $q
			mapq="Q${q}"


			if [ -f $aligned.pre.bowtie2.${human_ref_base}.position.bam ] ; then
	
				samtools view -q $q -F 20 $aligned.pre.bowtie2.${human_ref_base}.position.bam \
					| awk '{print $3"|"$4+length($10)}' \
					| sort > $aligned.pre.bowtie2.${human_ref_base}.$mapq.insertion_points

				samtools view -q $q -F 4 -f 16 $aligned.pre.bowtie2.${human_ref_base}.position.bam \
					| awk '{print $3"|"$4}' \
					| sort > $aligned.pre.bowtie2.${human_ref_base}.$mapq.rc_insertion_points

			fi

			if [ -f $aligned.post.bowtie2.${human_ref_base}.position.bam ] ; then

				samtools view -q $q -F 20 $aligned.post.bowtie2.${human_ref_base}.position.bam \
					| awk '{print $3"|"$4}' \
					| sort > $aligned.post.bowtie2.${human_ref_base}.$mapq.insertion_points

				samtools view -q $q -F 4 -f 16 $aligned.post.bowtie2.${human_ref_base}.position.bam \
					| awk '{print $3"|"$4+length($10)}' \
					| sort > $aligned.post.bowtie2.${human_ref_base}.$mapq.rc_insertion_points

			fi


			if [ -f $aligned.pre.bowtie2.${human_ref_base}.$mapq.insertion_points -a \
					 -f $aligned.post.bowtie2.${human_ref_base}.$mapq.insertion_points ] ; then

				awk -v distance=$distance -f $basedir/chimera_positions_within.awk \
					$aligned.*.bowtie2.${human_ref_base}.$mapq.insertion_points \
					| sort | uniq -c > $aligned.both.bowtie2.${human_ref_base}.$mapq.insertion_points.overlappers

			fi

			if [ -f $aligned.pre.bowtie2.${human_ref_base}.$mapq.rc_insertion_points -a \
					 -f $aligned.post.bowtie2.${human_ref_base}.$mapq.rc_insertion_points ] ; then

				awk -v distance=$distance -f $basedir/chimera_positions_within.awk \
					$aligned.*.bowtie2.${human_ref_base}.$mapq.rc_insertion_points \
					| sort | uniq -c > $aligned.both.bowtie2.${human_ref_base}.$mapq.rc_insertion_points.rc_overlappers

			fi
	
		done	#	for q in 20 10 00 ; do

	done	#	for human_ref in ${human/,/ } ; do

	echo
	echo "Finished at ..."
	date

} 1>>$base.$script.out 2>&1 

#	Don't run in the background by default
#	&


#		REFERENCES
#
#	Flags:
#		1    0x1   PAIRED        .. paired-end (or multiple-segment) sequencing technology
#		2    0x2   PROPER_PAIR   .. each segment properly aligned according to the aligner
#		4    0x4   UNMAP         .. segment unmapped
#		8    0x8   MUNMAP        .. next segment in the template unmapped
#		16   0x10  REVERSE       .. SEQ is reverse complemented
#		32   0x20  MREVERSE      .. SEQ of the next segment in the template is reversed
#		64   0x40  READ1         .. the first segment in the template
#		128  0x80  READ2         .. the last segment in the template
#		256  0x100 SECONDARY     .. secondary alignment
#		512  0x200 QCFAIL        .. not passing quality controls
#		1024 0x400 DUP           .. PCR or optical duplicate
#		2048 0x800 SUPPLEMENTARY .. supplementary alignment
#
#
#	Sam file columns
#	1 QNAME String Query template NAME
#	2 FLAG Int bitwise FLAG
#	3 RNAME String Reference sequence NAME
#	4 POS Int 1-based leftmost mapping POSition
#	5 MAPQ Int MAPping Quality
#	6 CIGAR String CIGAR string
#	7 RNEXT String Ref.  name of the mate/next read
#	8 PNEXT Int Position of the mate/next read
#	9 TLEN Int observed Template LENgth
#	10 SEQ String segment SEQuence
#	11 QUAL String
#
