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
threads=2
distance=10

function usage(){
	echo
	echo "Usage: (NO EQUALS SIGNS)"
	echo
	echo "Searches for reads that are part viral, part human."
	echo
	echo "$script [--human STRING] [--viral STRING] [--threads INTEGER] [--distance INTEGER] [-1 <lane 1 fast[aq] file(s)> -2 <lane 2 fast[aq] file(s)>] [--bam <sorted by name, laned bam file to be split>]"
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
	echo " but MUST be separated by commas without white space."
	echo "They will be passed directly to bowtie2 as such."
	echo "For example:"
	echo "  -1 a_1.fastq.gz,b_1.fastq.gz -2 a_2.fastq.gz,b_2.fastq.gz"
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
		-1)
			shift; lane_1=$1; shift ;;
		-2)
			shift; lane_2=$1; shift ;;
		-b|--b*)
			shift; bam=$1; shift ;;
		-*)
			echo ; echo "Unexpected args from: ${*}"; usage ;;
		*)
			break;;
	esac
done


#	Basically, this is TRUE AND DO ...
#	Should be nothing left.
[ $# -ne 0 ] && usage
#[ -z $lane_1 ] && usage
#[ -z $lane_2 ] && usage

if [ -n "$bam" ] ; then
	if [ -n "$lane_1" -o -n "$lane_2" ] ; then
		echo "Bam AND fasta/q files provided?"
		echo "Bam:$bam:1:$lane_1:2:$lane_2:"
		usage
	else
		#	lane_1=${bam%.*}.1.fastq
		#	lane_2=${bam%.*}.2.fastq
		lane_1=${bam%.*}.1.fasta
		lane_2=${bam%.*}.2.fasta
	fi
else
	if [ -z "$lane_1" ] || [ -z "$lane_2" ] ; then
		usage
	fi
fi

#echo "bam:$bam:"
#echo "lane_1:$lane_1:"
#echo "lane_2:$lane_2:"

base=$( basename $PWD )

#	Print a lot more stuff
set -x

samtools --version
bowtie2 --version

{
	echo "Starting at ..."
	date

	if [ -n "$bam" -a -r "$bam" ] ; then
		samtools fasta -N -1 $lane_1 -2 $lane_2 $bam
	fi


	#indexes=/Volumes/cube/working/indexes
	#	leading with the ": " stops execution
	#	just ${BOWTIE2_INDEXES:"/Volumes/cube/working/indexes"}
	#	would try to execute the result.  I just want the OR/EQUALS feature
	#: ${BOWTIE2_INDEXES:="/Volumes/cube/working/indexes"}
	: ${BOWTIE2_INDEXES:="$HOME/BOWTIE2_INDEXES"}

	#	they MUST be exported, to be picked up by bowtie2
	export BOWTIE2_INDEXES

	#	One line "if-then-else" to determine filetype by last character of first file name.
	[ "${lane_1:(-1)}" == "q" -o "${lane_1:(-5)}" == "q.bz2" -o "${lane_1:(-4)}" == "q.gz" ] \
		&& filetype='-q' || filetype='-f'

	base="$base.bowtie2.$viral.very_sensitive_local.paired"
	aligned="$base.aligned"

	#	gawk script filters those where at least 1 read aligned
	bowtie2 --very-sensitive-local --threads $threads -x $viral \
		$filetype -1 $lane_1 -2 $lane_2 \
		| gawk -F"\t" '
			( /^@/ ){ print; next; }
			( xor ( !and($2,4), !and($2,8) ) ){ print }
		' \
		| awk -v base=$aligned -f $basedir/chimera_paired_trim_aligned_to_fastas.awk



	for human_ref in ${human/,/ } ; do

		#	PRE means that it was BEFORE the reference and trimmed on the RIGHT (unless reversed)
		#	POST means that it was AFTER the reference and trimmed on the LEFT (unless reversed)
	
		for viral_direction in F R ; do
	
			for pre_or_post in pre post ; do
	
				if [ -f $aligned.$viral_direction.${pre_or_post}_1.fasta -a \
						 -f $aligned.$viral_direction.${pre_or_post}_2.fasta ] ; then
	
					#	Align the chimeric reads to the human reference.
					bowtie2 -x $human_ref --threads $threads -f \
						-1 $aligned.$viral_direction.${pre_or_post}_1.fasta \
						-2 $aligned.$viral_direction.${pre_or_post}_2.fasta \
						-S $aligned.$viral_direction.$pre_or_post.bowtie2.$human_ref.sam
					status=$?
					if [ $status -ne 0 ] ; then
						date
						echo "bowtie failed with $status"
						exit $status
					fi
	
					#	SORT BY NAME and convert to bam
					samtools sort -n -o $aligned.$viral_direction.$pre_or_post.bowtie2.$human_ref.name.bam \
						$aligned.$viral_direction.$pre_or_post.bowtie2.$human_ref.sam
	
					#	SORT and convert to bam and remove the sam.
					samtools sort -o $aligned.$viral_direction.$pre_or_post.bowtie2.$human_ref.position.bam \
						$aligned.$viral_direction.$pre_or_post.bowtie2.$human_ref.sam
					rm $aligned.$viral_direction.$pre_or_post.bowtie2.$human_ref.sam
	
					#	Index now so don't have to before running IGV
					samtools index $aligned.$viral_direction.$pre_or_post.bowtie2.$human_ref.position.bam
	
				else
					echo "$aligned.$viral_direction.${pre_or_post}_1.fasta or"
					echo "$aligned.$viral_direction.${pre_or_post}_2.fasta not found"
				fi
	
			done
	
		done
	
	
		echo "Seeking insertion points and overlaps"
	
		for q in 20 10 00 ; do
			echo $q
			mapq="Q${q}"
	
			#	use name sorted bam as awk script expecting them in order
			#	add -f 2 to get PROPER_PAIR
	
			#	-f = ALL of ...
			#	-F = NONE of ...
	
			#  2 = read and mate aligned a proper concordant pair (very desireable)
			#  4 = read UNMAPPED ( -F 4 = MAPPED )
			#  8 = read's MATE UNMAPPED
			# 16 = read REVERSED
			# 18 = 16,2 ( -f 18 = proper pair alignment, reversed )
			# 20 = 16,4 ( -F 20 = not unmapped, not reverse = MAPPED FORWARD )
			# 32 = read's MATE REVERSED
	
			#	This is getting funky.
			#	Can the reads be aligned proper paired, but in different directions? YES YES YES, sadly.
			#		unlikely to have REVERSE and MREVERSE
	
			#	chimera_paired_insertion_point.awk expecting BOTH READS! Remove the 4 flag. Ehh...
	
			#	As the flags are read specific, using them here and expecting samtools to
			#	return both mates is seeming to be impossible.
			#	The logic will need to be moved to the awk script.
			#	based on shortest (trimmed) read and given direction and pre_or_post
	
			# -f 2 -F 12 (-F 12 is redundant as -f 2 implies both aligned)
	
			#	Add -v logging=1 to chimera_paired_insertion_point.awk call when debugging
			#	Also chimera_paired_trim_aligned_to_fastas.awk above
	
			samtools view -q $q -f 2 $aligned.F.pre.bowtie2.$human_ref.name.bam \
				| awk -f $basedir/chimera_paired_insertion_point.awk -v direction=F -v pre_or_post=pre \
				| sort > $aligned.pre.bowtie2.$human_ref.$mapq.insertion_points
	
			samtools view -q $q -f 2 $aligned.F.post.bowtie2.$human_ref.name.bam \
				| awk -f $basedir/chimera_paired_insertion_point.awk -v direction=F -v pre_or_post=post \
				| sort > $aligned.post.bowtie2.$human_ref.$mapq.insertion_points
	
	
			#	This will likely find nothing.
			awk -v distance=$distance -f $basedir/chimera_positions_within.awk \
				$aligned.*.bowtie2.$human_ref.$mapq.insertion_points \
				| sort | uniq -c > $aligned.both.bowtie2.$human_ref.$mapq.insertion_points.overlappers
	
	
			samtools view -q $q -f 2 $aligned.R.pre.bowtie2.$human_ref.name.bam \
				| awk -f $basedir/chimera_paired_insertion_point.awk -v direction=R -v pre_or_post=pre \
				| sort > $aligned.pre.bowtie2.$human_ref.$mapq.rc_insertion_points
	
			samtools view -q $q -f 2 $aligned.R.post.bowtie2.$human_ref.name.bam \
				| awk -f $basedir/chimera_paired_insertion_point.awk -v direction=R -v pre_or_post=post \
				| sort > $aligned.post.bowtie2.$human_ref.$mapq.rc_insertion_points
	
	
			#	This will likely find nothing.
			awk -v distance=$distance -f $basedir/chimera_positions_within.awk \
				$aligned.*.bowtie2.$human_ref.$mapq.rc_insertion_points \
				| sort | uniq -c > $aligned.both.bowtie2.$human_ref.$mapq.rc_insertion_points.rc_overlappers
	
	
			#	I think that this one would imply that the virus is inserted reverse complement
			awk -v distance=$distance -f $basedir/chimera_positions_within.awk \
				$aligned.post.bowtie2.$human_ref.$mapq.rc_insertion_points \
				$aligned.pre.bowtie2.$human_ref.$mapq.insertion_points \
				| sort | uniq -c > $aligned.both.bowtie2.$human_ref.$mapq.frc_insertion_points.frc_overlappers
	
	
			#	And this one would imply that the virus is inserted forward.
			#	Lane 1 of the PRE is forward, Lane 2 is reverse complement
			#	This one is rather popular
			awk -v distance=$distance -f $basedir/chimera_positions_within.awk \
				$aligned.post.bowtie2.$human_ref.$mapq.insertion_points \
				$aligned.pre.bowtie2.$human_ref.$mapq.rc_insertion_points \
				| sort | uniq -c > $aligned.both.bowtie2.$human_ref.$mapq.rcf_insertion_points.rcf_overlappers
	
		done

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
