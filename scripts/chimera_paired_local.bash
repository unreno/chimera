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

script=`basename $0`
basedir=`dirname $0`
human='hg19'
viral='herv_k113'
threads=2

function usage(){
	echo
	echo "Usage: (NO EQUALS SIGNS)"
	echo
	echo "Searches for reads that are part viral, part human."
	echo
	echo "$script [--human STRING] [--viral STRING] [--threads INTEGER] -1 <lane 1 fast[aq] file(s)> -2 <lane 2 fast[aq] file(s)"
	echo
	echo "Defaults:"
	echo "  human ..... : $human"
	echo "  viral ..... : $viral"
	echo "  threads ... : $threads (for bowtie2)"
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
		-1)
			shift; lane_1=$1; shift ;;
		-2)
			shift; lane_2=$1; shift ;;
		-*)
			echo ; echo "Unexpected args from: ${*}"; usage ;;
		*)
			break;;
	esac
done


#	Basically, this is TRUE AND DO ...
#	Should be nothing left.
[ $# -ne 0 ] && usage

base=`basename $PWD`

#	Print a lot more stuff
set -x

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
	[ "${lane_1:(-1)}" == "q" -o "${lane_1:(-5)}" == "q.bz2" -o "${lane_1:(-4)}" == "q.gz" ] \
		&& filetype='-q' || filetype='-f'

	base="$base.bowtie2.$viral.very_sensitive_local.paired"
	aligned="$base.aligned"

	#	I think that using --all here would be a good idea, theoretically.
	#	Bowtie2 seems to prefer to soft clip over ends rather than over unknown bp though.
	#	I did try and compare and the final results were identical.
	bowtie2 --very-sensitive-local --threads $threads -x $viral \
		$filetype -1 $lane_1 -2 $lane_2 -S $base.sam

#		$filetype $files | samtools view -b -F 4 > $base.aligned.bam
#	samtools does not seem to process STDIN pipe, so can't do that.
#	Actually, it may but you might have to use - as the filename.

#	I could let the output go to STDOUT then pipe to samtools view -b -F 4 -o $base.bam
#	That would remove the need to convert and delete later.
#	Would save on disk space if that's an issue.
#	May require more memory to do the pipe processing though.

	status=$?
	if [ $status -ne 0 ] ; then
		date
		echo "bowtie failed with $status"
		exit $status
	fi

#	samtools view -b -F 4 -o $aligned.bam $base.sam
#	rm $base.sam


	#	requires bash >= 4.0
	#	${VARIABLE^^} converts to uppercase
	#	${VARIABLE,,} converts to lowercase

	#
	#	Find alignments that align past the appropriate end of the ends of the ltr.
	#
	#    f4 = unmapped
	#    F4 = NOT unmapped = mapped
	#    F8 = mate NOT unmapped = mate mapped
	#
	#	Older versions of awk do not directly support "interval expressions",
	#		ie ({4}, {4,}, {4,6})
	#	Need a newer version or add the --posix option

	samtools view -h $base.sam | awk -v base=$aligned -f $basedir/chimera_paired_trim_aligned_to_fastas.awk


	for pre_or_post in pre post ; do

		#	Align the chimeric reads to the human reference.
		bowtie2 -x $human --threads $threads -f \
			-1 $aligned.${pre_or_post}_1.fasta \
			-2 $aligned.${pre_or_post}_2.fasta \
			-S $aligned.$pre_or_post.bowtie2.$human.sam
		status=$?
		if [ $status -ne 0 ] ; then
			date
			echo "bowtie failed with $status"
			exit $status
		fi

		#	SORT and convert to bam and remove the sam.
		samtools sort -o $aligned.$pre_or_post.bowtie2.$human.bam \
			$aligned.$pre_or_post.bowtie2.$human.sam
		rm $aligned.$pre_or_post.bowtie2.$human.sam

		#	Index now so don't have to before running IGV
		samtools index $aligned.$pre_or_post.bowtie2.$human.bam

	done







#	Which one? How to identify the non-chimeric read without comparing?

#	Full length is a decent idea, but not all sets are the same.

#	Why not just dump both?

#	Only interested in the position of the trimmed end.

#	Assume that the shorter of the pair is the trimmed read.
#	Buffer both reads, output the position of the shortest end?

	echo "Seeking insertion points and overlaps"

	for q in 20 10 00 ; do
		echo $q
		mapq="Q${q}"

#		samtools view -q $q -F 20 $aligned.pre.bowtie2.$human.bam \
#			| awk '( $6 != "*" && $6 != "101M" ){print $3"|"$4+length($10)}' \
#			| sort > $aligned.pre.bowtie2.$human.$mapq.insertion_points
#		samtools view -q $q -F 20 $aligned.post.bowtie2.$human.bam \
#			| awk '( $6 != "*" && $6 != "101M" ){print $3"|"$4}' \
#			| sort > $aligned.post.bowtie2.$human.$mapq.insertion_points
#		chimera_positions_within_10bp.bash $aligned.*.bowtie2.$human.$mapq.insertion_points \
#			| sort | uniq -c > $aligned.both.bowtie2.$human.$mapq.insertion_points.overlappers
#
#		samtools view -q $q -F 4 -f 16 $aligned.pre.bowtie2.$human.bam \
#			| awk '( $6 != "*" && $6 != "101M" ){print $3"|"$4}' \
#			| sort > $aligned.pre.bowtie2.$human.$mapq.rc_insertion_points
#		samtools view -q $q -F 4 -f 16 $aligned.post.bowtie2.$human.bam \
#			| awk '( $6 != "*" && $6 != "101M" ){print $3"|"$4+length($10)}' \
#			| sort > $aligned.post.bowtie2.$human.$mapq.rc_insertion_points
#		chimera_positions_within_10bp.bash $aligned.*.bowtie2.$human.$mapq.rc_insertion_points \
#			| sort | uniq -c > $aligned.both.bowtie2.$human.$mapq.rc_insertion_points.rc_overlappers

		samtools view -q $q -F 20 $aligned.pre.bowtie2.$human.bam \
			| awk -f $basedir/chimera_paired_insertion_point.awk -v direction=F -v pre_or_post=pre \
			| sort > $aligned.pre.bowtie2.$human.$mapq.insertion_points
		samtools view -q $q -F 20 $aligned.post.bowtie2.$human.bam \
			| awk -f $basedir/chimera_paired_insertion_point.awk -v direction=F -v pre_or_post=post \
			| sort > $aligned.post.bowtie2.$human.$mapq.insertion_points
#		chimera_positions_within_10bp.bash $aligned.*.bowtie2.$human.$mapq.insertion_points \
		awk -f $basedir/chimera_positions_within_10bp.awk $aligned.*.bowtie2.$human.$mapq.insertion_points \
			| sort | uniq -c > $aligned.both.bowtie2.$human.$mapq.insertion_points.overlappers

		samtools view -q $q -F 4 -f 16 $aligned.pre.bowtie2.$human.bam \
			| awk -f $basedir/chimera_paired_insertion_point.awk -v direction=R -v pre_or_post=pre \
			| sort > $aligned.pre.bowtie2.$human.$mapq.rc_insertion_points
		samtools view -q $q -F 4 -f 16 $aligned.post.bowtie2.$human.bam \
			| awk -f $basedir/chimera_paired_insertion_point.awk -v direction=R -v pre_or_post=post \
			| sort > $aligned.post.bowtie2.$human.$mapq.rc_insertion_points
#		chimera_positions_within_10bp.bash $aligned.*.bowtie2.$human.$mapq.rc_insertion_points \
		awk -f $basedir/chimera_positions_within_10bp.awk $aligned.*.bowtie2.$human.$mapq.rc_insertion_points \
			| sort | uniq -c > $aligned.both.bowtie2.$human.$mapq.rc_insertion_points.rc_overlappers

	done




	echo
	echo "Finished at ..."
	date

} 1>>$base.$script.out 2>&1 &


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
