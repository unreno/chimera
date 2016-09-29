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

dir=`dirname $0`

human='hg19'
viral='herv_k113'
threads=2

function usage(){
	echo
	echo "Usage: (NO EQUALS SIGNS)"
	echo
	echo "`basename $0` [--human STRING] [--viral STRING] [--threads INTEGER] <fastq file(s)>"
	echo
	echo "Defaults:"
	echo "  human ..... : $human"
	echo "  viral ..... : $viral"
	echo "  threads ... : $threads (for bowtie2)"
	echo
	echo "Note: all output files will be based on the working directory's name"
	echo "Ignores laning"
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
		-*)
			echo ; echo "Unexpected args from: ${*}"; usage ;;
		*)
			break;;
	esac
done


#       Basically, this is TRUE AND DO ...
[ $# -gt 2 -o $# -eq 0 ] && usage


function positions_within_10bp(){
	for line in `cat $1` ; do
		#	echo "line:" $line
		chr=${line%%|*}	#	remove everything after first pipe (including pipe)
		#	echo "chr: " $chr
		line=${line#*|}	#	remove everything before first pipe (including pipe)
		#	echo "line:" $line
		pos=${line%%|*}	#	remove everything after first pipe (including pipe)
		#	echo "pos: " $pos

		#	Expecting file content format like so ...
		#	chrY|6616930|

		awk -F\| -v chr="$chr" -v pos="$pos" '
			( ( $1 == chr ) && ( (pos-10) < $NF ) && ( (pos+10) > $NF ) ){
				print
			}' $2

	done
}


base=`basename $PWD`

set -x

#	{
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

	#bowtie2
	#-x <bt2-idx> The basename of the index for the reference genome. The basename is the name of any of the index files up to but not including the final .1.bt2 / .rev.1.bt2 / etc. bowtie2 looks for the specified index first in the current directory, then in the directory specified in the BOWTIE2_INDEXES environment variable.


	[[ ${1:(-1)} == 'q' ]] && filetype='-q' || filetype='-f'

	if [ $# -eq 1 ] ; then
		files="-U $1"
	else
		files="-U $1,$2"
	fi

	base="$base.bowtie2.$viral.__very_sensitive_local"

	#	I think that using --all here would be a good idea, theoretically.
	#	Bowtie2 seems to prefer to soft clip over ends rather than over unknown bp though.
	#	I did try and compare and the final results were identical.
	bowtie2 --very-sensitive-local --threads $threads -x $viral \
		$filetype $files -S $base.sam

	status=$?
	if [ $status -ne 0 ] ; then
		date
		echo "bowtie failed with $status"
		exit $status
	fi

	#samtools view -b -S -F 4 -o $base.aligned.bam $base.sam
	samtools view -b -F 4 -o $base.aligned.bam $base.sam
	rm $base.sam

	base="$base.aligned"

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
	#		ie ({4}, {4,}, {4,6])
	#	Need a newer version or add the --posix option

#	samtools view -h -F 4 $base.bam | \
#		gawk -v base=$base -v out="fasta" \
#			-f "$dir/chimera_samtools_extract_and_clip_chimeric_reads.gawk"
	samtools view -h -F 4 $base.bam | gawk -v base=$base \
		'BEGIN {
			pre_out=sprintf("%s.pre.fasta",base)
			post_out=sprintf("%s.post.fasta",base)
		}
		( ( NR % 10000 ) == 0 ){ print "Read "NR" records" }

		( /^@SQ/ ){ ref[substr($2,4)] = substr($3,4) }

		#	Ensure at least 2-digit soft clip and ensure matches near the beginning of the reference.
		( ( $6 ~ /^[0-9]{2,}S[0-9IDM]*$/ ) && ( $4 <= 5 ) ){
			split($6,a,"S")
			clip=a[1]-$4+1
			if( out == "fastq" ){
				print "@"$1"_pre" >> pre_out
				print substr($10,1,clip) >> pre_out
				print "+" >> pre_out
				print substr($11,1,clip) >> pre_out
			} else {
				print ">"$1"_pre" >> pre_out
				print substr($10,1,clip) >> pre_out
			}
		}

		#	Ensure at least 2-digit soft clip and ensure matches near the end of the reference.
		( ( $6 ~ /^[0-9IDM]*[0-9]{2,}S$/ ) && ( $4 >= ( ref[$3] - length($10) + 5 ) ) ){
			clip=ref[$3]-$4+2
			if( out == "fastq" ){
				print "@"$1"_post" >> post_out
				print substr($10,clip) >> post_out
				print "+" >> post_out
				print substr($11,clip) >> post_out
			} else {
				print ">"$1"_post" >> post_out
				print substr($10,clip) >> post_out
			}
		}'
	#	-> pre.fasta
	#	-> post.fasta


	bowtie2 -x $human --threads $threads -f $base.pre.fasta \
		-S $base.pre.bowtie2.$human.sam
	status=$?
	if [ $status -ne 0 ] ; then
		date
		echo "bowtie failed with $status"
		exit $status
	fi

	bowtie2 -x $human --threads $threads -f $base.post.fasta \
		-S $base.post.bowtie2.$human.sam
	status=$?
	if [ $status -ne 0 ] ; then
		date
		echo "bowtie failed with $status"
		exit $status
	fi

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

		samtools view -q $q -F 20 $base.pre.bowtie2.$human.sam \
			| awk '{print $3"|"$4+length($10)}' \
			| sort > $base.pre.bowtie2.$human.$mapq.insertion_points
		samtools view -q $q -F 20 $base.post.bowtie2.$human.sam \
			| awk '{print $3"|"$4}' \
			| sort > $base.post.bowtie2.$human.$mapq.insertion_points
#		$dir/chimera_positions_within_10bp.bash $base.*.bowtie2.$human.$mapq.insertion_points \
		positions_within_10bp $base.*.bowtie2.$human.$mapq.insertion_points \
			| sort | uniq -c > $base.both.bowtie2.$human.$mapq.insertion_points.overlappers

		samtools view -q $q -F 4 -f 16 $base.pre.bowtie2.$human.sam \
			| awk '{print $3"|"$4}' \
			| sort > $base.pre.bowtie2.$human.$mapq.rc_insertion_points
		samtools view -q $q -F 4 -f 16 $base.post.bowtie2.$human.sam \
			| awk '{print $3"|"$4+length($10)}' \
			| sort > $base.post.bowtie2.$human.$mapq.rc_insertion_points
#		$dir/chimera_positions_within_10bp.bash $base.*.bowtie2.$human.$mapq.rc_insertion_points \
		positions_within_10bp $base.*.bowtie2.$human.$mapq.rc_insertion_points \
			| sort | uniq -c > $base.both.bowtie2.$human.$mapq.rc_insertion_points.rc_overlappers

	done

#	samtools view -S -b -o $base.pre.bowtie2.$human.bam  $base.pre.bowtie2.$human.sam
	samtools view -b -o $base.pre.bowtie2.$human.bam  $base.pre.bowtie2.$human.sam
	rm $base.pre.bowtie2.$human.sam

#	samtools view -S -b -o $base.post.bowtie2.$human.bam $base.post.bowtie2.$human.sam
	samtools view -b -o $base.post.bowtie2.$human.bam $base.post.bowtie2.$human.sam
	rm $base.post.bowtie2.$human.sam

	echo
	echo "Finished at ..."
	date

#	} 1>>$base.`basename $0`.out 2>&1
