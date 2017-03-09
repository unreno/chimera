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


function positions_within_10bp(){
	#
	#	When called like ...
	#	positions_within_10bp $base.*.bowtie2.$human.$mapq.insertion_points ....
	#		$1 will be $base.post.bowtie2.$human.$mapq.insertion_points
	#		$2 will be $base.pre.bowtie2.$human.$mapq.insertion_points
	#
	#	Both can be of the format ...
	#		chr8|99063786
	#		chr9|105365407
	#		chrX|74554211
	#		chrY|14564844
	#		... OR ...
	#		chrX:154433612-155433612|68228
	#		chrX:154433612-155433612|790825
	#		chrX:93085499-94085499|110644
	#		chrX:93085499-94085499|112146
	#		... OR EVEN ...
	#		chr8|99063786|F
	#		chr9|105365407|F
	#		chrX|74554211|R
	#		chrY|14564844|R
	#
	for line in `cat $1` ; do
		#	echo "line:" $line
		chr=${line%%|*}	#	remove everything after first pipe (including pipe)
		#	echo "chr: " $chr    #	chrX or chrX:154433612-155433612
		line=${line#*|}	#	remove everything before first pipe (including pipe)
		#	echo "line:" $line   #	74554211 or 790825
		pos=${line%%|*}	#	remove everything after first pipe (including pipe)
		#	echo "pos: " $pos    #	74554211 or 790825

		#	Print out the lines with the same reference chromosome
		#		and a position within 10bp in either direction.
		awk -F\| -v chr="$chr" -v pos="$pos" '
			( ( $1 == chr ) && ( (pos-10) < $2 ) && ( (pos+10) > $2 ) ){
				print
			}' $2

	done
}


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
#	[ "${1:(-1)}" == 'q' ] && filetype='-q' || filetype='-f'
	[ "${lane_1:(-1)}" == "q" -o "${lane_1:(-5)}" == "q.bz2" -o "${lane_1:(-4)}" == "q.gz" ] \
		&& filetype='-q' || filetype='-f'

	base="$base.bowtie2.$viral.__very_sensitive_local"

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

	aligned="$base.aligned"
	samtools view -b -F 4 -o $aligned.bam $base.sam
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






#			PAIRED NOW.  NEED TO UPDATE FROM HERE DOWN.

#		GIVEN THE EXTENSIVE LOGIC NEEDED, LIKELY NEED TO CHECK FULL FILE NOT JUST THOSE ALIGNED.

#		As likely won't be getting insertion points from this, 
#		we won't need to differentiate pre or post?

#	Eventually, the database will contain multiple elements.
#	HERVs, SINEs, LINEs, SVAs, ...
#	Include this in the output read? Separate files?



	samtools view -h $base.sam | gawk -v base=$aligned \
		'function print_to_fasta(a){
			lane=(and(a[2],64))?"1":"2";
			print ">"a[1]"/"lane >> base"_"lane".fasta"
			print a[10]          >> base"_"lane".fasta"
		}
		function trim(r){
			#	Ensure at least 2-digit soft clip and ensure matches near the beginning of the reference.
			#	"near the beginning" means starts at position <= 5
			if( ( r[6] ~ /^[0-9]{2,}S[0-9IDM]*$/ ) && ( r[4] <= 5 ) ){
print "Trimming beginning"
				split(r[6],a,"S")
				clip=a[1]-r[4]+1
#				print ">"$1"_pre" >> pre_out
#				print substr($10,1,clip) >> pre_out
				r[10]=substr(r[10],1,clip)
			}

			#	Ensure at least 2-digit soft clip and ensure matches near the end of the reference.
			#	"near the end" means starts at position >= 5 more than the reference minus the length of the read
			if( ( r[6] ~ /^[0-9IDM]*[0-9]{2,}S$/ ) && ( r[4] >= ( ref[r[3]] - length(r[10]) + 5 ) ) ){
print "trimming end"
				clip=ref[r[3]]-r[4]+2
#				print ">"$1"_post" >> post_out
#				print substr($10,clip) >> post_out
				r[10]=substr(r[10],clip)
			}
#			return r;	#	cannot return arrays. Arrays passed as reference so mods made are actual.
		}
		BEGIN {
#			out[1]=sprintf("%s.1.fasta",base)
#			out[2]=sprintf("%s.2.fasta",base)
			split("",b);split("",l);
		}
		#	Store all of the reference lengths
		#	Ex. @SQ	SN:chr1	LN:249250621
		#	... ref[chr1] = 249250621
		( /^@SQ/ ){ ref[substr($2,4)] = substr($3,4); next; }

		#	Simply for progress
		( ( !/@SQ/ ) && ( ( NR % 10000 ) == 0 ) ){ print "Read "NR" records" }

		#	Non-sequence reference lines, with a mapped reference, matching previous sequence name
		( ( !/@SQ/ ) && ( $3 != "*" ) && ( b[1] == $1 ) ){
			for(i=0;i<=NF;i++)l[i]=$i;

			#	1 and only 1 read aligned. Many ways to check this.
			if( ( and(l[2],4) && !and(l[2],8) ) || ( !and(l[2],4) && and(l[2],8) ) ){
				#	If this read unmapped and mate not unmapped (mate mapped) ...
				if( and(l[2],4) && !and(l[2],8) ){
					trim(b);
				#	If mate read unmapped and this not unmapped (this mapped) ...
				}else if( !and(l[2],4) && and(l[2],8) ){
					trim(l);
				}
				print_to_fasta( b )
				print_to_fasta( l )
			}

			delete b; delete l;
			next; #	do not buffer this line
		}

		#	Non-sequence reference lines, with a mapped reference, not matching previous sequence name
		#	Buffer first occurence of sequence name.
		( ( !/@SQ/ ) && ( $3 != "*" ) && ( b[1] != $1 ) ){
			for(i=0;i<=NF;i++)b[i]=$i;
		}'
	#	-> .1.fasta
	#	-> .2.fasta


#	for pre_or_post in pre post ; do
#
#		#	Align the chimeric reads to the human reference.
#		bowtie2 -x $human --threads $threads -f -U $aligned.$pre_or_post.fasta \
#			-S $aligned.$pre_or_post.bowtie2.$human.sam
#		status=$?
#		if [ $status -ne 0 ] ; then
#			date
#			echo "bowtie failed with $status"
#			exit $status
#		fi
#
#		#	Convert to bam and remove the sam.
#		samtools view -b -o $aligned.$pre_or_post.bowtie2.$human.bam \
#			$aligned.$pre_or_post.bowtie2.$human.sam
#		rm $aligned.$pre_or_post.bowtie2.$human.sam
#
#
##	Sadly, this won't work as is.
##	The TCGA data uses sequence names like "@SRR448740.1 1/1"
##	When read as unlaned by bowtie, the name truncates at the space.
##	By doing this the lane is unattainable. We would need to do some renaming.
##	That would take a bit of time and disk space. What a waste.
#
##	Or perhaps run each lane separately into clearly marked files?
#
##		gawk -v base=$base -v pre_or_post=$pre_or_post '
##		BEGIN {
##			out[1]=sprintf("%s.%s.1.fasta",base,pre_or_post)
##			out[2]=sprintf("%s.%s.2.fasta",base,pre_or_post)
##		}
##		( ( NR == FNR ) && ( /^>/ ) ){
##			sub(/^>/,"");
##			split($0,name,"/");
##			names[name[1]]=name[2];
##			seqs[name[1]]="";
##			last=name[1];
##		}
##		( ( NR == FNR ) && ( !/^>/ ) ){
##			seqs[last]=seqs[last] $0
##		}
##		( ( NR != FNR ) && ( ! /^@/ ) ) {
##			split($1,name,"/");
##			#	sam file: sequence name in fasta, not the same lane, and unmapped
##			if( ( name[1] in names ) && ( name[2] != names[name[1]] ) && ( and($2,4) ) ) ){
##
##				print name[1] "/" names[name[1]] >> out[names[name[1]]]
##				print seqs[name[1]] >> out[names[name[1]]]
##				print name[1] "/" name[2] >> out[name[2]]
##				print $10 >> out[name[2]]
##
##			}
##		}' $aligned.$pre_or_post.fasta $base.sam
#
#	done

	#	find insertion points
	#	then find those with the signature overlap

	#	 f = ALL/YES
	#	 F = NONE/NOT	(results in double negatives)
	#	 4 = not aligned
	#	 8 = mate not aligned
	#	16 = reverse complement

#	echo "Seeking insertion points and overlaps"
#
#	for q in 20 10 00 ; do
#		echo $q
#		mapq="Q${q}"
#
#		samtools view -q $q -F 20 $aligned.pre.bowtie2.$human.bam \
#			| awk '{print $3"|"$4+length($10)}' \
#			| sort > $aligned.pre.bowtie2.$human.$mapq.insertion_points
#		samtools view -q $q -F 20 $aligned.post.bowtie2.$human.bam \
#			| awk '{print $3"|"$4}' \
#			| sort > $aligned.post.bowtie2.$human.$mapq.insertion_points
#		positions_within_10bp $aligned.*.bowtie2.$human.$mapq.insertion_points \
#			| sort | uniq -c > $aligned.both.bowtie2.$human.$mapq.insertion_points.overlappers
#
#		samtools view -q $q -F 4 -f 16 $aligned.pre.bowtie2.$human.bam \
#			| awk '{print $3"|"$4}' \
#			| sort > $aligned.pre.bowtie2.$human.$mapq.rc_insertion_points
#		samtools view -q $q -F 4 -f 16 $aligned.post.bowtie2.$human.bam \
#			| awk '{print $3"|"$4+length($10)}' \
#			| sort > $aligned.post.bowtie2.$human.$mapq.rc_insertion_points
#		positions_within_10bp $aligned.*.bowtie2.$human.$mapq.rc_insertion_points \
#			| sort | uniq -c > $aligned.both.bowtie2.$human.$mapq.rc_insertion_points.rc_overlappers
#
#	done

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
