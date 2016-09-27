#
#	Older versions of awk do not directly support "interval expressions", 
#		ie ({4}, {4,}, {4,6])
#	Need a newer version or add the --posix option

#	--posix NEEDS to be AFTER any -v settings!

#	$cmd | awk -v base=$base -v out=$out --posix '


#	UNTESTED SINCE EXTRACTION


BEGIN {
	pre_out=sprintf("%s.pre.%s",base,out)
	post_out=sprintf("%s.post.%s",base,out)
}
( ( NR % 10000 ) == 0 ){ print "Read "NR" records" }

( /^@SQ/ ){ ref[substr($2,4)] = substr($3,4) }

#( ( $6 ~ /^[0-9]{2,}S[0-9IDM]*$/ ) && ( $4 <= 5 ) && ( $3 ~ /beginning/ ) ){
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

#( ( $6 ~ /^[0-9IDM]*[0-9]{2,}S$/ ) && ( $4 >= ( ref[$3] - (length($10)*0.9) ) ) && ( $3 ~ /ending/ ) ){
#( ( $6 ~ /^[0-9IDM]*[0-9]{2,}S$/ ) && ( $4 >= ( ref[$3] - (length($10)*0.9) ) ) ){
( ( $6 ~ /^[0-9IDM]*[0-9]{2,}S$/ ) && ( $4 >= ( ref[$3] - length($10) - 5 ) ) ){
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
