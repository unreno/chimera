
#	Expecting direction (F or R), pre_or_post (pre or post)

BEGIN{
	#	input is a sam file so specifying tab may be unnecessary, but clarity is nice.
	FS="\t"
	log="chimera_paired_insertion_point.awk."
}

#  Non-reference lines, with a mapped reference, not matching previous sequence name
( ( !/^@/ ) && ( $3 != "*" ) && ( b[1] != $1 ) ){
	#  Buffer first occurence of sequence name.
	for(i=0;i<=NF;i++)b[i]=$i;

	print "Buffer new" >> log
	print $0 >> log
}  
#  Non-reference lines, with a mapped reference, MATCHING previous sequence name
( ( !/^@/ ) && ( $3 != "*" ) && ( b[1] == $1 ) ){
	for(i=0;i<=NF;i++)l[i]=$i;

	print "Match" >> log
	print $0 >> log

#	For the shortest read ...

	print "shortest" >> log
	if( length(l[10]) < length(b[10]) ){
		print l[10] >> log
		for(i=0;i<=NF;i++)s[i]=l[i];
	} else {
		print b[10] >> log
		for(i=0;i<=NF;i++)s[i]=b[i];
	}

	print direction >> log
	print pre_or_post >> log
	print "insertion point" >> log

	if( \
		( direction == "F" && pre_or_post == "pre" ) ||
		( direction == "R" && pre_or_post == "post" ) ){
			print s[3]"|"s[4]+length(s[10])
			print "position plus length"
			print s[4]+length(s[10]) >> log
	} else if( \
		( direction == "R" && pre_or_post == "pre" ) ||
		( direction == "F" && pre_or_post == "post" ) ){
			print s[3]"|"s[4]
			print "position"
			print s[4] >> log
	}

}
