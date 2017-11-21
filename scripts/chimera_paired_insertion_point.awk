
#	Expecting direction (F or R), pre_or_post (pre or post)

function lg(s){ if( logging ) print s >> log_file }

BEGIN{
	#	input is a sam file so specifying tab may be unnecessary, but clarity is nice.
	FS="\t"
	log_file="chimera_paired_insertion_point.awk"
	if( logging ) print "\nNew file\n\n" >> log_file
}

#  Non-reference lines, with a mapped reference, not matching previous sequence name
( ( !/^@/ ) && ( $3 != "*" ) && ( b[1] != $1 ) ){
	#  Buffer first occurence of sequence name.
	for(i=0;i<=NF;i++)b[i]=$i;

	lg( "Buffer new" )
	lg( $0 )
	next;	#	WILL run next block too, unless next
		#	previously, this resulted in BOTH points included
}  
#  Non-reference lines, with a mapped reference, MATCHING previous sequence name
#
#		Does the pair matching the same reference chromosome matter here?
#		If so, would need to check extended columns for concordant pairing (YT:Z:CP, not YT:Z:DP, YT:Z:UP, ...?)
#		PROPER_PAIR flag is good (same as CP), but would result is very small numbers. Trying.
#
( ( !/^@/ ) && ( $3 != "*" ) && ( b[1] == $1 ) ){
	for(i=0;i<=NF;i++)l[i]=$i;

	lg( "Match" )
	lg( $0 )

#	For the shortest read ...

	lg( "shortest" )
	if( length(l[10]) < length(b[10]) ){
		lg( l[10] )
		for(i=0;i<=NF;i++)s[i]=l[i];
	} else {
		lg( b[10] )
		for(i=0;i<=NF;i++)s[i]=b[i];
	}

	lg( "Human aligned direction: " (and(s[2],16))? "reverse" : "forward" )
	lg( "Viral aligned direction: " direction )
	lg( pre_or_post )
	lg( "insertion point" )

#	s[2] is the flag field
#	 and(s[2],16) = REVERSE
#	!and(s[2],16) = NOT REVERSE = FORWARD

	#	PRE means that it was BEFORE the reference and trimmed on the RIGHT
	#	POST means that it was AFTER the reference and trimmed on the LEFT

#	As i've now reversed the trimmed read, pre and post are irrelevant.  no!

#	PRE - viral align FORWARD - human align FORWARD - insertion point on RIGHT (pos + length)
#	PRE - viral align FORWARD - human align REVERSE - insertion point on LEFT (pos)
#	PRE - viral align REVERSE - TRIMMED/REVERSED - human align FORWARD - insertion point on LEFT (pos)
#	PRE - viral align REVERSE - TRIMMED/REVERSED - human align REVERSE - insertion point on RIGHT (pos + length)
#	POST - viral align FORWARD - human align FORWARD - insertion point on LEFT (pos)
#	POST - viral align FORWARD - human align REVERSE - insertion point on RIGHT (pos + length)
#	POST - viral align REVERSE - TRIMMED/REVERSED - human align FORWARD - insertion point on RIGHT (pos + length)
#	POST - viral align REVERSE - TRIMMED/REVERSED - human align REVERSE - insertion point on LEFT (pos)

	#	direction = VIRAL alignment direction

	if( \
		( !and(s[2],16) && direction == "F" && pre_or_post == "pre" ) ||
		(  and(s[2],16) && direction == "R" && pre_or_post == "pre" ) ||
		(  and(s[2],16) && direction == "F" && pre_or_post == "post" ) ||
		( !and(s[2],16) && direction == "R" && pre_or_post == "post" ) ){
#		( !and(s[2],16) && direction == "F" && pre_or_post == "pre" ) ||
#		(  and(s[2],16) && direction == "R" && pre_or_post == "post" ) ){
			print s[3]"|"s[4]+length(s[10])
			lg( "position plus length" )
			lg( s[4]+length(s[10]) )
	} else if( \
		(  and(s[2],16) && direction == "F" && pre_or_post == "pre" ) ||
		( !and(s[2],16) && direction == "R" && pre_or_post == "pre" ) ||
		( !and(s[2],16) && direction == "F" && pre_or_post == "post" ) ||
		(  and(s[2],16) && direction == "R" && pre_or_post == "post" ) ){
#		(  and(s[2],16) && direction == "R" && pre_or_post == "pre" ) ||
#		( !and(s[2],16) && direction == "F" && pre_or_post == "post" ) ){
			print s[3]"|"s[4]
			lg( "position" )
			lg( s[4] )
	} else {
		lg( "Skipping as does not match alignment orientation" )
	}
	next;
}
