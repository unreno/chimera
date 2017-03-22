
#	Expecting direction (F or R), pre_or_post (pre or post)

BEGIN{
	#	input is a sam file so specifying tab may be unnecessary, but clarity is nice.
	FS="\t"
}

#  Non-reference lines, with a mapped reference, not matching previous sequence name
( ( !/^@/ ) && ( $3 != "*" ) && ( b[1] != $1 ) ){
	#  Buffer first occurence of sequence name.
	for(i=0;i<=NF;i++)b[i]=$i;
}  
#  Non-reference lines, with a mapped reference, MATCHING previous sequence name
( ( !/^@/ ) && ( $3 != "*" ) && ( b[1] == $1 ) ){
	for(i=0;i<=NF;i++)l[i]=$i;

#	For the shortest read ...

	if( length(l[10]) < length(b[10]) ){
		for(i=0;i<=NF;i++)s[i]=l[i];
	} else {
		for(i=0;i<=NF;i++)s[i]=b[i];
	}

	if( 
		( direction == "F" && pre_or_post == "pre" ) ||
		( direction == "R" && pre_or_post == "post" ) ){
			print s[3]"|"s[4]+length(s[10])
	} else if(
		( direction == "R" && pre_or_post == "pre" ) ||
		( direction == "F" && pre_or_post == "post" ) ){
			print s[3]"|"s[4]
	}

}
