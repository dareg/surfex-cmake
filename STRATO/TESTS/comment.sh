run=`echo $1 | cut -c1`
name=`echo $1 | cut -c3-`

if [ $run -eq 1 ]; then

	awk -F " " '{\
		if ('"/\"$name\"/"') { \
			if ( $1=="#") {$1=""; print $0} \
	       		else {print $0}} \
		else {print $0}}	' file1.sh > file2.sh

		mv -f file2.sh file1.sh

elif [ $run -eq 0 ]; then

	awk -F " " '{\
		if ('"/\"$name\"/"') { \
			if ( $1!="#") { print "# "$0} \
	       		else {print $0}} \
		else {print $0}}	' file1.sh > file2.sh

		mv -f file2.sh file1.sh

fi
