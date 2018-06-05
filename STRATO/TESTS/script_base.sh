#!/bin/ksh
 
#set -x 

ycas=$1



for ext in lfi txt fa nc
do
	if [ -f PGD_BASE_$ycas.$ext ] 
	then
		cp -f PGD_BASE_$ycas.$ext PGD_BASE.$ext
	fi
	if [ -f PREP_BASE_$ycas.$ext ] 
	then
		cp -f PREP_BASE_$ycas.$ext PREP_BASE.$ext
	fi
done


