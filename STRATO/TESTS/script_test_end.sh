#!/bin/ksh
#set -x

ycas=$1
ymod=$2
yexp=$3

yname=$(echo $ymod | cut -c 6-)
yend="END_$yname"
yplante="PLANTE_$yname"

if [ ! -f $yend ]
then
	touch $yend
fi
if [ ! -f $yplante ]
then
	touch $yplante
fi

if [ -f LISTING_SODA0.txt ]
then
	ends=$(grep ENDS\ CORRECTLY LISTING_SODA0.txt)
	if [ ! "$ends" = "" ] 
	then
		echo "$yexp: $ycas" >> $yend
	else
		echo "$yexp: $ycas" >> $yplante
	fi
else
	if [ -f LISTING_OFFLINE0.txt ]
	#if [ -f LISTING_PREP0.txt ]
	then
		ends=$(grep ENDS\ CORRECTLY LISTING_OFFLINE0.txt)
		#ends=$(grep ENDS\ CORRECTLY LISTING_PREP0.txt)
		if [ ! "$ends" = "" ] 
		then
			echo "$yexp: $ycas" >> $yend
		else
			echo "$yexp: $ycas" >> $yplante
		fi
	else
		echo "$yexp: $ycas" >> $yplante
	
	#elif [ -f LISTING_PREP.txt ]
	#then
#		#ends=$(grep ENDS\ CORRECTLY LISTING_OFFLINE0.txt)
#		ends=$(grep ENDS\ CORRECTLY LISTING_PREP.txt)
#		if [ ! "$ends" = "" ] 
#		then
#			echo "$yexp: $ycas" >> $yend
#		else
#			echo "$yexp: $ycas" >> $yplante
#		fi		
#	else
#		echo "$yexp: $ycas" >> $yplante
	fi
fi

rm -f LISTING*

