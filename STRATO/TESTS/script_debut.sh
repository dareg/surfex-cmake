script=$(basename $2 .txt)
script=`echo $script | cut -c6- | tr [:upper:] [:lower:]`

for file in _script_$script.sh
do
	for cas in `cat $file`
	do
		nam=`echo $cas | cut -c3-`
		if [ "$nam" = "$1" ]
		then
			run=`echo $cas | cut -c1`
			if [ "$run" = "1" ]
			then

