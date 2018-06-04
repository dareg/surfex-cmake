set -x

tes=$1

for file in script_exec*_0.sh script*_exte_0.sh 
do
	file2=$(basename $file _0.sh)
	if [ "$tes" = "1" ]
	then	
		cat script_debut.sh $file script_fin.sh > ${file2}.sh
	else
		cat $file > ${file2}.sh
	fi
done

