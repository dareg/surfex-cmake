file=$1
script=$2

awk '{if ('"/$script/"') {print $0}}' $file >> _$file	
sed -e 's/\"\ $fname\ $2\ $3//g' _$file > __$file
sed -e 's/.\/'$script'\ \"OPTIONS.nam\"\ \"/\.-/g' __$file > _$file

rm -f __$file
