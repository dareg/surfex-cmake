file=$1
script=$2

awk '{if ('"/$script/"') {print $0}}' $file >> _$file	
sed -e 's/\"\ $3\ $4\ $5//g' _$file > __$file
sed -e 's/.\/'$script'\ \"OPTIONS.nam\"\ \"/\.-/g' __$file > _$file
sed -e 's/$2//g' _$file > __$file
mv -f __$file _$file

rm -f __$file
