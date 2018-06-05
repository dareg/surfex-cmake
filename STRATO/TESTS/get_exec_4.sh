script1=_$1
script2=_$2
script3=_$3
script4=_$4

#set -x

cp -f $script1 file2.sh
#cas de n1
for cas1 in `cat $script1`
do

  run1=`echo $cas1 | cut -c1`
  nam1=`echo $cas1 | cut -c3-`
  #cas d'un cas qui en appelle d'autres n2
  if [ "$run1" != "0" -a "$run1" != "1" ]
  then
    #on supprime la ligne dans le fichier racine
    sed -e '/'$cas1'$/d' file2.sh > file1.sh

    #pour tous les cas appelés n2
    for cas2 in `cat $script2`
    do

      run2=`echo $cas2 | cut -c1`
      nam2=`echo $cas2 | cut -c3-`
      #si ce cas n2 en appelle d'autres n3
      if [ "$run2" = "." ]; then

	#pour tous les cas appelés n3
	for cas3 in `cat $script3`
	do

          run3=`echo $cas3 | cut -c1`
	  nam3=`echo $cas3 | cut -c3-`
	  #si ce cas n3 en appelle d'autres n4
	  if [ "$run3" = "." ]; then

            #pour tous les cas appelés n4
            for cas4 in `cat $script4`
            do

              run4=`echo $cas4 | cut -c1`
	      nam4=`echo $cas4 | cut -c3-`
	      #si ce cas n'en appelle pas d'autres
	      if [ "$run4" = "0" -o "$run3" = "1" ]; then
		    
		#on écrit le cas concaténé dans le fichier racine
	        echo "0-"$nam1$nam2$nam3$nam4 >> file1.sh 

              fi
            done

	  #si ce cas n3 n'en appelle pas d'autres
          else

            echo "0-"$nam1$nam2$nam3 >>  file1.sh 
#
          fi	  
        done

      #si ce cas n2 n'en appelle pas d'autres
      else

	echo "0-"$nam1$nam2 >> file1.sh

      fi
    done
    mv -f file1.sh file2.sh
  fi
done

mv -f file2.sh $script1

