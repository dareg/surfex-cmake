#set -x

#script_pgd_physio.sh*
#script_pgd_ecoclimap.sh*
#script_prep.sh*
#script_offline.sh*
#script_cas_parts.sh*
for script in script_pgd_physio.sh script_pgd_ecoclimap.sh script_prep.sh script_offline.sh script_cas_parts.sh
do
	echo $script
	if [ -f _${script} ]; then

		cp -f $script file1.sh
		for cas in `cat _${script}`
		do
			./comment.sh $cas
		done
		mv -f file1.sh .${script}
	fi
done

#script_water.sh* script_water_canopy.sh*
#script_flake.sh* script_flake_canopy.sh*
#script_csts.sh* script_csts_canopy.sh*
#script_teb.sh* script_teb_prep_simple.sh*


#script_sea.sh* script_sea_canopy.sh* script_1d_ocean.sh*
#script_isba.sh* script_isba_photo.sh* script_isba_canopy_parall.sh*
#script_isba.sh* script_isba_photo_dif.sh* script_isba_canopy_parall.sh*
#script_isba.sh* script_isba_phys.sh* script_isba_canopy.sh*
#script_isba.sh* script_isba_phys.sh* script_isba_canopy_parall.sh*


#script_isba.sh* script_isba_photo.sh* script_isba_neige.sh* script_isba_canopy_parall.sh*
#script_isba.sh* script_isba_photo_dif.sh* script_isba_neige_dif.sh* script_isba_canopy_parall.sh*
#script_teb.sh* script_teb_pgd.sh* script_teb_prep.sh* script_teb_run.sh*





