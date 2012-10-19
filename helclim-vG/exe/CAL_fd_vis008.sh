 # shell sript for DWD Offenbach
 # eff. cloud albedo calculation for VIS008 XPIF images 
 # Dr. Richard Mueller
 # Free University of Palatnum
 # August 2012
 #######################################################
 
HELHOME=/home/rmueller/archive/software/magicsol-v01/helclim/exe

 GZPFAD=/cmsaf/cmsaf-rad2/data/MSGraw # here are the compressed Meteosat images
 #GZPFAD=$HELHOME/in # here are the compressed Meteosat images
 PFAD=/cmsaf/cmsaf-rad2/temp1 # the uncompressed images will be mv to this directory
 #CIPFAD=$HELHOME/out/CI_VIS008  # here are the cloud index output files
CIPFAD=/cmsaf/cmsaf-rad2/data/heliosatOUT/VIS008/2007
GPFAD=/cmsaf/cmsaf-rad2/data/heliosatOUT/VIS008/2007 
# here are the surface albedo images
 
 
for YEAR in 2007 #2005 2006 2007 2008
do
 #MON=11
    for  MON in 01 02 03 04 05 06 07 08 09 10 11 12   # 01 02 03 04 05   
    do
	for HH in 13 05 06 07 08 09 10 11 12 14 15 16 17  #10 16 #11 12 13 14 15 #16 #13 #12 #09 12 15  
# 05 06 07 08 09 10 11 12 13 14 15 16 17 
	do
	    for MM in 00  #00 15 30 45 #30
	    do
	#	echo "# Year Mon hh mm " $YEAR $MON $HH $MM >> ../temp/rhomax90.out
                echo "# Year Mon hh mm " $YEAR $MON $HH $MM >> ../temp/rhostat.out
#15 30 45
		for DAY in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31
#31
		do
		    NAME=$YEAR$MON$DAY$HH$MM"_MSG_VIS008.XPIF"
		    GZNAME=$GZPFAD/$YEAR$MON/$NAME.gz
		    BZNAME=$GZPFAD/$YEAR$MON/$NAME.bz2
    # echo $BZNAME
    # die Bilder auspacken, falls vorhanden
		    if test -s $GZNAME
		    then
			echo auspacken $GZNAME nach $PFAD/$NAME
			gunzip -c $GZNAME > $PFAD/$NAME
		    fi
		    if test -s $BZNAME
		    then
			echo auspacken $BZNAME nach $PFAD/$NAME
			bunzip2 -c $BZNAME > $PFAD/$NAME
		    fi
		done
 # Liste erzeugen. diese enthaelt die Bildnamen
		cd $PFAD ;  ls $YEAR$MON??$HH$MM*.XPIF >> $PFAD/liste 
# orig -s 23
#$HELHOME/Cloudindex/heliosat0x8 -c $CIPFAD -g $GPFAD -l $PFAD/liste -b $PFAD -s 23 -z 0 
# nach Eumetsat konf wieder herstellen !!
#orig 3.9.2008		$HELHOME/Cloudindex/tmp300708.exe -c $CIPFAD -g $GPFAD -l $PFAD/liste -b $PFAD -s -3 -z 0
# $HELHOME/Cloudindex/tmp300708.exe -c $CIPFAD -g $GPFAD -l $PFAD/liste -b $PFAD -s -3 -z 0
# $HELHOME/Cloudindex/heliosat-32b-v1x23.exe -c $CIPFAD -g $GPFAD -l $PFAD/liste -b $PFAD -s -3 -z 0
$HELHOME/magicsol-v0.exe -c $CIPFAD -g $GPFAD -l $PFAD/liste -b $PFAD -s -3 -z 0

 rm  $PFAD/liste
 rm $PFAD/*"_MSG_VIS008.XPIF"
	    done
	    sync
# done
	done
       rm $PFAD/rhomax.out
    done
done



