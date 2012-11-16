 # shell sript for DWD Offenbach
 # eff. cloud albedo calculation for VIS008 XPIF images 
 # Dr. Richard Mueller
 # Free University of Palatnum
 # August 2012
 #######################################################
 
HELHOME=/home/lee/Documents/eclipse-workspace/gnu_magic/helclim-vG/exe

 GZPFAD=/home/lee/Downloads/goes/062003_regridded_decompressed # here are the compressed SAT images
 #GZPFAD=$HELHOME/in # here are the compressed Meteosat images
 PFAD=/home/lee/Downloads/goes/magicsol_tmp # the uncompressed images will be mv to this directory
 #CIPFAD=$HELHOME/out/CI_VIS008  # here are the cloud index output files
CIPFAD=/home/lee/Downloads/goes/magicsol_output
GPFAD=$CIPFAD
# here are the surface albedo images
 
sat=goes09 

for YEAR in 2003 #2005 2006 2007 2008
do
 #MON=11
    for  MON in 06 #01 02 03 04 05 06 07 08 09 10 11 12   # 01 02 03 04 05   
    do
	for HH in 01 # 13 05 06 07 08 09 10 11 12 14 15 16 17  #10 16 #11 12 13 14 15 #16 #13 #12 #09 12 15  
# 05 06 07 08 09 10 11 12 13 14 15 16 17 
	do
	    for MM in 25  #00 15 30 45 #30
	    do
	#	echo "# Year Mon hh mm " $YEAR $MON $HH $MM >> ../temp/rhomax90.out
                echo "# Year Mon hh mm " $YEAR $MON $HH $MM >> ../temp/rhostat.out
#15 30 45
		for DAY in {152..181} #153 154 155 156 157 158 159 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 181
#01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31
#31
		do
		    NAME=$GZPFAD/$sat.$YEAR.$DAY.$HH$MM".out"
		    GZNAME=$GZPFAD/$NAME.gz
		    BZNAME=$GZPFAD/$NAME.bz2
    # echo $BZNAME
    # die Bilder auspacken, falls vorhanden
                    if test -s $NAME
                    then
                    echo  cp $NAME nach $PFAD/$NAME
                    cp $NAME $PFAD
                    fi
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
		cd $PFAD ;  ls $sat.$YEAR.*.$HH$MM".out" >> $PFAD/liste 
                echo $sat.$YEAR.*.$HH$MM".out"

# orig -s 23
#$HELHOME/Cloudindex/heliosat0x8 -c $CIPFAD -g $GPFAD -l $PFAD/liste -b $PFAD -s 23 -z 0 
# nach Eumetsat konf wieder herstellen !!
#orig 3.9.2008		$HELHOME/Cloudindex/tmp300708.exe -c $CIPFAD -g $GPFAD -l $PFAD/liste -b $PFAD -s -3 -z 0
# $HELHOME/Cloudindex/tmp300708.exe -c $CIPFAD -g $GPFAD -l $PFAD/liste -b $PFAD -s -3 -z 0
# $HELHOME/Cloudindex/heliosat-32b-v1x23.exe -c $CIPFAD -g $GPFAD -l $PFAD/liste -b $PFAD -s -3 -z 0
# gms position 155 E
$HELHOME/magicsol-v0goes.exe -c $CIPFAD -g $GPFAD -l $PFAD/liste -b $PFAD -s 155 -z 99

rm  $PFAD/liste
# rm $PFAD/*"_MSG_VIS008.XPIF"
	    done
	    sync
# done
	done
       rm $PFAD/rhomax.out
    done
done



