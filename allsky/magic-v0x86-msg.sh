# shell script for conversion of cloud index to solar irradiance
 
EXEHOME=/home/rmueller/archive/software/magicsol-v01/allsky
INPATH=satimages
OPATH=out

 
for YEAR in 2007
do
 OUTPATH=$OPATH
 #MON=11
    for  MON in 03 #01 02 03 04 05 06 07 08 09 10 11 12   # 01 02 03 04 05   
    do
       # 01 02 03 04 
 	for DAY in 07 10 # 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31
        do
            
	    doy=$(date -d $MON/$DAY/$YEAR  +"%j")
            echo doy = $doy
            for HH in 12 # 05 06 07 08 09 10 11 12 13 14 15 16 17 18 
	    do
		for MM in 00 # 00 15 30 45 
		do
	#	echo "# Year Mon hh mm " $YEAR $MON $HH $MM >> ../temp/rhomax90.out
                    echo "# Year Mon hh mm " $YEAR $MON $HH $MM >> $OPATH/stat.out
		   
# CINAME has to be changed for MVIRI, Meteosat1-7: 
		    CINAME=$YEAR$MON$DAY$HH$MM"_MSG_VIS008.CI.XPIF"
                    OUTNAME=$YEAR$MON$DAY$HH$MM"G.out"
		    
		    echo   "write config file"
                    echo "##/year/month/day/hour/minute(GMT)" > magic-config.inp #> magic-config.inp
		    echo $YEAR >> magic-config.inp 
		    echo $MON >> magic-config.inp 
		    echo $doy >> magic-config.inp
		    echo $HH >> magic-config.inp
		    echo $MM >> magic-config.inp
                    echo "##/Dimension/and/name/of/aerosol/climatology" >> magic-config.inp
                    echo 360 >> magic-config.inp
		    echo 180 >> magic-config.inp
	     echo climatologies/aeronet_modmed.dat >> magic-config.inp
             echo "##/Dimension/and/name/of/H20/climatology" >> magic-config.inp
             echo 1441 >> magic-config.inp
	     echo 721 >> magic-config.inp
             echo climatologies/h2oclim-era.dat >> magic-config.inp
	     echo "##/name/and/dimension/of/clear/sky/LUT/DOnotCHANGE!" >> magic-config.inp
	     echo luts/magic-clear.lut >> magic-config.inp
	     echo climatologies/albedo.map >> magic-config.inp
	     echo climatologies/IGBPa_2006.map >> magic-config.inp
             echo 345 >> magic-config.inp 
             echo "##the|name|of|the|output|files|binary|ASCIIcdl"  >> magic-config.inp
             echo $OPATH"/"$OUTNAME >> magic-config.inp
	     echo $OPATH"/"$OUTNAME".cdl" >> magic-config.inp
             echo "##/CLOUDinterface/1CloudYeSNo/..." >> magic-config.inp
             echo 1 >> magic-config.inp
             echo 2 256  >> magic-config.inp
# next line 5000 5000 for MVIRI !
	     echo 3712 3712 >> magic-config.inp
	     echo satimages/$CINAME  >> magic-config.inp
# next line 0.0 18.0 12.0 for MVIRI
             echo  0.0  17.83 12.0 >> magic-config.inp
             echo -0.2 1.2 >> magic-config.inp
             echo  1 >> magic-config.inp
             echo "#definition|of|the|region|lon-W|lon-E|lat-N|lat-S|resolution" >> magic-config.inp
             echo -10  25  33.0 55 0.1 0.1 >> magic-config.inp
             echo geo/nlon.dat >> magic-config.inp
             echo geo/nlat.dat >> magic-config.inp
             echo geo/nsza.dat >> magic-config.inp
             $EXEHOME/magic-v0x86.exe
             rm $EXEHOME/magic-config.inp
             ncgen $OPATH"/"$OUTNAME".cdl" -o $OPATH"/"$OUTNAME".nc"
             rm $OPATH"/"$OUTNAME".cdl"
		done
		sync
# done
	    done
	done
    done
done
