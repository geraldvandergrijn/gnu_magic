# shell script for conversion of cloud index to solar irradiance
 
# Wo die executables leben
EXEHOME=/home/lee/Downloads/Sitzung_Richard/MAGOCSOL-v0/allsky
# Path to cloud albedo
INPATH=/home/lee/Downloads/Sitzung_Richard/Beispieldaten/cloud_index_output/200506
# Output path
OPATH=/home/lee/Downloads/Sitzung_Richard/Beispieldaten/ergebnisse
# Name of satellite
SATELLITE=M7

# ?
OLDPWD=`pwd`
 

## Achtung: Schleife immer anpassen!
for YEAR in 2005
do
 #MON=11
    for  MON in 06 #01 02 03 04 05 06 07 08 09 10 11 12   # 01 02 03 04 05   
    do
 	OUTPATH=$OPATH/$YEAR$MON
	mkdir $OUTPATH
       # 01 02 03 04 
 	for DAY in  15 #01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31
        do
            
	    doy=$(date -d $MON/$DAY/$YEAR  +"%j")
            echo doy = $doy
            for HH in 13 # 00 01 02 03 04 05 06 07 08 10 11 12 13 14 15 16 17 18 19 20 21 22 23 
	    do
		for MM in 00 # 30
		do
	#	echo "# Year Mon hh mm " $YEAR $MON $HH $MM >> ../temp/rhomax90.out
                    echo "# Year Mon hh mm " $YEAR $MON $HH $MM >> $OPATH/stat.out
		    
		    CINAME=$YEAR$MON$DAY$HH$MM"_${SATELLITE}_VISSN.CI.XPIF"
                     if test -s $INPATH/$YEAR/$MON/$CINAME".gz"
                     then
 			    echo 'I unpack.'
                            zcat $INPATH/$YEAR/$MON/$CINAME".gz" > $CINAME
                     fi

                    OUTNAME=$YEAR$MON$DAY$HH$MM"G.out"
                    OUTCNAME=$YEAR$MON$DAY$HH$MM
		    
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
             echo $OUTPATH/output-bin >> magic-config.inp
	     echo $OUTPATH/$OUTCNAME >> magic-config.inp
             echo "##/CLOUDinterface/1CloudYeSNo/..." >> magic-config.inp
             echo 1 >> magic-config.inp
             echo 1 256  >> magic-config.inp
	     echo 5000 5000 >> magic-config.inp
	     echo $INPATH"/"$CINAME  >> magic-config.inp
	    ## Here aufpassen - stimmen die Koordinaten???
             echo  0.0  18.0 25.0 >> magic-config.inp
             echo -0.2 1.2 >> magic-config.inp
             echo  1 >> magic-config.inp
             echo "#definition|of|the|region|lon-W|lon-E|lat-N|lat-S|resolution" >> magic-config.inp
	    ## Hier aufpassen - Dimensionen des Ausschnitts angeben mit Auflösung (NN resampling auf ein reguläres Gitter)
             echo -90  90  -90.0 90 0.05 0.05 >> magic-config.inp
             echo geo/nlon.dat >> magic-config.inp
             echo geo/nlat.dat >> magic-config.inp
             echo geo/nsza.dat >> magic-config.inp
             $EXEHOME/magic-v0x86.exe
             #rm ./magic-config.inp
	     rm $CINAME
	     cd $OUTPATH
	    ## Hier macht das Programm ein netCDF aus den Daten.
	     ncgen $OUTCNAME -o ${OUTCNAME}.nc
	     gzip -f ${OUTCNAME}.nc
	     rm $OUTCNAME
	     cd -

		done
		sync
# done
	    done
	done
    done
done
