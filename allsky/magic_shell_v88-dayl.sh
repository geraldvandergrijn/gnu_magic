# shell script for conversion of cloud index to solar irradiance
 
EXEHOME=/home/akniffka/software/heliosat-hdf/Magic/magic_v88   # path to executable
INPATH=/cmsaf/cmsaf-cld3/RAD_SEVIRI_HelMag/test               # path to the Cloud Albedo files
OPATH=out            # path to the output files
SATELLITE=M7

OLDPWD=`pwd` 
 
# CINMANE  -> file containing the cloud albedo, input.
# GFNMANE  -> file containing the global irradiance SIS, output.
# BFNMANE  -> file containing the direct (beam) irradiance SIS, output.
# ILFNANE  -> file containing the daylight in kLUX, output.
# CFSNAME  -> file containint the cloud forcing solar, output 
# all this file names have to be defined within the for loops!

OUTCNAME="cdl.out"   # only for testing of the output !! regrid flag has to be 1, 
                     # else cdl.out is not produced.

# for script test fixed names are given, names have to 
# be defined in for loops for real runs !!! example GFNAME 
BFNAME=bfile.dat
ILFNAME=il.dat
CFSNAME=cfs.dat


 
for YEAR in 2010
do
    for  MON in 05 #01 02 03 04 05 06 07 08 09 10 11 12   # 01 02 03 04 05   
    do
 	
	# orig OUTPATH=$OPATH$YEAR/$MON
	OUTPATH=$OPATH
	
       # 01 02 03 04 
 	for DAY in 31 #01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31
        do
            
	    doy=$(date -d $MON/$DAY/$YEAR  +"%j")
            echo doy = $doy
            for HH in 13 # 00 01 02 03 04 05 06 07 08 0910 11 12 13 14 15 16 17 18 19 20 21 22 23 
	    do
		for MM in 00  #30
		do
		GFNAME=SISin$YEAR$MON$DAY$HH$MM		
	#	echo "# Year Mon hh mm " $YEAR $MON $HH $MM >> ../temp/rhomax90.out
                    echo "# Year Mon hh mm " $YEAR $MON $HH $MM >> $OPATH/stat.out
	#	    CINAME="STAin"$YEAR$MON$DAY$HH$MM"_${SATELLITE}_VISSN.CI.XPIF"
                    CINAME=$INPATH/"STAin"$YEAR$MON$DAY$HH$MM".h5.CI.XPIF"
        # orig loop
        #            if test -s $INPATH/$YEAR/$MON/$CINAME".gz"
        #            then
        #                   zcat $INPATH/$YEAR/$MON/$CINAME".gz" > $CINAME
        #            fi
        # for script testing
                    if test -s $INPATH/$CINAME".gz"
                    then
                           zcat $INPATH/$CINAME".gz" > $CINAME
                    fi
                    OUTNAME=$YEAR$MON$DAY$HH$MM"G.out"
                    OUTCNAME=$YEAR$MON$DAY$HH$MM".cdl"
		    
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
	     echo input/aeronet_modmed.dat >> magic-config.inp
             echo "##/Dimension/and/name/of/H20/climatology" >> magic-config.inp
             echo 1441 >> magic-config.inp
	     echo 721 >> magic-config.inp
             echo input/h2oclim-era.dat >> magic-config.inp
	     echo "##/name/and/dimension/of/clear/sky/LUT/DOnotCHANGE!" >> magic-config.inp
	     echo input/magic-clear.lut >> magic-config.inp
	     echo input/albedo.map >> magic-config.inp
	     echo input/IGBPa_2006.map >> magic-config.inp
             echo 345 >> magic-config.inp 
             echo "##the|name|of|the|output|files|binary|ASCIIcdl"  >> magic-config.inp
             echo $OUTPATH/$GFNAME >> magic-config.inp
	     echo $OUTPATH/$OUTCNAME >> magic-config.inp
             echo "##/CLOUDinterface/1CloudYeSNo/..." >> magic-config.inp
             echo 1 >> magic-config.inp
             echo 1 256  >> magic-config.inp
	     echo 3712 3712 >> magic-config.inp
	     echo $CINAME  >> magic-config.inp
             echo  0.0  17.8 13.0 >> magic-config.inp
             echo -0.2 1.2 >> magic-config.inp
             echo  0 >> magic-config.inp
             echo "#definition|of|the|region|lon-W|lon-E|lat-N|lat-S|resolution" >> magic-config.inp
             echo 0  120  -40.0 60 0.1 0.1 >> magic-config.inp
             echo $OUTPATH/$BFNAME >> magic-config.inp
             echo $OUTPATH/$ILFNAME >> magic-config.inp
             echo $OUTPATH/$CFSNAME >> magic-config.inp
	     echo "sarbsal" >> magic-config.inp
             # only for testing exit
             #exit
             $EXEHOME/magic-v0x88-dayl.exe
             #rm ./magic-config.inp
     	     #rm $CINAME
	     cd $OUTPATH
	     # ncgen $OUTCNAME -o ${OUTCNAME}.nc  #only needed if regrid flag=1
	     # gzip -f ${OUTCNAME}.nc
             # rm $OUTCNAME
	     #A.K. was soll gezippt werden (die weiteren 4 zeilen)?
             #gzip -f
             #gzip -f 
             #gzip -f
             #gzip -f
	     cd -

		done
		sync
# done
	    done
	done
    done
done
