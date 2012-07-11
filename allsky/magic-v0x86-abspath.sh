# shell script for conversion of cloud index to solar irradiance
 
EXEHOME=/home/rmueller/archive/software/hel4eum/n2g
INPATH=/cmsaf/cmsaf-rad1/data/heliosatOUT/cloud
OPATH=/cmsaf/cmsaf-rad1/data/heliosatOUT/G4PI/

 
for YEAR in 2006
do
 OUTPATH=$OPATH$YEAR
 #MON=11
    for  MON in 01 08 #01 02 03 04 05 06 07 08 09 10 11 12   # 01 02 03 04 05   
    do
       # 01 02 03 04 
 	for DAY in 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31
        do
            
	    doy=$(date -d $MON/$DAY/$YEAR  +"%j")
            echo doy = $doy
            for HH in 05 06 07 08 09 10 11 12 13 14 15 16 17  #10 16 #11 12 13 14 15 #16 #13 #12 #09 12 15  
# 05 06 07 08 09 10 11 12 13 14 15 16 17 
	    do
		for MM in 00 15 30 45 #30
		do
	#	echo "# Year Mon hh mm " $YEAR $MON $HH $MM >> ../temp/rhomax90.out
                    echo "# Year Mon hh mm " $YEAR $MON $HH $MM >> $OPATH/stat.out
		    
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
	     echo /home/rmueller/archive/software/hel4eum/n2g/input/aeronet_modmed.dat >> magic-config.inp
             echo "##/Dimension/and/name/of/H20/climatology" >> magic-config.inp
             echo 1441 >> magic-config.inp
	     echo 721 >> magic-config.inp
             echo /home/rmueller/archive/software/hel4eum/n2g/input/h2oclim-era.dat >> magic-config.inp
	     echo "##/name/and/dimension/of/clear/sky/LUT/DOnotCHANGE!" >> magic-config.inp
	     echo /home/rmueller/archive/software/hel4eum/n2g/input/magic-clear.lut >> magic-config.inp
	     echo /home/rmueller/archive/software/hel4eum/n2g/input/albedo.map >> magic-config.inp
	     echo /home/rmueller/archive/software/hel4eum/n2g/input/IGBPa_2006.map >> magic-config.inp
             echo 345 >> magic-config.inp 
             echo "##the|name|of|the|output|files|binary|ASCIIcdl"  >> magic-config.inp
             echo /cmsaf/cmsaf-rad1/data/heliosatOUT/G4PI/$YEAR$MON/$OUTNAME >> magic-config.inp
	     echo /cmsaf/cmsaf-rad1/data/heliosatOUT/G4PI/output-cdl >> magic-config.inp
             echo "##/CLOUDinterface/1CloudYeSNo/..." >> magic-config.inp
             echo 1 >> magic-config.inp
             echo 2 256  >> magic-config.inp
	     echo 3712 3712 >> magic-config.inp
	     echo /cmsaf/cmsaf-rad1/data/heliosatOUT/cloud/$CINAME  >> magic-config.inp
             echo  0.0  17.83 12.0 >> magic-config.inp
             echo -0.2 1.2 >> magic-config.inp
             echo  1 >> magic-config.inp
             echo "#definition|of|the|region|lon-W|lon-E|lat-N|lat-S|resolution" >> magic-config.inp
             echo -60  60  -40.0 45 0.15 0.15 >> magic-config.inp
             echo geo/nlon.dat >> magic-config.inp
             echo geo/nlat.dat >> magic-config.inp
             echo geo/nsza.dat >> magic-config.inp
             $EXEHOME/v0x85.exe
             rm $EXEHOME/magic-config.inp
		done
		sync
# done
	    done
	done
    done
done
