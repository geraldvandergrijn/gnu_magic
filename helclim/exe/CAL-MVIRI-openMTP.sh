 # shell sript for DWD Offenbach
 # cloud index calculation 
 # Dr. Annette Hammer 
 # University of Oldenburg
 # 11. August 2005
 #######################################################

# Dieses Skript sieht immer gleich aus und muss angepasst werden - v.a. die Pfade beachten!
 
# Wo die Executables leben
HELHOME=/home/lee/Downloads/Sitzung_Richard/MAGOCSOL-v0/helclim/exe
## !! Später wird er nach einem Unterordner such mit Pfad ......$YEAR$MONTH.
GZPFAD=/home/lee/Downloads/Sitzung_Richard/Beispieldaten/ausgangsdaten	 # here are the compressed Meteosat images
 #GZPFAD=$HELHOME/in # here are the compressed Meteosat images
 PFAD=/home/lee/Downloads/Sitzung_Richard/Beispieldaten/workdir   # the uncompressed images will be mv to this directory
 #CIPFAD=$HELHOME/out/CI_VIS008  # here are the cloud index output files
CIPFAD=/home/lee/Downloads/Sitzung_Richard/Beispieldaten/cloud_index_output  # the path of cloud output
GPFAD=/home/lee/Downloads/Sitzung_Richard/Beispieldaten/surface_albedo_output # here are the surface albedo images

# Anpassen um die Zeiten zu machen

cp OpenMTP_convert $PFAD
for YEAR in 2005 #1995 1997
#1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 2004
do
 #MON=11
    for  MON in 06 # 02 03 04 05 06 07 08 09 10 11 12   
    do
       mkdir  $CIPFAD/$YEAR$MON	
	# Hier MUSS local noon von dem Gebiet, für das ich kalibriere, drin sind. Die Stunden sind in UTC.
       for HH in 13 #04 05 06 07 08 09 10 11 12 14 15 16 17 18 19 20 
# 05 06 07 08 09 10 11 12 13 14 15 16 17 
	do
	    for MM in 00 # 30
	    do
	#	echo "# Year Mon hh mm " $YEAR $MON $HH $MM >> ../temp/rhomax90.out
                echo "# Year Mon hh mm " $YEAR $MON $HH $MM >> $PFAD/rhostat.out
#15 30 45
		for DAY in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31
#31
		do
		    NAME=*$YEAR$MON$DAY$HH$MM*".tar"
                    echo "looking for the files" $GZPFAD/$YEAR$MON/$NAME 
		    GZNAME=$GZPFAD/$YEAR$MON/$NAME.gz
		    BZNAME=$GZPFAD/$YEAR$MON/$NAME.bz2
                    UNAME=$GZPFAD/$YEAR$MON/$NAME
    # echo $BZNAME
    # die Bilder auspacken, falls vorhanden
		    if test -s $GZNAME
		    then
			#
			# Daten werden vorbereitet, überflüssige Kanäle rausgeschmissen, etc.
			# 
			echo auspacken $GZNAME nach $PFAD/$NAME
			gunzip -c $GZNAME > $PFAD/$NAME
                        cd $PFAD
                        tar xvf $NAME 
                        rm -f *.IR *.WV
                        mv *VISSN $YEAR$MON$DAY$HH$MM".vissn"
                        $HELHOME/OpenMTP_convert -x $YEAR$MON$DAY$HH$MM".vissn"
                        rm $YEAR$MON$DAY$HH$MM".vissn"
		    fi
		    if test -s $BZNAME
		    then
			echo auspacken $BZNAME nach $PFAD/$NAME
			bunzip2 -c $BZNAME > $PFAD/tmp.tar
                        cd $PFAD
                        tar xvf tmp.tar
                        rm tmp.tar
                        rm -f *.IR *.WV
                        mv *VISSN $YEAR$MON$DAY$HH$MM".vissn"
                        $HELHOME/OpenMTP_convert -x $YEAR$MON$DAY$HH$MM".vissn"
                        rm $YEAR$MON$DAY$HH$MM".vissn"
                    fi
                    if test -s $UNAME
                    then
                        echo $UNAME available
			cp $UNAME $PFAD/tmp.tar
			cd $PFAD
			tar xvf tmp.tar
			rm tmp.tar
			rm -f *.IR *.WV
			mv *VISSN $YEAR$MON$DAY$HH$MM".vissn"
			$HELHOME/OpenMTP_convert -x $YEAR$MON$DAY$HH$MM".vissn"
			rm $YEAR$MON$DAY$HH$MM".vissn"	
                    fi
   		done
        # Liste erzeugen. diese enthaelt die Bildnamen
		ls $YEAR$MON??$HH$MM*VISSN.xpif >> $PFAD/liste 
# orig -s 23
#$HELHOME/Cloudindex/heliosat0x8 -c $CIPFAD -g $GPFAD -l $PFAD/liste -b $PFAD -s 23 -z 0 
# nach Eumetsat konf wieder herstellen !!
#orig 3.9.2008		$HELHOME/Cloudindex/tmp300708.exe -c $CIPFAD -g $GPFAD -l $PFAD/liste -b $PFAD -s -3 -z 0
# $HELHOME/Cloudindex/tmp300708.exe -c $CIPFAD -g $GPFAD -l $PFAD/liste -b $PFAD -s -3 -z 0
# works only with version before Ma0x4 !!!
# $HELHOME/Cloudindex/heliosat-32b-Ma0x3.exe -c $CIPFAD -g $GPFAD -l $PFAD/liste -b $PFAD -s -0 -z 0
#  $HELHOME/heliosat-byte-v01-alcor.exe -c $CIPFAD/$YEAR$MON -g $GPFAD -l $PFAD/liste -b $PFAD -s -0 -z 1


# $CIPFAD/$YEAR$MON kann ich anpassen --> Das ist der Outputpfad.
# -c --> Wo kommt cloud index hin?
# -g --> Rmin (Bodenalbedo)
# -l --> Die Liste der Daten, die ich einlese
# -b --> Arbeitspfad zu den entpackten Dateien
# -s --> Der Nadir-view vom Satelliten
# -z --> Welcher Satellit ist es? War ursprünglich für Meteosat West / East gedacht. Hat mit der Pixelgröße zu tun wg. Sonnenstand, Geolokation. Ist momentan binär, 1=MFG. Auch wg. dark offset.
$HELHOME/magicsol-v0.exe -c $CIPFAD/$YEAR$MON -g $GPFAD/ -l $PFAD/liste -b $PFAD -s -0 -z 1
 
## Aifpassen - diese Pfade anlegen!
rm  $PFAD/liste
 rm $PFAD/*xpif
	    done
	    sync
# done
	done
       rm $PFAD/rhomax.out
    done
done



