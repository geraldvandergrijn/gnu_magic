# Set month and year
moyear=062003
# Set days to be processed
for day in {152..181}
do
	# Set hours
	for hour in 01 #01 02 03 04 05 06 07 08 09 11 12 13 14 15 16 17 18 19 20 21 22 23
	do
		# Set input paths
		infile="goes09.2003."$day"."$hour"2514.BAND_01.nc"
		outcdl="goes09.2003."$day"."$hour"25.cdl"
		out="goes09.2003."$day"."$hour"25.out"
		echo "\"/home/lee/Downloads/goes/raw/$infile\"" > input-cdl.asc
		echo "\"/home/lee/Downloads/goes/raw/info.log\"" >> input-cdl.asc
		# Set xdim and ydim (can be retrieved with ncdump -h)
		echo "20812 10819" >> input-cdl.asc
		# Set number of regions, do not change
		echo "1" >> input-cdl.asc
		# Scan time, highest possible value (dependent on bit length) and year of scan
		echo "26 65535 2003" >> input-cdl.asc
		# No data value, do not change
		echo "-1000" >> input-cdl.asc
		# Set output paths
		echo "\"/home/lee/Downloads/goes/raw_processed/$outcdl\"" >>  input-cdl.asc
		echo  "\"/home/lee/Downloads/goes/raw_processed/$out\"" >>  input-cdl.asc
		# Set dimensions: xmin, xmax, ymin, ymax, and resolution in lat/long
		# Note: The program can't deal with values going over the 180° border
		# For <180° use GMS coverage, >-180° uses GOES-W coverage
		echo "95 180 -62 62  0.01  0.01" >>  input-cdl.asc
		# These values are only defined for testing
		echo "160 0.0" >> input-cdl.asc
		GOESnc2Hel.exe
		# Remove ASCII CDL file
		rm "/home/lee/Downloads/goes/062003_regridded_decompressed/goes09.2003."$day"."$hour"25.cdl"
	done
done
