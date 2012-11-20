# Set month and year
moyear=062003
for day in {152..181} # 153 154 155 156 157 158 159 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 181 
do 
	# Set hours
	for hour in 01 # 00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23
	do
		infile="goes09.2003."$day"."$hour"2514.BAND_01.nc"
		outcdl="goes09.2003."$day"."$hour"25.cdl"
		out="goes09.2003."$day"."$hour"25.out"
		echo "\"/home/lee/Downloads/goes/062003_full_res/$infile\"" > input-cdl.asc
		echo "\"/home/lee/Downloads/goes/062003_full_res/info.log\"" >> input-cdl.asc
		# Set xdim and ydim (can be retrieved with ncdump -h)
		echo "20836 10827" >> input-cdl.asc
		# Set number of regions, do not change
		echo "1" >> input-cdl.asc
		# Scan time, highest possible value (dependent on bit length) and year of scan
		echo "26 65535 2003" >> input-cdl.asc
		# No data value, do not change
		echo "-1000" >> input-cdl.asc
		# Set output paths
		echo "\"/home/lee/Desktop/$outcdl\"" >> input-cdl.asc
		echo  "\"/home/lee/Desktop/$out\"" >>  input-cdl.asc
		# Set dimensions: xmin, xmax, ymin, ymax, and resolution in lat/long
		# Note: The program can't deal with values going over the 180° border
		# For < 180° use GMS coverage, > -180° uses GOES-W coverage
		echo "95 180 -62 62  0.02 0.02" >>  input-cdl.asc
		# These values are only defined for testing
		echo "160 0.0" >> input-cdl.asc
		/home/lee/Documents/eclipse-workspace/GOES2Hel/GOESnc2Hel.exe
	done
done
