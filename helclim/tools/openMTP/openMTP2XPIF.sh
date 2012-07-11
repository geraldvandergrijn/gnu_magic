year=1994
month=02
slot=26
path=/cmsaf/cmsaf-rad3/data/MSGraw/openMTP/import/m5/199402-VIS/
for day in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31
do
OpenMTP_convert -x $path"OpenMTP_"$year"-"$month"-"$day"_"$slot"_512345_1_1_1.VISSN"
done 