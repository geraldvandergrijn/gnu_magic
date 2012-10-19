moyear=062003
for day in 152 153 154 155 156 157 158 159 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 181 
do 
infile="goes09.2003."$day".012514.BAND_01.nc"
outcdl="goes09.2003."$day".0125.cdl"
out="goes09.2003."$day".0125.out"
echo "\"/cmsaf/cmsaf-rad2/data/GOESraw/gms-062003_low_res/$infile\"" > input-cdl.asc
echo "\"/cmsaf/cmsaf-rad2/data/GOESraw/gms-062003_low_res/info.log\"" >> input-cdl.asc
echo "5208 2706" >> input-cdl.asc
echo "1" >> input-cdl.asc
echo "26 255 2003" >> input-cdl.asc
echo "-1000" >> input-cdl.asc
echo "\"/cmsaf/cmsaf-rad2/data/GOESraw/regnc/$outcdl\"" >>  input-cdl.asc
echo  "\"/cmsaf/cmsaf-rad2/data/GOESraw/regnc/$out\"" >>  input-cdl.asc
echo "95 180 -62 62  0.04  0.04" >>  input-cdl.asc
echo "160 0.0" >> input-cdl.asc
GOESnc2Hel.exe
done