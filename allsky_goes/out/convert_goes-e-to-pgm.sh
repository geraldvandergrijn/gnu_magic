# Converts all files into a PGM matching the input dimensions

# Cloud Index und Grundalbedo Fulldisk vis008
rawtopgm -bpp 1 -maxval 255 -headerskip 256 5851 6201 $1 > $1.pgm
