//---------------------------------------------------------------------------
//
//  File        :   OpenMTP_to_NetCDF.cpp
//  Description :   Export Meteosat OpenMTP format in NetCDF format
//  Project     :   Cetemps 2003
//  Author      :   Graziano Giuliani
//  Source      :   n/a
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

// Modified by Reto Stockli Dec 2008 for:
// - compilation with gcc/g++ 4.3+
// - NetCDF 3.6.3+
// - Handle ncbyte- unsigned char netcdf variable write
// - write XPIF header + data according to VCS standard
// - optionally save NetCDF or XPIF
//  
// Todo:
// - find out cres/lres for M1-M7
// - optionally specify an output filename
// - check for handling NetCDF unsigned byte data for M1-7 (currently: -127..127
// - include calibration of M7 in XPIF (already present in NetCDF)
// - check for correct West-East orientation in NetCDF, add dimension variables ... maybe
//---------------------------------------------------------------------------

#define PACKAGE_STRING "OpenMTP_convert by Reto Stockli, MeteoSwiss. Code is based on meteosatlib 0.5.5"

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <string.h>

// Unidata NetCDF

#ifdef HAS_NETCDF
#include <netcdfcpp.h>
#include <netcdf.h>
#endif

// OpenMTP format interface

#include <OpenMTP.h>

#define INSTITUTION "MeteoSwiss Climate Analysis Zurich Switzerland"
#define RECLENG (256) /* Recordlaenge fuer XPIF und PIF-Format */

//
// Creates NetCDF product
//

#ifdef HAS_NETCDF
bool NetCDFProduct( char *inpath )
{
  OpenMTP omtp;
  struct tm tmtime;
  char NcName[1024];
  char title[64];
  char reftime[64];
  int cfac, lfac, coff, loff;
  float scale;
  int wd, hg;
  int bpp;
  float *cal;
  NcVar *ivar;
  NcVar *tvar;
  NcDim *tdim;
  NcDim *ldim;
  NcDim *cdim;
  NcDim *caldim;
  ncbyte *image;
  char chname[10];

  char *fileptr1;
  char *fileptr2;
  char filename[256];
  char filepath[256];
  char filetemp[256];
  
  // first: separate filename and path, since we want to save the 
  // output file in the same path as the input file

  // extract filename from full path
  strcpy(filetemp,inpath);
  fileptr2 = filetemp;
  for (fileptr1=fileptr2;*fileptr2;fileptr2++)
    {
      if (*fileptr2 == '/') 
	{
	  fileptr1 = fileptr2+1;
	}
    }
  strcpy(filename,fileptr1);

  // remove filename from full path
  for (fileptr2=fileptr1;fileptr2<fileptr1+strlen(filename);fileptr2++)
    {
      *fileptr2 = '\0';
    }
  strcpy(filepath,filetemp);

  omtp.open( inpath );

  printf("Converting to NetCDF ... \n");

  memset(chname,0,sizeof(chname));
  if (omtp.is_ir_data( )) {
    printf("saving: IR channel \n");
    strcat(chname,"IR");
  } else if (omtp.is_wv_data( )) {
    printf("saving: WV channel \n");
    strcat(chname,"WV");
  } else if (omtp.is_vis_data( )) {
    printf("saving: VIS channel \n");
    strcat(chname,"VISSN");
  } else {
    printf("saving: Unknown channel \n");
    strcat(chname,"NA");
  }
   
  tmtime = omtp.get_datetime( );

  // Build up output NetCDF file name and open it
  sprintf(NcName,"%s%06d%04d%s%s%s%s%s",filepath,(tmtime.tm_year+1900)*10000+(tmtime.tm_mon+1)*100+tmtime.tm_mday,
	  tmtime.tm_hour*100+tmtime.tm_min,"_",omtp.get_satellite_name( ),"_",chname,".nc");

  NcFile ncf ( NcName , NcFile::Replace );
  if (! ncf.is_valid()) return false;

  // Fill arrays on creation
  ncf.set_fill(NcFile::Fill);

  // Add Global Attributes
  if (! ncf.add_att("Satellite", omtp.get_satellite_name( ))) return false;
  if (! ncf.add_att("Antenna", "Meteosat Archive")) return false;
  if (! ncf.add_att("Receiver", "Meteosat Archive")) return false;
  sprintf(reftime, "%04d-%02d-%02d %02d:%02d:00 UTC",
      tmtime.tm_year + 1900, tmtime.tm_mon + 1, tmtime.tm_mday,
      tmtime.tm_hour, tmtime.tm_min);
  if (! ncf.add_att("Time", reftime) ) return false;

  if (omtp.is_A_format( ))
  {
    if (! ncf.add_att("Area_Name", "AFormat" ) ) return false;
    if (omtp.is_ir_data( ) || omtp.is_wv_data( ))
    {
      cfac = -9102222;
      lfac = -9102222;
      coff = 1248;
      loff = 1249;
      scale = 1.0;
    }
    else if (omtp.is_vis_data( ))
    {
      cfac = -18204444;
      lfac = -18204444;
      coff = 2500;
      loff = 2500;
      scale = 1.0;
    }
    else
      throw "Unsupported AFormat image (?)\n";
  }
  else if (omtp.is_B_format( ))
  {
    if (! ncf.add_att("Area_Name", "BFormat" ) ) return false;
    cfac = -18204444;
    lfac = -18204444;
    coff = 1248;
    loff = -1118;
    scale = 2.0;
  }
  else
    {
      if (! ncf.add_att("Area_Name", "Custom Format" ) ) return false;
      if (omtp.is_ir_data( ) || omtp.is_wv_data( ))
	{
	  // add: line/pixel offsets
	  cfac = -9102222;
	  lfac = -9102222;
	  coff = 1248;
	  loff = 1249;
	  scale = 1.0;
	}
      else if (omtp.is_vis_data( ))
	{
	  cfac = -18204444;
	  lfac = -18204444;
	  coff = 2500;
	  loff = 2500;
	  scale = 1.0;
	}
      else
	throw "Unsupported Custom Format image (?)\n";
    }

  if (! ncf.add_att("Projection", "GEOS(0.0)") ) return false;
  if (! ncf.add_att("Columns", omtp.npixels( ) * (int) scale) ) return false;
  if (! ncf.add_att("Lines", omtp.nlines( ) * (int) scale) ) return false;
  if (! ncf.add_att("SampleX", scale) ) return false;
  if (! ncf.add_att("SampleY", scale) ) return false;
  if (! ncf.add_att("Column_Scale_Factor", cfac) ) return false;
  if (! ncf.add_att("Line_Scale_Factor", lfac) ) return false;
  if (! ncf.add_att("Column_Offset", coff) ) return false;
  if (! ncf.add_att("Line_Offset", loff) ) return false;
  if (! ncf.add_att("Orbit_Radius", omtp.orbit_radius( )) ) return false;
  if (! ncf.add_att("Longitude", omtp.subsatellite_point( )) ) return false;
//  if (! ncf.add_att("NortPolar", omtp.is_NortPolar( )) ) return false;
//  if (! ncf.add_att("NorthSouth", omtp.is_NorthSouth( )) ) return false;
  if (! ncf.add_att("Institution", INSTITUTION) ) return false;
  if (! ncf.add_att("Conventions", "COARDS") ) return false;
  if (! ncf.add_att("history", "Created from OpenMTP format data") )
    return false;
  sprintf(title, "OpenMTP %s product for %04d-%02d-%02d %02d:%02d:00 UTC",
            omtp.get_field_name( ), tmtime.tm_year + 1900, tmtime.tm_mon + 1,
            tmtime.tm_mday, tmtime.tm_hour, tmtime.tm_min);
  if (! ncf.add_att("title", title) ) return false;


  // Dimensions
  wd = omtp.npixels( );
  hg = omtp.nlines( );
  bpp = (int) pow(2.0, omtp.bits_per_pixel( ));

  // Define dimensions
  tdim = ncf.add_dim("time");
  if (!tdim->is_valid()) return false;
  ldim = ncf.add_dim("line", hg);
  if (!ldim->is_valid()) return false;
  cdim = ncf.add_dim("column", wd);
  if (!cdim->is_valid()) return false;
  caldim = ncf.add_dim("calibration", bpp);
  if (!caldim->is_valid()) return false;

  // Get calibration values
  cal = omtp.get_calibration( );

  // Add Calibration values
  NcVar *cvar = ncf.add_var("calibration", ncFloat, caldim);
  if (!cvar->is_valid()) return false;
  cvar->add_att("long_name", "Calibration coefficients");
  cvar->add_att("variable", omtp.get_chname( ));
  cvar->add_att("units", omtp.get_chunit( ));
  if (!cvar->put(cal, bpp)) return false;

  tvar = ncf.add_var("time", ncDouble, tdim);
  if (!tvar->is_valid()) return false;
  tvar->add_att("long_name", "Time");
  tvar->add_att("units", "seconds since 2000-01-01 00:00:00 UTC");
  double atime;
  time_t ttime;
  extern long timezone;
  ttime = mktime(&tmtime);
  atime = ttime - 946684800 - timezone;
  if (!tvar->put(&atime, 1)) return false;

  // Define channel values variable
  ivar = ncf.add_var(omtp.get_chname( ), ncByte, tdim, ldim, cdim);
  if (!ivar->is_valid()) return false;
  if (!ivar->add_att("add_offset", 0.0)) return false;
  if (!ivar->add_att("scale_factor", 1.0)) return false;

  image = new ncbyte[hg*wd];

  //image2 = omtp.get_image( );

  memcpy(&image[0],omtp.get_image( ),hg*wd);

  // Write output values
  //if (!ivar->put(omtp.get_image( ), 1, hg, wd)) return false;
  if (!ivar->put(&image[0], 1, hg, wd)) return false;

  // Close NetCDF output
  (void) ncf.close( );

  return( true );
}
#endif

//
// Creates XPIF product
//
bool XPIFProduct( char *inpath )
{

  // XPIF header structure
  struct xpif_header_structure
  { 
    /* VCS Parameters are read from XPIF-HEADER  
       - do not change order of structure elements ! */
    unsigned char  nrbytes;    /* Number of Bytes per Pixel */
    unsigned short nrrecords;  /* Number of Records per Line 1REC=256byte */
    unsigned long  nrlines;    /* Number of Lines */
    unsigned long  nrcolumns;  /* Number of Columns with Information */
    unsigned char  add_header; /* Number of additional header Records */ 
    unsigned char  pif_comp;   /* Compression function? 0=no compression */
    unsigned short pif_ahx;    /* Number of additional header records following image da
ta*/ 
    unsigned char  pif_nb;    /* Number of bits per Pixel */
    unsigned char  nav_func;   /* Navigation function code */
    long line_off;             /* Line offset */
    long col_off;              /* Column offset */
    long nav_lres;             /* Line resolution */
    long nav_cres;             /* Column resolution */
    char nav_data[95];         /* specific navigation data, set to 0 if unused */
    unsigned char dat_f;       /* greylevel function */
    char dat_nam[16];          /* data_name ascii */
    char dat_enh[16];          /* default enhancement function */
    char dat_sp[29];           /* spare set to zero */
    short dat_toff;            /* temperature offset */
    char src_obj[16];          /* image source object ascii */
    short src_year;            /* year UT of nominal image time */
    short doy;                 /* day of year of nominal image time  */
    short aq_time;             /* aquisition time = hour*100+minute(UT) of NIT */
    char src_sp[34];           /* spare set to zero */
    char src_version[8];       /* File format version  */
    /* UNI-Oldenburg parameters, not directly given in XPIF-HEADER */
    // unsigned long nrcolbig;    /* Number of Columns with black frame */
  };
  
  char xpif_filename[1024];  // output filename
  FILE *xpif_fileid;        // output file id
  const int xpif_header_length = 256;  // XPIF header length (standard VCS length is 256bytes
  unsigned char xpif_hdr[xpif_header_length];  // XPIF header (raw binary)
  struct xpif_header_structure xpif_header;  // XPIF header (as structure)

  // Image information
  int wd, hg;  // width and height
  int bpp;     // bytes per pixel
  int cfac, lfac, coff, loff;
  float scale;
  struct tm tmtime;
  char chname[10];

  char *fileptr1;
  char *fileptr2;
  char filename[256];
  char filepath[256];
  char filetemp[256];
  
  // first: separate filename and path, since we want to save the 
  // output file in the same path as the input file

  // extract filename from full path
  strcpy(filetemp,inpath);
  fileptr2 = filetemp;
  for (fileptr1=fileptr2;*fileptr2;fileptr2++)
    {
      if (*fileptr2 == '/') 
	{
	  fileptr1 = fileptr2+1;
	}
    }
  strcpy(filename,fileptr1);

  // remove filename from full path
  for (fileptr2=fileptr1;fileptr2<fileptr1+strlen(filename);fileptr2++)
    {
      *fileptr2 = '\0';
    }
  strcpy(filepath,filetemp);

  OpenMTP omtp;
  
  omtp.open(inpath);

  // Dimensions
  wd = omtp.npixels( );
  hg = omtp.nlines( );
  bpp = (int) omtp.bits_per_pixel( ) / 8;
  tmtime = omtp.get_datetime( );

  // line/column offsets etc.
  if (omtp.is_A_format( ))
  {
    if (omtp.is_ir_data( ) || omtp.is_wv_data( ))
    {
      cfac = -9102222;
      lfac = -9102222;
      coff = 1248;
      loff = 1249;
      scale = 1.0;
    }
    else if (omtp.is_vis_data( ))
    {
      cfac = -18204444;
      lfac = -18204444;
      coff = 2500;
      loff = 2500;
      scale = 1.0;
    }
    else
      throw "Unsupported AFormat image (?)\n";
  }
  else if (omtp.is_B_format( ))
  {
    cfac = -18204444;
    lfac = -18204444;
    coff = 1248;
    loff = -1118;
    scale = 2.0;
  }
  else
    {
      if (omtp.is_ir_data( ) || omtp.is_wv_data( ))
	{
	  // add: line/pixel offsets
	  cfac = -9102222;
	  lfac = -9102222;
	  coff = 1248;
	  loff = 1249;
	  scale = 1.0;
	}
      else if (omtp.is_vis_data( ))
	{
	  cfac = -18204444;
	  lfac = -18204444;
	  coff = 2500;
	  loff = 2500;
	  scale = 1.0;
	}
      else
	throw "Unsupported Custom Format image (?)\n";
    }


  printf("Converting to XPIF ... \n");

  // collect XPIF header information from OpenMTP information
  /* nrbytes      1 or 2
     nrlines      lines
     nrrecords    number of 256 byte records
     nrcolumns    columns
     nrcolbig     nrrecords*recleng/nrbytes
     add_header   0
     nav_lres     667 for MSG
     nav_cres     667 for MSG
     doy          DDD
     aq_time      HHMM  */

  xpif_header.nrbytes = bpp;
  memcpy(&xpif_hdr[0], &bpp,1);
  xpif_header.nrrecords = wd / RECLENG * bpp;
  memcpy(&xpif_hdr[1], &xpif_header.nrrecords,2);
  xpif_header.nrlines = hg;
  memcpy(&xpif_hdr[3], &xpif_header.nrlines,4);
  xpif_header.nrcolumns = wd;
  memcpy(&xpif_hdr[7], &xpif_header.nrcolumns,4);
  xpif_header.add_header = 0;
  memcpy(&xpif_hdr[11], &xpif_header.add_header,1);
  xpif_header.pif_comp = 0;
  memcpy(&xpif_hdr[12], &xpif_header.pif_comp,1);
  xpif_header.pif_ahx = 0;
  memcpy(&xpif_hdr[13], &xpif_header.pif_ahx,2);
  xpif_header.pif_nb = bpp*8; // number of bits in signal (10: MSG, 8: M1-M7)
  memcpy(&xpif_hdr[15], &xpif_header.pif_nb,1);
  xpif_header.nav_func = 42;  // ????
  memcpy(&xpif_hdr[16], &xpif_header.nav_func,1);
  xpif_header.line_off = loff;
  memcpy(&xpif_hdr[17], &xpif_header.line_off,4);
  xpif_header.col_off = coff;
  memcpy(&xpif_hdr[21], &xpif_header.col_off,4);
  xpif_header.nav_lres = -999;   // needs to be changed for M1-M7!
  memcpy(&xpif_hdr[25], &xpif_header.nav_lres,4);
  xpif_header.nav_cres = -999;   // needs to be changed for M1-M7!
  memcpy(&xpif_hdr[29], &xpif_header.nav_cres,4);
  memset(&xpif_header.nav_data,0,sizeof(xpif_header.nav_data));
  memcpy(&xpif_hdr[33], &xpif_header.nav_data[0],95);
  xpif_header.dat_f = 20;   // what is that? unused by cloudindex routine
  memcpy(&xpif_hdr[128], &xpif_header.dat_f,1);
  xpif_header.dat_toff = 0; // what is that? unused by cloudinex routine
  memcpy(&xpif_hdr[190], &xpif_header.dat_toff,2);
  xpif_header.src_year = tmtime.tm_year+1900;
  memcpy(&xpif_hdr[208], &xpif_header.src_year,2);
  xpif_header.doy = tmtime.tm_yday+1;
  memcpy(&xpif_hdr[210], &xpif_header.doy,2);
  xpif_header.aq_time = tmtime.tm_hour*100 + tmtime.tm_min;
  memcpy(&xpif_hdr[212], &xpif_header.aq_time,2);

  memset(chname,0,sizeof(chname));
  if (omtp.is_ir_data( )) {
    printf("saving: IR channel \n");
    strcat(chname,"IR");
  } else if (omtp.is_wv_data( )) {
    printf("saving: WV channel \n");
    strcat(chname,"WV");
  } else if (omtp.is_vis_data( )) {
    printf("saving: VIS channel \n");
    strcat(chname,"VISSN");
  } else {
    printf("saving: Unknown channel \n");
    strcat(chname,"NA");
  }
   
  sprintf(xpif_filename,"%s%06d%04d%s%s%s%s%s",filepath,(tmtime.tm_year+1900)*10000+(tmtime.tm_mon+1)*100+tmtime.tm_mday,
	  tmtime.tm_hour*100+tmtime.tm_min,"_",omtp.get_satellite_name( ),"_",chname,".xpif");
  if (!(xpif_fileid=fopen(xpif_filename,"wb"))) return false;
  fwrite(&xpif_hdr[0],xpif_header_length,1,xpif_fileid);    
  fwrite(omtp.get_image( ),wd*hg,1,xpif_fileid);  
  fclose(xpif_fileid);
  
  return( true );
}

/* ****************************************************************************** */
/* Reads Data Images, performs calibration and store output in NetCDF and/or XPIF */
/* ****************************************************************************** */

int main( int argc, char* argv[] )
{

  int i;
  char *cp;

  if (argc == 2) {
    if (!strcmp(argv[1], "-V"))
      {
	std::cout << argv[0] << " " << PACKAGE_STRING << std::endl;
	return 0;
      }
  }
  if (argc < 3)
  {
    std::cerr << PACKAGE_STRING << std::endl;
    std::cerr << "Syntax: " << argv[0] << " [options] <OpenMTP file>" << std::endl;
    std::cerr << "Options: " << std::endl;
    std::cerr << "-n: convert OpenMTP to NetCDF (default: false)" << std::endl;
    std::cerr << "-x: convert OpenMTP to XPIF   (default: false)" << std::endl;
    std::cerr << "Example: " << argv[0] << " -n -x OpenMTP_2004-05-01_24.VISSN " << std::endl;
    std::cerr << "creates both a NetCDF and a XPIF file from the OpenMTP input file" << std::endl;
    return 1;
  }

  for (i=1; i<argc-1; i++)  /* loop through arguments without last argument (should be the filename) */
 
   {

     cp=argv[i];
     if (*cp!='-') 
       {
	 std::cerr << "No valid option specified, use" << argv[0] << " -h for help" << std::endl;
	 return 1;
       }
     
     cp++;

     if (*cp=='n') 
       {
#ifdef HAS_NETCDF
	 cp=argv[argc-1];
	 if ((*cp=='-')) // || (*cp=="*") || (*cp=="?"))
	   {
	     std::cerr << "Invalid Filename specified, use " << argv[0] << " -h for help" << std::endl;	 
	     return 1;
	   }

	 if (!NetCDFProduct(argv[argc-1]))
	   {
	     std::cerr << "Created NetCDF product NOT OK" << std::endl;
	     return 1;
	   }
#else
	 std::cerr << "NetCDF Support was not enabled during compilation. Exiting." << std::endl;	 
	 return 1;
#endif
       }
     
     else if (*cp=='x') 
       {
	 if (!XPIFProduct(argv[argc-1]))
	   {
	     std::cerr << "Created XPIF product NOT OK" << std::endl;
	     return 1;
	   } 
       }
     else 
       {
	 std::cerr << "No valid option specified, use" << argv[0] << " -h for help" << std::endl;
	 return 1;
       }
     
   }
  
  return(0);
}
//---------------------------------------------------------------------------
