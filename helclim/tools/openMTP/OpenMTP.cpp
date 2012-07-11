//-----------------------------------------------------------------------------
//
//  File        : OpenMTP.cpp
//  Description : Meteosat OpenMTP format interface
//  Project     : Meteosatlib
//  Author      : Graziano Giuliani
//  References  : Eumetsat EUM-FG-1 Format Guide: OpenMTP format
//                Revision 2.1 April 2000
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
//
//-----------------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <time.h>
#include <stdlib.h>
#include <OpenMTP_ascii_header.h>
#include <OpenMTP_binary_header.h>
#include <OpenMTP_image.h>
#include <OpenMTP.h>

OpenMTP::OpenMTP( ) { nominal_orbit_radius = 42164.0; }
OpenMTP::~OpenMTP( ) { }

void OpenMTP::open( char *inname )
{
  std::ifstream infile(inname, (ios::binary | ios::in));
  if (infile.fail())
  {
    cerr << "Cannot open input OpenMTP file " << inname << endl;
    throw;
  }
  read(infile);
  return;
}

void OpenMTP::read( std::ifstream &file )
{
  ascii_header.read(file);
  binary_header.read(file);
  image.read(file, binary_header);
  return;
}

struct tm &OpenMTP::get_datetime( )
{
  // 2008/12/10 Reto Stockli: added yday calculation
  static struct tm itm;
  int tmp;
  bool leap;
  int monthsum[13] = {0,31,59,90,120,151,181,212,243,273,304,334,365};

  itm.tm_year = binary_header.year( ) - 1900;  // year a.d. - 1900

  // check if it is a leap year
  leap = ( (((binary_header.year() % 4) == 0) && 
	    ((binary_header.year() % 100) != 0)) || 
	   ((binary_header.year() % 400) == 0) ); 

  tmp = binary_header.date( );
  tmp = tmp-((tmp/10000)*10000);
  itm.tm_mon  = tmp/100 - 1;  // 0-11
  itm.tm_mday = tmp-(itm.tm_mon+1)*100; // 1-31
  if (leap && ((itm.tm_mon+1)>2)) {
    /* leap year with march - december */
    itm.tm_yday = monthsum[itm.tm_mon] + 1 + itm.tm_mday - 1; // [0-365] days since january 1 [time.h]
  } else {
    /* regular year */
    itm.tm_yday = monthsum[itm.tm_mon] + itm.tm_mday - 1; // [0-364] days since january 1 [time.h]
  }
  tmp = binary_header.time( );
  itm.tm_hour = tmp/100;  // hour since midnight 0-23
  itm.tm_min = tmp-(itm.tm_hour*100);  // minute after the hour 0-59

  return itm;
}

char * OpenMTP::get_datatimestr( )
{
  struct tm itm = get_datetime( );
  static char tmp[32];
  strftime(tmp, 32, "%Y-%m-%d %H:%M:%S %Z", &itm);
  return tmp;
}

char * OpenMTP::get_satellite_code( )
{
  return binary_header.satellite_name( );
}

char * OpenMTP::get_satellite_name( )
{
  return binary_header.satellite_name( );

// Reto Stockli: we cannot return a string constant with a char* function
// please adjust for current gcc. Fix should work right now.

  /*
  if(!strcmp(binary_header.satellite_name( ), "M7")) return "Meteosat7";
  if(!strcmp(binary_header.satellite_name( ), "M6")) return "Meteosat6";
  if(!strcmp(binary_header.satellite_name( ), "M5")) return "Meteosat5";
  if(!strcmp(binary_header.satellite_name( ), "M4")) return "Meteosat4";
  if(!strcmp(binary_header.satellite_name( ), "M3")) return "Meteosat3";
  if(!strcmp(binary_header.satellite_name( ), "M2")) return "Meteosat2";
  if(!strcmp(binary_header.satellite_name( ), "M1")) return "Meteosat1";
  throw "Can only manage METEOSAT OpenMTP files\n"; 
  */
}

char * OpenMTP::get_timestr( )
{
  static char timestr[14];
  int tmp;
  tmp = binary_header.date( );
  tmp = tmp-((tmp/10000)*10000);

  sprintf(timestr, "%04d%04d_%04d", binary_header.year( ),
     tmp, binary_header.time( ));
  return timestr;
}

bool OpenMTP::is_A_format( )
{
  if ((binary_header.is_ir_data( ) || binary_header.is_wv_data( )) &&
     (binary_header.nlines( ) == 2500 && binary_header.npixels( ) == 2500))
    return true;
  else if (binary_header.is_visible_composite( ) &&
     (binary_header.nlines( ) == 5000 && binary_header.npixels( ) == 5000))
    return true;
  return false;
}

bool OpenMTP::is_B_format( )
{
  if ((binary_header.is_ir_data( ) || binary_header.is_wv_data( ))    &&
     (binary_header.nlines( )==625 && binary_header.npixels( )==1250) &&
     (binary_header.first_line( )==1810))
    return true;
  else if (binary_header.is_visible_composite( )                        &&
      (binary_header.nlines( )==1250 && binary_header.npixels( )==2500) &&
      (binary_header.first_line( )==3620))
    return true;
  return false;
}

bool OpenMTP::is_subwindow( )
{
  return (!is_A_format( ) || !is_B_format( ));
}

int OpenMTP::nlines( ) { return binary_header.nlines( ); }
int OpenMTP::npixels( ) { return binary_header.npixels( ); }
int OpenMTP::bits_per_pixel( ) { return OpenMTP_bpp; }

bool OpenMTP::is_ir_data( ) { return binary_header.is_ir_data( ); }
bool OpenMTP::is_wv_data( ) { return binary_header.is_wv_data( ); }
bool OpenMTP::is_vis_data( ) { return binary_header.is_vis_data( ); }

float OpenMTP::orbit_radius( ) { return nominal_orbit_radius; }
float OpenMTP::subsatellite_point( )
{
  return binary_header.subsatellite_point( );
}

float *OpenMTP::get_calibration( ) { return image.cal( ); }
unsigned char *OpenMTP::get_image( ) { return image.data( ); }

const char *OpenMTP::get_chname( )
{
  if (binary_header.is_vis_data( )) return "VIS";
  else if (binary_header.is_ir_data( )) return "IR";
  else if (binary_header.is_wv_data( )) return "WV";
  else
    throw "Invalid Channel in OpenMTP format\n";
}

const char *OpenMTP::get_chunit( )
{
  if (binary_header.is_vis_data( )) return "%";
  else if (binary_header.is_ir_data( )) return "K";
  else if (binary_header.is_wv_data( )) return "K";
  else
    throw "Invalid Channel in OpenMTP format\n";
}

const char *OpenMTP::get_field_name( )
{
  static char tmp[32];
  strncpy(tmp, binary_header.field_name( )+1, 32);
  return tmp;
}

#ifdef TESTME

int main(int argc, char *argv[])
{
  OpenMTP omtp;

  std::ifstream infile(argv[1], (std::ios::binary | std::ios::in));
  if (infile.fail())
  {
    std::cerr << "Cannot open input OpenMTP file " << argv[1] << std::endl;
    throw;
  }

  omtp.read(infile);

  std::cout << omtp;

  infile.close( );
  return 0;
}

#endif
