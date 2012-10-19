! The program is based on the example program provided 
! as part of the netCDF package.

! The oroginal program relies to 
! Copyright 2006 University Corporation for Atmospheric Research/Unidata.
! See COPYRIGHT file for conditions of use.

! This program reads SIS, SID and n from the netcdf files provided by MeteoSwiss
! Programmer of modifications: R. Mueller

! Full documentation of the netCDF Fortran 90 API can be found at:
! http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90

! $Id: sfc_pres_temp_rd.f90,v 1.7 2006/12/09 18:44:58 russ Exp $

program read_ms_hm
  use netcdf
  implicit none

  ! This is the name of the data file we will read.
  !character (len = *), parameter :: FILE_NAME = "delme.nc"
  character (len = *), parameter :: FILE_NAME = "/cmsaf/cmsaf-hcp1/rad_tmp/netcdf/2000/m7_20000126_hm_0x03.nc"
  ! character (len = *), parameter :: FILE_NAME = "pres_temp_4D.nc"
  integer :: ncid

  ! We are reading 3D data, a NLATS x NLONS lat-lon grid, 
  ! Each parameter covers 24 hours 
  integer, parameter :: NDIMS = 3
  integer, parameter :: NLATS = 4667, NLONS = 4667
  integer, parameter :: NTIME = 24 
  character (len = *), parameter :: LAT_NAME = "lat"
  character (len = *), parameter :: LON_NAME = "lon"
  character (len = *), parameter :: TIME_NAME = "time"
  integer :: lat_dimid, lon_dimid, time_dimid

  ! For the lat lon coordinate netCDF variables.
  real :: lats(NLATS), lons(NLONS)
  integer :: lat_varid, lon_varid, time_varid

  ! We will read solar surface irradiance and direct irradiande. 
  character (len = *), parameter :: SIS_NAME = "SIS"
  character (len = *), parameter :: SID_NAME = "SID"
  integer :: sis_varid, sid_varid
  integer :: dimids(NDIMS)

  ! To check the units attributes.
  character (len = *), parameter :: UNITS = "units"
  character (len = *), parameter :: SIS_UNITS = "W m-2"
  character (len = *), parameter :: SID_UNITS = "W m-2"
  character (len = *), parameter :: LAT_UNITS = "degrees_north"
  character (len = *), parameter :: LON_UNITS = "degrees_east"
  integer, parameter :: MAX_ATT_LEN = 80
  integer :: att_len
  character*(MAX_ATT_LEN) :: sis_units_in, sid_units_in
  character*(MAX_ATT_LEN) :: lat_units_in, lon_units_in
   character*(MAX_ATT_LEN) :: TMP_UNITS

  ! Read the data into these arrays.
  integer(kind=2) :: sis_in(NLONS, NLATS, NTIME), sid_in(NLONS, NLATS, NTIME)

  ! These are used to calculate the values we expect to find.
  real, parameter :: START_LAT = 25.0, START_LON = -125.0
  real, parameter :: SAMPLE_PRESSURE = 900.0
  real, parameter :: SAMPLE_TEMP = 9.0

  ! We will learn about the data file and store results in these
  ! program variables.
  integer :: ndims_in, nvars_in, ngatts_in, unlimdimid_in

  ! Loop indices
  integer :: lat, lon, time

  write(*,*) "I start now"

  ! Open the file. 
  call check( nf90_open(FILE_NAME, nf90_nowrite, ncid) )

  ! There are a number of inquiry functions in netCDF which can be
  ! used to learn about an unknown netCDF file. NF90_INQ tells how many
  ! netCDF variables, dimensions, and global attributes are in the
  ! file; also the dimension id of the unlimited dimension, if there
  ! is one.
  call check( nf90_inquire(ncid, ndims_in, nvars_in, ngatts_in, unlimdimid_in) )

  write(*,*) ndims_in, nvars_in, ngatts_in, unlimdimid_in

  ! In this case we know that there are 2 netCDF dimensions, 4 netCDF
  ! variables, no global attributes, and no unlimited dimension.
!  if (ndims_in /= 3 .or. nvars_in /= 8 .or. ngatts_in /= 6 &
!       .or. unlimdimid_in /= 3) stop "1st"

  write(*,*) ncid
  ! Get the varids of the latitude and longitude coordinate variables.
  call check( nf90_inq_varid(ncid, LAT_NAME, lat_varid) )
  call check( nf90_inq_varid(ncid, LON_NAME, lon_varid) )
  write(*,*) "after lat lon check"
  call check( nf90_inq_varid(ncid, TIME_NAME, time_varid) )
  ! Read the latitude and longitude data.
  call check( nf90_get_var(ncid, lat_varid, lats) )
  call check( nf90_get_var(ncid, lon_varid, lons) )

 

  ! Check to make sure we got what we expected.
  do lat = 1, NLATS
!     if (lats(lat) /= START_LAT + (lat - 1) * 5.0) stop 2
  write(*,*) lats(lat)
  end do
  do lon = 1, NLONS
!     if (lons(lon) /= START_LON + (lon - 1) * 5.0) stop 2
  write(*,*) lons(lon)
  end do

  ! Get the varids of the pressure and temperature netCDF variables.
  call check( nf90_inq_varid(ncid, SID_NAME, sis_varid) )
  call check( nf90_inq_varid(ncid, SIS_NAME, sid_varid) )

  ! Read the surface pressure and temperature data from the file.
  ! Since we know the contents of the file we know that the data
  ! arrays in this program are the correct size to hold all the data.
  call check( nf90_get_var(ncid, sid_varid, sid_in) )
  call check( nf90_get_var(ncid, sis_varid, sis_in) )

  ! Check the data. It should be the same as the data we wrote.
  do lon = 1, NLONS
     do lat = 1, NLATS
        do time = 1, NTIME
!        write(*,*) sis_in(lon,lat,time)
!           if (sis_in(lon, lat, time) /= SAMPLE_PRESSURE + &
!                (lon - 1) * NLATS + (lat - 1)) stop 2
!           if (sid_in(lon, lat, time) /= SAMPLE_TEMP + &
!                .25 * ((lon - 1) * NLATS + (lat - 1))) stop 2
        end do
     end do
  end do
  
  ! Each of the netCDF variables has a "units" attribute. Let's read
  ! them and check them.
  write (*,*)
  call check( nf90_get_att(ncid, lat_varid, UNITS, lat_units_in) )
  call check( nf90_inquire_attribute(ncid, lat_varid, UNITS, len = att_len) )
  write (*,*) lat_units_in(1:att_len),"Hallo",LAT_UNITS,"Hallo",lat_units_in(1:att_len),"Hallo", att_len
  TMP_UNITS=lat_units_in(1:att_len)
  IF (TMP_UNITS .EQ. LAT_UNITS) WRITE (*,*) 'The strings are equal! '
  if (LGE(lat_units_in(1:att_len),LAT_UNITS)) stop 10

  call check( nf90_get_att(ncid, lon_varid, UNITS, lon_units_in) )
  call check( nf90_inquire_attribute(ncid, lon_varid, UNITS, len = att_len) )
  write (*,*) LON_UNITS," Hallo ",lon_units_in(1:att_len)," Hallo ", att_len
 if (trim(lon_units_in(1:att_len)) /= trim(LON_UNITS)) stop 11

  call check( nf90_get_att(ncid, sis_varid, UNITS, sis_units_in) )
  call check( nf90_inquire_attribute(ncid, sis_varid, UNITS, len = att_len) )
  if (sis_units_in(1:att_len) /= SIS_UNITS) stop 12

  call check( nf90_get_att(ncid, sid_varid, UNITS, sid_units_in) )
  call check( nf90_inquire_attribute(ncid, sid_varid, UNITS, len = att_len) )
  if (sid_units_in(1:att_len) /= SID_UNITS) stop 13

  ! Close the file. This frees up any internal netCDF resources
  ! associated with the file.
  call check( nf90_close(ncid) )

  ! If we got this far, everything worked as expected. Yipee! 
  print *,"*** SUCCESS reading example file sfc_pres_temp.nc!"

contains
  subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check  
end program read_ms_hm

