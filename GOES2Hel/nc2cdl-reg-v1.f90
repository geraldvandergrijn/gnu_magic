program nc2cdl
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Name: nc2cdl
!
! compile with
! gfortran -o my.exe -m32 pres_temp_4D_rd.f90 /usr/lib/libnetcdff.a /usr/lib/libnetcdf.a 
! -I/usr/include   -I/usr/includecompile with xlf90 -O3 
!
! PURPOSE: convert MeteoSwiss netcdf data to CDL ascii files, as 1st stage for
! netcdf files containing only 1 parameter. 
!
! LANGUAGE: Fortran 90
!
!
! PROGRAMMER: R. Mueller DWD
!
!
!  INPUT FILES:
!       
!    MeteoSwiss Heliosat netcdf hm files :
! 
!  COMMENT quick and dirty solution
!!
!
! OUTPUT FILES:
!         - 
! 
! Version: 1.0  based on mediomedia.f90, but different output format.
!               use of different regione not needed, could be skiped
!               in a later version
! 
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  use netcdf
  IMPLICIT NONE
 

  INTEGER               :: year, month, day, hour, minute, prod,ii,jj, iostat,xdim,ydim,xmdim,ymdim,rdim
  REAL                  :: lonstart,latstart,multi,min,vundef,ixp,iyp
  REAL                  :: dummy
  REAL, DIMENSION(:,:,:), ALLOCATABLE :: value
  REAL, DIMENSION(:), ALLOCATABLE ::  lonrev
  REAL, DIMENSION(:), ALLOCATABLE ::  latrev
  REAL :: blat,blon,elon,elat,dx,dy
 ! This is the name of the data file we will read.
  CHARACTER*256           :: ncfname, logfname, fname
  CHARACTER*256 ,DIMENSION(:), ALLOCATABLE ::  outsisfn, outsidfn
  CHARACTER*256 ::  matchinput='input-cdl.asc'  
  INTEGER                   :: lunit_inp=10,iu,lunit1=11,iu_lat,iu_lon,iu_msg,iu_out,iu_log
  LOGICAL             :: readable
  INTEGER  :: nlonmax=0,nlonmin=10000,nlatmax=0,nlatmin=10000,latdim=0,londim=0
  
 ! This is the name of the data file we will read.
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
 !      xdim=1150
 !      ydim=600
  

  iu=lunit_inp
  INQUIRE (UNIT=iu,OPENED=readable)
  DO WHILE ( readable )
     iu = iu+1
     INQUIRE (UNIT=iu,OPENED=readable)
  END DO
  OPEN(lunit_inp,file=trim(matchinput),form='formatted',status='OLD',iostat=iostat)
  IF(iostat /= 0) THEN
     CALL printlog('E',iostat,'matchup',trim('ERROR OPENING FILE '//matchinput))
     STOP 202
  ENDIF
  write (*,*) "used config file ", trim(matchinput)
  
   ! READ AND OPEN FILENAMES FOR INPUT AND OUTPUT FILES

  write(*,*)  "READ AND OPEN FILENAMES FOR INPUT AND OUTPUT FILES: unit",lunit_inp
  
!  READ(lunit_inp,iostat=iostat,fmt='(a)') msgfname 
  READ(lunit_inp,*) ncfname ! filename of CM-SAF data in MSG proj
     print*, 'FILE: ', ncfname
  !   READ(lunit_inp,*) msgfname
  IF(iostat /= 0) THEN
     CALL printlog('E',iostat,'matchup',trim('ERROR READING FILE '//ncfname))
     STOP 202
  END IF
  
   
  READ(lunit_inp,*) logfname ! name of output file containing the name of the logfile
  IF(iostat /= 0) THEN
     CALL printlog('E',iostat,'matchup',trim('ERROR READING FILE '//logfname))
     STOP 202
  END IF
  
  ! Proceed with reading, now the dimension and limits
  
  READ(lunit_inp,*) xdim, ydim ! x,y-dimension of the allocatable fields
  IF(iostat /= 0) THEN
     CALL printlog('E',iostat,'matchup ERROR READING XDIM ', xdim)
     STOP 202
  END IF
 
 
  READ(lunit_inp,*) rdim ! number of regions
  IF(iostat /= 0) THEN
     CALL printlog('E',iostat,'matchup Error Reading RDIM', rdim)
     STOP 202
  END IF
  ALLOCATE (outsisfn(rdim),outsidfn(rdim))
 ! write(*,*) rdim
  
   READ(lunit_inp,*) year, month, day ! date of the input/output 
  IF(iostat /= 0) THEN
     CALL printlog('E',iostat,'matchup Error Reading date', minute)
     STOP 202
  END IF
  

! vmax and vmin are currently not needed  
!  READ(lunit_inp,*) vmin ! minimal value of parameter 
!  IF(iostat /= 0) THEN
!     CALL printlog('E',iostat,'matchup ERROR READING VMIN', vmin)
!     STOP 202
!  END IF
  
!  READ(lunit_inp,*) vmax ! maximum value of parameter
!  IF(iostat /= 0) THEN
!     CALL printlog('E',iostat,'matchup ERROR READING VMIN',vmax)
!     STOP 202
!  END IF

    READ(lunit_inp,*) vundef ! undefined value for parameter
  IF(iostat /= 0) THEN
     CALL printlog('E',iostat,'matchup ERROR READING vundef',vundef)
     STOP 202
  END IF
  
  ! outfname: filename characteristic of output data for the region
  ! outpfname: name of output file containing the data 

  ! finally read the information of the grids and output files...
  DO time=1,NTIME
    READ(lunit_inp,*) outsisfn(time)
    write(*,*) "sisfname", outsisfn(time)
    READ(lunit_inp,*) outsidfn(time)
    write(*,*) "sidfname", outsidfn(time)
  ! start value (lon, lat) of the regular grids  
  IF(iostat /= 0) THEN
     CALL printlog('E',iostat,'matchup ERROR READING output information: filenames ore grid')
     STOP 202
  END IF
  END DO
  READ(lunit_inp,*) blon,elon, blat, elat, dx, dy 

  CLOSE(lunit_inp,iostat=iostat)
  IF(iostat /= 0) THEN
     CALL printlog('E',iostat,'SIS',trim('ERROR CLOSING FILE '//matchinput))
     STOP 202
  ENDIF
  !write (*,*) "filename", logfname,latfname,lonfname  
    
  write(*,*) "read config file ready"
  

     write(*,*) "I start now with the netcdf interface part"

  ! Open the file. 
  call check( nf90_open(ncfname, nf90_nowrite, ncid) )

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
!  write(*,*) lats(lat)
  end do
  do lon = 1, NLONS
!     if (lons(lon) /= START_LON + (lon - 1) * 5.0) stop 2
!  write(*,*) lons(lon)
  end do

  ! Get the varids of the pressure and temperature netCDF variables.
  call check( nf90_inq_varid(ncid, SID_NAME, sid_varid) )
  call check( nf90_inq_varid(ncid, SIS_NAME, sis_varid) )

  ! Read SIS or SID from the file.
  ! Since we know the contents of the file we know that the data
  ! arrays in this program are the correct size to hold all the data.
  call check( nf90_get_var(ncid, sid_varid, sid_in) )
  call check( nf90_get_var(ncid, sis_varid, sis_in) )

  ! Check the data. It should be the same as the data we wrote.
!!! only for testing ...
!  do lon = 1, NLONS
!     do lat = 1, NLATS
!        do time = 1, NTIME
!        write(*,*) sis_in(lon,lat,time)
!!           if (sis_in(lon, lat, time) /= SAMPLE_PRESSURE + &
!!                (lon - 1) * NLATS + (lat - 1)) stop 2
!!           if (sid_in(lon, lat, time) /= SAMPLE_TEMP + &
!!                .25 * ((lon - 1) * NLATS + (lat - 1))) stop 2
!        end do
!     end do
!  end do
  
  ! Each of the netCDF variables has a "units" attribute. Let's read
  ! them and check them.
  call check( nf90_get_att(ncid, lat_varid, UNITS, lat_units_in) )
  call check( nf90_inquire_attribute(ncid, lat_varid, UNITS, len = att_len) )
! rm string comparison seems not do work properly, strings are equal, but program stops
  write (*,*) lat_units_in(1:att_len),LAT_UNITS, att_len
!  if (lat_units_in(1:att_len) /= LAT_UNITS) stop 10

  call check( nf90_get_att(ncid, lon_varid, UNITS, lon_units_in) )
  call check( nf90_inquire_attribute(ncid, lon_varid, UNITS, len = att_len) )
   write (*,*) lon_units_in(1:att_len),LON_UNITS, att_len
!  if (lon_units_in(1:att_len) /= LON_UNITS) stop 11

  call check( nf90_get_att(ncid, sis_varid, UNITS, sis_units_in) )
  call check( nf90_inquire_attribute(ncid, sis_varid, UNITS, len = att_len) )
   write (*,*) sis_units_in(1:att_len),SIS_UNITS, att_len
!  if (sis_units_in(1:att_len) /= SIS_UNITS) stop 12

  call check( nf90_get_att(ncid, sid_varid, UNITS, sid_units_in) )
  call check( nf90_inquire_attribute(ncid, sid_varid, UNITS, len = att_len) )
!  if (sid_units_in(1:att_len) /= SID_UNITS) stop 13
    write (*,*) sid_units_in(1:att_len),SID_UNITS, att_len
  ! Close the file. This frees up any internal netCDF resources
  ! associated with the file.
  call check( nf90_close(ncid) )

  ! If we got this far, everything worked as expected. Yipee! 
  print *,"*** SUCCESS reading nc file "


  ! now allocate the fields 
  
!  ALLOCATE (field1(xdim),field2(xdim),field3(xdim)) 
!  ALLOCATE (var1(xdim,ydim))
 
  ! open the input files containing the data and lat,lon fields ------------> of
  
  iu=lunit1

!  INQUIRE (UNIT=iu,OPENED=readable,iostat=iostat)
!  !write(*,*) "verdammt", iu, readable
!  DO WHILE ( readable )
!     iu = iu+1
!     INQUIRE (UNIT=iu,OPENED=readable)
!  END DO
!  iu_log=iu
!  OPEN(UNIT=iu_log,file=logfname)
!  IF(iostat /= 0) THEN
!     CALL printlog('E',iostat,'matchup',trim('ERROR OPENING FILE '//logfname))
!     STOP 202
!  ENDIF
 
 ! open logfile

   INQUIRE (UNIT=iu,OPENED=readable,iostat=iostat)
   !write(*,*) "verdammt", iu, readable
   DO WHILE ( readable )
   iu = iu+1
   INQUIRE (UNIT=iu,OPENED=readable)
   END DO
   iu_log=iu
   OPEN(UNIT=iu_log,file=logfname)
   IF(iostat /= 0) THEN
      CALL printlog('E',iostat,'merging',trim('ERROR OPENING FILE '//logfname))
      STOP 202
   ENDIF
      
! -> only needed for regridding        
    xmdim=nint(abs((blon-elon)/dx))
    ymdim=nint(abs((blat-elat)/dy))
   
 !   ixp=nint(abs((blon-lon)/dx))
 !   iyp=nint(abs((blat-lat)/dy))

   !   write(*,*) xmdim,ymdim,ixp,iyp,lat,blat

    
     ALLOCATE  (latrev(ymdim),lonrev(xmdim),value(rdim,xmdim,ymdim))
     value=-100
   !  xmdim=0
   !  ymdim=0
 
  ! now define the match-up regions, intlon west starting point, blat=south starting point
  write(*,*) "now define the match-up regions"

 
  DO ii=1,xmdim
     lonrev(ii)=blon+ii*dx
  END DO
  
  DO jj=1,ymdim
     latrev(jj)=blat+jj*dy
  END DO
  
  !k=1
 ! write (*,*) "latrev",latrev(1), latrev(k,ymdim), lonrev(k,1),lonrev(k,xmdim)  
  write (*,*) "start to search for the near neighbour"
  
 ! DO I=1,XDIM
 !    DO J=1,YDIM
 !       ! if (latitude(i,j) .GT. (latrev(k,1)-10*dy(k)) .AND. latitude(i,j) .LT. (latrev(k,ymdim)+10*dy(k))  & 
 !       ! .AND. longitude(i,j) .GT. (lonrev(k,1)-10*dx(k)) .AND. longitude(i,j) .LT. (lonrev(k,xmdim)+10*dx(k))) then  
 !       
 !       if (latitude(i,j) .GE. (latrev(1)) .AND. latitude(i,j) .LE. (latrev(ymdim))  & 
 !            .AND. longitude(i,j) .GE. (lonrev(1)) .AND. longitude(i,j) .LE. (lonrev(xmdim))) then  
 !          
 !          !       write(*,*) "now call the neighbours"
 !          call neighbours(latrev,ymdim,latitude(i,j),jj)
 !          call neighbours(lonrev,xmdim,longitude(i,j),ii)
 !       END IF
 !    END DO
 ! END DO


  
  !! <- only for regridding ....

  !! cut out region and define the start and end index 
  
  do lon = 1, NLONS
     if ((lons(lon) .GE. blon) .AND. (lons(lon) .LE. elon) )  then
        if (lon .LE. nlonmin) then 
           nlonmin=lon
        end if
        if (lon .GE. nlonmax) then
           nlonmax=lon
        end if
        londim=londim+1
     end if
  end do
  do lat = 1, NLATS
     if ((lats(lat) .GE. blat) .AND. (lats(lat) .LE. elat) )  then
        if (lat .LE. nlatmin) then 
           nlatmin=lat
        end if
        if (lat .GE. nlatmax) then
           nlatmax=lat
        end if
        latdim=latdim+1   
     end if
  end do
         
  write(*,*) blat,elat,latdim,lats(nlatmin),lats(nlatmax),nlatmin,nlatmax
  write(*,*) blon,elon,londim,lons(nlonmin),lons(nlonmax),nlonmin,nlonmax


  do prod=1,2
     do time = 1, NTIME
         ! open the respective output file
         write(*,*) "open the output files"
         
         product:  select case (prod)
         case (1) 
            fname=outsisfn(time)
         case (2) 
            fname=outsidfn(time)
         end select product

         INQUIRE (UNIT=iu,OPENED=readable,iostat=iostat)
         !write(*,*) "verdammt", iu, readable
         DO WHILE ( readable )
            iu = iu+1
            INQUIRE (UNIT=iu,OPENED=readable)
         END DO
         iu_out=iu
         OPEN(UNIT=iu_out,file=fname)
         IF(iostat /= 0) THEN
            CALL printlog('E',iostat,'matchup',trim('ERROR OPENING FILE '//fname))
            STOP 202
         ENDIF
         


        write(iu_out,*) "netcdf tmp {" 
        write(iu_out,*)  "dimensions:"
        write(iu_out,*) "lat =", latdim, ", lon =",londim,", time = unlimited ;"
        write(iu_out,*) "variables:"
        write(iu_out,*) "float lat(lat), lon(lon);"
        write(iu_out,*) 	"float Z(lat,lon);"
        write(iu_out,*)  "// variable attributes"
        write(iu_out,*) 'lat:long_name = "latitude";'
        write(iu_out,*) 'lat:units = "degree";'
        write(iu_out,*)  'lon:long_name = "longitude";'
        write(iu_out,*) 'lon:units = "degrees";'
        write(iu_out,*) 'Z:units = "Watt";'
        write(iu_out,*) "Z:valid_range = 0., 1400.;"
        write(iu_out,*)  "data:"
        write(iu_out,*) "lon="
        do lon = nlonmin, nlonmax
       !!    if ((lons(lon) .GE. blon) .AND. (lons(lon) .LE. elon) )  then
              if (lon .EQ. nlonmax) then
                 write(iu_out,*) lons(lon),";"
              else          
                 write(iu_out,*) lons(lon),","
              end if
       !!    end if
        end do
        write(iu_out,*) "lat="
        do lat = nlatmin, nlatmax
        !!   if ((lats(lat) .GE. blat) .AND. (lats(lat) .LE. elat) ) then 
              if (lat .EQ. nlatmax) then
                 write(iu_out,*) lats(lat),";"
              else          
                 write(iu_out,*) lats(lat),","
              end if          
        !!   end if
        end do
        
        ! write product array without a 2 bixel border in order to avoid lots of undefined values 
        write(iu_out,*) "Z="
        DO lat=nlatmin,nlatmax
           DO lon=nlonmin,nlonmax
              ! WRITE(iu_out,'(F7.3,X,F7.3,X,F12.5)') lonrev(ii),latrev(jj),(value(k,ii,jj))
        !!      if ((lons(lon) .GE. blon) .AND. (lons(lon) .LE. elon) .AND. (lats(lat) .GE. blat) .AND. (lats(lat) .LE. elat) ) then  
              product2:  select case (prod)
                 case (1)
                    if(lat .EQ. nlatmax .AND. lon .EQ. nlonmax) then
                       write(iu_out,*) 0.1*sis_in(lon,lat,time),";"
                       write(iu_out,*) "}"
                     else                      
                       write(iu_out,*) 0.1*sis_in(lon,lat,time),","
                    end if
                   ! write(*,*) "in sis", iu_out, fname
                    ! Anweisungsblock fuer den Fall, dass die Antwort positiv war 
                 case (2)
                    if(lat .EQ. nlatmax .AND. lon .EQ. nlonmax) then
                       write(iu_out,*) 0.1*sid_in(lon,lat,time),";"
                        write(iu_out,*) "}"
                    else                       
                       write(iu_out,*) 0.1*sid_in(lon,lat,time),","
                    end if
                 end select product2
          !!    endif
           END DO
        END DO
        CLOSE(iu_out,iostat=iostat)
        IF(iostat /= 0) THEN
           CALL printlog('E',iostat,'matchup',trim('ERROR CLOSING FILE '//fname))
           STOP 202
        ENDIF
     END DO
  END DO
! close first output file
  

  write(*,*) "Hurra, erster fertig"
 
  
 ! write(iu_outp,*) "# value [W/m2] nearest neighbour of original MSG products"
  CLOSE(iu_log,iostat=iostat)
  IF(iostat /= 0) THEN
     CALL printlog('E',iostat,'matchup',trim('ERROR CLOSING FILE '//logfname))
     STOP 202
  ENDIF

  DEALLOCATE  (latrev,lonrev)
  DEALLOCATE (value) 
 
 ! write(*,*) "timte=",k

  
  
  DEALLOCATE (outsisfn,outsidfn, STAT=iostat )
  IF (iostat .NE. 0) THEN
     STOP 'deallocation error one dimensional fields'
  end if
  
  



contains
  subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check  

  
end program nc2cdl


!end program trafo_dirty             



!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     MODULE:         neighbours
!
!     PURPOSE:        searching for the location of x and provide the 
!                     nearest neighbour of x within the discret x-array xa(n)
!
!     LANGUAGE:       FORTRAN 90
!
!     CALLING MODE:   subroutine call
!
!     FORTRAN CALLING INSTRUCTION:
!     CALL borders(xa,n,x,i)
!
!     PARAMETER DESCRIPTION:
!     INTEGER  klo, khi, k   -> counters
!     REAL xa(n)   -> array containing the discrete x values   
!     REAL x       -> value for which the borders xa(k) and
!                    xa(k-1) are scanned
!     INTEGER n    -> Dimension of the array xa() 
!     INTEGER k    -> k=khi upper counter of the border
!                     the interval
!
!     PROGRAMMER:    based on numerical recipes, modifications 
!                    Richard Mueller/Rainer Hollman - DWD
!
!     VERSION: 1.0         
!     DATE:    24.11.2004
!
!
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE neighbours(xa,n,x,k)
  
  IMPLICIT NONE
  
  INTEGER ::  n
  REAL    ::  x, xa(n)
  INTEGER ::  klo, khi,k
  klo=1
  khi=n
!  if ( xa(1) .GT. xa(n) ) then
!  write(*,*) "Hallo"
!  khi=1
!  klo=n
!  end if

1 IF (khi-klo.GT.1) THEN
     k=(khi+klo)/2
     IF (xa(k).GT. x)THEN
        khi=k
     ELSE
        klo=k
     ENDIF
     GOTO 1
  ENDIF
  IF ((khi-1) .NE. klo) THEN 
      print*, 'WARNING problems in borders subroutine khi-1 NEQ klo'
  ENDIF
  if (abs(xa(klo)-x) .GT. abs(xa(khi)-x)) then
  k=khi
 ! write(*,*) "will use khi",xa(khi),xa(klo),x  
  else
  k=klo
 !  write(*,*) "will use klo",xa(khi),xa(klo),x
  endif
! only for testing ------------------------------------->
!  if ((k .NE. khi) .AND. (k .NE. klo )) then
!   print*, 'WARNING problems in borders k is not eq klo or khi'
!   print*, 'k,khi,klo', k,khi,klo 
! endif
! <---------------------------------------------------    
  return
END SUBROUTINE neighbours



!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     MODULE:         printlog
!
!     PURPOSE:        printlog writes logging information to stdout
!
!     LANGUAGE:       FORTRAN 90
! 
!     CALLING MODE:   subroutine call
!
!     FORTRAN CALLING INSTRUCTION: 
!     CALL PRINTLOG (type, code,pge, message)
!
!     PARAMETER DESCRIPTION:
!     Variable       Dim/Length  Type  I/O    Meaning
!
!       type          variable     C    I     type of notification
!       code            1          I    I     numerical code of the message
!       pge           varable      C    I     record number
!       message       variable     C    I     length of datfield to read
!
!
!     PROGRAMMER:    Annegret Gratzki  /DWD         DATE:    08.11.2001
!
!
!     MODIFICATIONS: 
!
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
   SUBROUTINE printlog(type, code,pge, message)
   
   IMPLICIT NONE
   
   CHARACTER*(*)          :: type
   CHARACTER*(*)          :: pge
   CHARACTER*(*)          :: message
   INTEGER                :: code
  
   CHARACTER                   :: ydate*8, ytime*10, yzone*5 ! current UTC time   
   INTEGER ,DIMENSION (8)      :: inow_time   ! current UTC date and time
   
   ! get current date
   CALL DATE_AND_TIME (ydate, ytime, yzone, inow_time)
   
   WRITE(*,1000) type,code,pge,ydate,inow_time(5),inow_time(6),inow_time(7),pge,message
   1000 FORMAT(a,1x,i7.6,1x,a,1x,a8,1x,i2.2,':',i2.2,':',i2.2,1x,a,1x,a)
   
   RETURN
   END subroutine printlog
   
