program GOESnc2Hel
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Name: GOESnc2Hel
!
!  compile with 
!  gfortran -m32 -g GOESnc2Hel.f90 /usr/lib/libnetcdff.a /usr/lib/libnetcdf.a -I/usr/include -o GOESnc2Hel.exe  
!
! PURPOSE: convert netcdf GOES data to binary files in regular lat lon, as 1st stage for
! netcdf files . 
!
! LANGUAGE: Fortran 90
!
!
! PROGRAMMER: R. Mueller DWD
!
!
!  INPUT FILES:
!       
!    GOES VIS data together with latitude and longitude files 
!    in GOES projection
!
!   OUTPUT file:
!   GOES VIS data on a regular latitude longitude grid.
! 
!  COMMENT quick and dirty solution
!!
!
! OUTPUT FILES:
!         - 
! 
! Version: 1.0  ! 
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  use netcdf       ! use netcdf in order to be able to read/write netcdf files
  IMPLICIT NONE
  
  INTEGER               :: i,ii,j,jj,ilat,ilon,iostat,in,jn,iref,jref,kk,ll,k
  INTEGER               :: xdim,ydim,recno,reclength,xmdim,ymdim,kr,ixp,iyp
  INTEGER               :: imax=0,imin=30000,jmax=0,jmin=30000,rdim
  INTEGER               :: year, month, day, hour, minute
  REAL                  :: lonstart,latstart,multi,min,vundef
  REAL                  :: dummy
  REAl, DIMENSION(:,:), ALLOCATABLE  :: delta,delta1
 ! INTEGER*2, DIMENSION(:,:), 
  byte, DIMENSION(:,:), ALLOCATABLE :: image
  REAL, DIMENSION(:,:,:), ALLOCATABLE :: value,var1
  REAL, DIMENSION(:), ALLOCATABLE ::  lonrev
  REAL, DIMENSION(:), ALLOCATABLE ::  latrev
  REAL, DIMENSION(:), ALLOCATABLE ::  blon,blat,elon,elat,dx,dy,lat,lon
  REAL, DIMENSION(:,:), ALLOCATABLE :: longitude,latitude
  CHARACTER*256           :: goesfname, logfname
  CHARACTER*256 ,DIMENSION(:), ALLOCATABLE ::  outfname,outbinfname
  CHARACTER*256 ::  matchinput='input-cdl.asc'  
  INTEGER                   :: lunit_inp=10,iu,lunit1=11,iu_lat,iu_lon,iu_msg,iu_out,iu_outp,iu_log
  INTEGER :: ncid
  LOGICAL             :: readable
  byte :: bla
  
! needed for netcdf files
  
  character (len = *), parameter :: LAT_NAME = "lat"
  character (len = *), parameter :: LON_NAME = "lon"
  character (len = *), parameter :: TIME_NAME = "imageTime"
  character (len = *), parameter :: VAR_NAME = "data"
  character (len = *), parameter :: DATE_NAME = "imageDate"


  integer :: lat_varid, lon_varid, time_varid, var_varid,  date_varid, date
  integer :: idummy, time
  character (len = 7) cdummy

! type def (c struct) needed for definition of the header 
  type :: header
     integer*4   :: nrlines
     integer*4   :: nrcolumns
     real*4        :: latbegin
     real*4        :: lonbegin   
     real*4        :: dlat
     real*4        :: dlon
     real*4        :: latend
     real*4        :: lonend
     real*4   :: lat_off=999
     real*4   :: lon_off=999
     integer*4   :: aq_time
     integer*4   :: doy 
     integer*4   :: scandur ! duration of the scan
     integer*4   :: vmax    ! maximum value of the data
     integer*4   :: dum1 ! dummy value
     integer*4   :: dum2 ! dummy value 
     integer*4    :: fillv(48)  ! needed to get 256 byte header

!     byte, DIMENSION(:,:), ALLOCATABLE :: image 
  end type header
  type (header) ::  head
  DO i=1,48
  head%fillv(i)=999
  END DO
!  write(*,*) sizeof(i) 

 ! head%nrlines
  !      xdim=1150
  !      ydim=600


  ! call check( nf90_open(goesfname, nf90_nowrite, ncid) )

  write(*,*) "start of program"

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
  
!  READ(lunit_inp,iostat=iostat,fmt='(a)') goesfname 
  READ(lunit_inp,*) goesfname ! filename of CM-SAF data in MSG proj
  !   print*, 'FILE: ', goesfname
  !   READ(lunit_inp,*) goesfname
  IF(iostat /= 0) THEN
     CALL printlog('E',iostat,'matchup',trim('ERROR READING FILE '//goesfname))
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
   ALLOCATE (blon(rdim), blat(rdim), dx(rdim), dy(rdim), lon(rdim), lat(rdim)) 
   ALLOCATE (outfname(rdim),outbinfname(rdim),elat(rdim),elon(rdim))
 ! write(*,*) rdim
  
   READ(lunit_inp,*) head%scandur, head%vmax, head%dum1 ! date of the input/output 
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
  ! outbinfname: name of output file containing the data 

  ! finally read the information of the grids and output files...
  DO K=1,RDIM
    READ(lunit_inp,*) outfname(k)
    write(*,*) outfname(k)
    READ(lunit_inp,*) outbinfname(k) 
    write(*,*) outbinfname(k)
    READ(lunit_inp,*) blon(k),elon(k), blat(k), elat(k), dx(k), dy(k) 
    READ(lunit_inp,*) lon(k), lat(k)
  ! start value (lon, lat) of the regular grids  
  IF(iostat /= 0) THEN
     CALL printlog('E',iostat,'matchup ERROR READING output information: filenames ore grid' &
   ,lat(k))
     STOP 202
  END IF
  END DO
  

  
  CLOSE(lunit_inp,iostat=iostat)
  IF(iostat /= 0) THEN
     CALL printlog('E',iostat,'SIS',trim('ERROR CLOSING FILE '//matchinput))
     STOP 202
  ENDIF
  !write (*,*) "filename", logfname,latfname,lonfname  
    
  write(*,*) "read config file ready"
  
 
 ! reclength=xdim*4
 ! write(*,*) "here I am",reclength 

  ! now allocate the fields 
  
  ALLOCATE (longitude(xdim,ydim),latitude(xdim,ydim),var1(xdim,ydim,1))

 
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
      
  
  INQUIRE (UNIT=iu,OPENED=readable,iostat=iostat)
  !write(*,*) "verdammt", iu, readable
  DO WHILE ( readable )
     iu = iu+1
     INQUIRE (UNIT=iu,OPENED=readable)
  END DO

! replace this read command by netcdf read ...
 
! Open the file. 
  call check( nf90_open(goesfname, nf90_nowrite, ncid) )

 ! check and get the variables
 ! Get the varids of the latitude and longitude coordinate variables.
  call check( nf90_inq_varid(ncid, LAT_NAME, lat_varid) )
  call check( nf90_inq_varid(ncid, LON_NAME, lon_varid) )
  write(*,*) "after lat lon check"
  call check( nf90_inq_varid(ncid, TIME_NAME, time_varid) )
  ! Read the latitude and longitude data.
  call check( nf90_get_var(ncid, lat_varid, latitude) )
  call check( nf90_get_var(ncid, lon_varid, longitude) )
  ! check and read the data of the visible channel
  call check( nf90_inq_varid(ncid, VAR_NAME, var_varid) )  
  call check( nf90_get_var(ncid, var_varid, var1) )
  call check ( nf90_inq_varid(ncid, DATE_NAME, date_varid) )
  call check ( nf90_get_var(ncid, date_varid, date) )
  call check ( nf90_inq_varid(ncid, TIME_NAME, time_varid) ) 
  call check ( nf90_get_var(ncid, time_varid, time) )

  call check( nf90_close(ncid) )


  ! <---------------of
  
  ! now read the original CM-SAF data
  write(*,*) "now read the original CM-SAF data"
  
 
  ! close all files relvant for reading the msg data
  write(*,*) "start to  close all files relvant for reading the msg data"

 ! allocate the fields needed for the next step and deallocate the fields not 
 ! needed any longer

    
  
    !initialise the fields
  
  DO K=1,RDIM
         
    xmdim=nint(abs((blon(k)-elon(k))/dx(k)))+5
    ymdim=nint(abs((blat(k)-elat(k))/dy(k)))+5
   
    ixp=nint(abs((blon(k)-lon(k))/dx(k)))+3
    iyp=nint(abs((blat(k)-lat(k))/dy(k)))+3

      write(*,*) xmdim,ymdim,ixp,iyp,lat(k),blat(k)

    
     ALLOCATE  (latrev(ymdim),lonrev(xmdim),value(rdim,xmdim,ymdim))
     ALLOCATE (delta(xmdim,ymdim),delta1(xmdim,ymdim))
     value=-100.0
     delta=1000.0
     delta1=1000.0
   !  xmdim=0
   !  ymdim=0

    

  
  ! now define the match-up regions, intlon west starting point, blat=south starting point
  write(*,*) "now define the match-up regions"
 

 
  
  DO ii=1,xmdim
     lonrev(ii)=blon(k)+(ii-3)*dx(k)
  END DO
  
  DO jj=1,ymdim
     latrev(jj)=blat(k)+(jj-3)*dy(k)
  END DO
  
  !k=1
 ! write (*,*) "latrev",latrev(1), latrev(k,ymdim), lonrev(k,1),lonrev(k,xmdim)  
  write (*,*) "start to search for the near neighbour"
  
  DO I=1,XDIM
     DO J=1,YDIM
  ! if (latitude(i,j) .GT. (latrev(k,1)-10*dy(k)) .AND. latitude(i,j) .LT. (latrev(k,ymdim)+10*dy(k))  & 
  ! .AND. longitude(i,j) .GT. (lonrev(k,1)-10*dx(k)) .AND. longitude(i,j) .LT. (lonrev(k,xmdim)+10*dx(k))) then  
          
 if (latitude(i,j) .GE. (latrev(1)) .AND. latitude(i,j) .LE. (latrev(ymdim))  & 
   .AND. longitude(i,j) .GE. (lonrev(1)) .AND. longitude(i,j) .LE. (lonrev(xmdim))) then  

!       write(*,*) "now call the neighbours"
         call neighbours(latrev,ymdim,latitude(i,j),jj)
         call neighbours(lonrev,xmdim,longitude(i,j),ii)
        ! value(k,ii,jj)=var1(i,j)
        ! write (8,*) latrev(k,:), lonrev(k,:)
             !  write(*,*) ii,jj, var1(i,j), longitude(i,j),latitude(i,j) 
           if (ii .GT. 2 .AND. ii .LT. (xmdim-1) .AND. jj .GT. 2 .AND. jj .LT. (ymdim-1)) then
               DO kk=ii-2,ii+2
                  DO ll=jj-2,jj+2
               delta(kk,ll)=sqrt(abs(longitude(i,j)-lonrev(kk))**2+abs(latitude(i,j)-latrev(ll))**2)
                      !write(8,*) "outif", ii,jj,delta1(ii,jj),delta(ii,jj), var1(i,j)
                      if (delta(kk,ll) .LT. delta1(kk,ll)) THEN
                         delta1(kk,ll)=delta(kk,ll)
                         value(k,kk,ll)=var1(i,j,1)
                !     write(8,*) "inif", ii,jj,delta1(ii,jj),delta(ii,jj), var1(i,j), value(k,kk,ll)
                      end if
                   END DO
                END DO
            else 
            delta(ii,jj)=sqrt(abs(longitude(i,j)-lonrev(ii))**2+abs(latitude(i,j)-latrev(jj))**2)
                !write(*,*) "outif", ii,jj,delta1(ii,jj),delta(ii,jj), sis(i,j)
                if (delta(ii,jj) .LT. delta1(ii,jj)) THEN
                   delta1(ii,jj)=delta(ii,jj)
                   value(k,ii,jj)=var1(i,j,1)
             ! write(8,*) "edge", ii,jj,delta1(ii,jj),delta(ii,jj), var1(i,j), value(k,ii,jj)
                end if
             end if
              !  write (*,*) "latitude, longitude", latitude,longitude    
              !if (sis(i,j) .LT. 1368 .AND. sis(i,j) .GE. 0) then
              !              WRITE(19,*) longitude(i,j), latitude(i,j), latitude(i,j), longitude(i,j), sis(i,j)
              !ENDIF
           END IF
        END DO
     END DO
  

! open the respective output file
  write(*,*) "open the output file"
    
  INQUIRE (UNIT=iu,OPENED=readable,iostat=iostat)
  !write(*,*) "verdammt", iu, readable
  DO WHILE ( readable )
     iu = iu+1
     INQUIRE (UNIT=iu,OPENED=readable)
  END DO
  iu_out=iu
  OPEN(UNIT=iu_out,file=outfname(k))
  IF(iostat /= 0) THEN
     CALL printlog('E',iostat,'matchup',trim('ERROR OPENING FILE '//outfname(k)))
     STOP 202
  ENDIF
  
 ! write header and data to first output  file


  write(iu_out,*) "netcdf tmp {" 
  write(iu_out,*)  "dimensions:"
  write(iu_out,*) "lat =", ymdim-4, ", lon =",xmdim-4,", time = unlimited ;"
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
  write(iu_out,*) "lat="
  DO jj=3,ymdim-2
   if (jj .EQ. (ymdim-2)) then
      WRITE(iu_out,*) latrev(jj),";"
   else
      WRITE(iu_out,*) latrev(jj),","
   endif
  END DO

  write(iu_out,*) "lon="
  DO ii=3,xmdim-2
   if (ii .EQ. (xmdim-2)) then
      WRITE(iu_out,*) lonrev(ii),";"
   else
      WRITE(iu_out,*) lonrev(ii),","
   endif
  END DO
   
   

  ! write product array without a 2 bixel border in order to avoid lots of undefined values 

  

  write(iu_out,*) "Z="
  DO jj=3,ymdim-2
     DO ii=3,xmdim-2
       ! WRITE(iu_out,'(F7.3,X,F7.3,X,F12.5)') lonrev(ii),latrev(jj),(value(k,ii,jj))
        if (ii .EQ. (xmdim-2) .AND. jj .EQ. (ymdim-2)) then
           WRITE(iu_out,*) value(k,ii,jj),";"
            WRITE(iu_out,*) "}"
        else
           WRITE(iu_out,*) value(k,ii,jj),","
        endif
    END DO
  END DO


! close first output file
  
  CLOSE(iu_out,iostat=iostat)
  IF(iostat /= 0) THEN
     CALL printlog('E',iostat,'matchup',trim('ERROR CLOSING FILE '//outfname(k)))
     STOP 202
  ENDIF
  

  write(*,*) "Hurra, erster fertig"
  ! open second output file
 !!!! WRITEIT

 ! define the header values

  head%nrlines=ymdim-4  
  head%nrcolumns=xmdim-4
  head%latbegin=blat(1)
  head%lonbegin=blon(1)
  head%dlat=dy(1)
  head%dlon=dx(1)
  head%latend=elat(1)
  head%lonend=elon(1)
  head%doy=date-head%dum1*1000
  write(*,*) "head=", head%doy, date, head%dum1*1000
  head%aq_time=time/100
  write(*,*) "aq_time=", head%aq_time  
  
  !print *, SIZEOF(dummy)
  ! reclength=xdim ! * byte per var
  ALLOCATE (image(xmdim-4,ymdim-4))
  write(*,*) head

  reclength=head%nrlines*head%nrcolumns+256
  OPEN (UNIT=99,file=outbinfname(1),form='unformatted',access='direct',recl=reclength) 
  DO i=3,xmdim-2
     DO j=3,ymdim-2   
! image(i,j)=transfer(bla,var1(i,j,1))  
    ! image(i,j)=varvalue(1,i,j)1(i,j,1) 
    image(i-2,j-2)=(255/head%vmax)*value(1,i,j)
!transfer(bla,value(1,i,j))
!(255/head%vmax)*value(1,i,j) ! in order to provide always values in the byte range 0-256
    END DO
  END DO
 


 WRITE(99,rec=1,iostat=iostat) head ,image(:,:)
 !DO i=1,xdim
 !   WRITE (99,rec=i+1,iostat=iostat) (image(i,j),j=1,ydim)
 !END DO
 
 !write(99,*) "hallo" 

 ! INQUIRE (UNIT=iu,OPENED=readable,iostat=iostat)
 ! !write(*,*) "verdammt", iu, readable
 ! DO WHILE ( readable )
 !    iu = iu+1
 !    INQUIRE (UNIT=iu,OPENED=readable)
 ! END DO
 ! iu_outp=iu
 ! OPEN(UNIT=iu_outp,file=outbinfname(k),POSITION='APPEND')
  IF(iostat /= 0) THEN
     CALL printlog('E',iostat,'matchup',trim('ERROR OPENING FILE '//outbinfname(k)))
     STOP 202
  ENDIF
  CLOSE(iu_outp,iostat=iostat)
  IF(iostat /= 0) THEN
     CALL printlog('E',iostat,'matchup',trim('ERROR CLOSING FILE '//outbinfname(k)))
     STOP 202
  ENDIF

 
  
 ! write(iu_outp,*) "# value [W/m2] nearest neighbour of original MSG products"
  CLOSE(iu_log,iostat=iostat)
  IF(iostat /= 0) THEN
     CALL printlog('E',iostat,'matchup',trim('ERROR CLOSING FILE '//logfname))
     STOP 202
  ENDIF

  DEALLOCATE  (latrev,lonrev)
  DEALLOCATE (delta,delta1,value) 
  DEALLOCATE (image) 

  write(*,*) "k=",k
   
END DO
  
  
  DEALLOCATE (outfname,outbinfname, STAT=iostat )
  IF (iostat .NE. 0) THEN
     STOP 'deallocation error one dimensional fields'
  end if
  
  
  DEALLOCATE (longitude, latitude, var1, STAT=iostat )
  IF (iostat .NE. 0) THEN
     STOP 'deallocation error two dimensional fields'
  end if
    

contains
  subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check  


  
end program GOESnc2Hel


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
   
