program netcdfchecker

! ---------------------------------
! --> Purpose: Check NetCDF file
!
! Marco Rohrer, Sep 2015

! v0.2:   Open NetCDF file
! v0.3:   Save NetCDF file, minor first file manipulations
! v0.4:   Handle missing values and invert axis, handle packed data
! v0.5:   Adapt output routine (only Z200), handle wrong time axis, leap year and corrupt files
! v0.6:   Added several measures to check netcdf files
! v0.7.1: Added mm=-1 switch for yearly files
! v0.7.x: Added several expansions, for example support of 3D variables (x,y,t) 
! v0.8:   Major revisions to modularize and tidy up  everything
! v0.8.4: NetCDF read finished
! v0.8.5: First version of NetCDF write finished
! v0.8.6: NetCDF check part implemented, subroutine invar, handle non-existant data
! v0.8.7: Implement corrgrid (part of NetCDF check, bugfix for invertlat
! v0.8.8: Implement fort.fails and netcdf.check now outputs wrong timesteps, dataset type must now be given
! v0.8.9: Bug fixes and style
! v0.8.10:Implement loading of two variables (only from two different variables so far) and bugfixes
! v0.8.11:Memory leaks...bug fixes and existing netcdf file will be overwritten
! v0.8.12:Include fort.fails routine in main program, minor bug fixes
! v0.8.14:Version clarity, corrected some false alarms, correct handling of not existing ccc400 files, bugfixes
! v0.8.15:Fix bug which sets ccc400 non GPH/SLP values not to NA, better documentation, code cleanup, leap day bug fix
! v0.8.16:Bug fixes (leap years and missing files)and now able to unpack and correct packed lon/lat/lev/time axis, now able to use ccc400 NO_VOLC data
! v0.8.17:Bug fixes (leap years in 20CR)
! v0.8.18:Bug fixes (ERAint) and make ERA20C compatible
! v0.8.19:Bug fixes (ERAint), new data (netcdf) directly downloaded from apps.ecmwf.int do not work
! v0.8.20:Bug fixes (20CR) and JRA55 can now be read
!         --> Attention: Use /home/marco/prog NetCDF library, /usr/lib does not work!

! Still missing: - Arithmetics of two levels of the same variable, e.g. z1000-z500)
!                - Lower memory footprint
!                - Modularize input/output for future programs
!                - Check for more datatypes
!                - OpenMP support
!                - Code cleanup
! Requirements:  - A functional NetCDF library 
!                - F2003 compatible compiler (e.g. gfortran4.8, should work since gfortran4.3 or gfortran4.4... gcc.gnu.org/wiki/Fortran2003Status)
!                - datafails.txt for variables other than slp/geopotential
USE netcdf

IMPLICIT NONE 


! Variable declarations
TYPE attributes                             ! Container that holds all important attributes from the CF infile
    character(32)                           :: units, long_name, standard_name
    character(32)                           :: calendar, positive
    real(4)                                 :: FillVal,scalef,offset
END TYPE attributes ! type to store attributes

! Read variable string thing
CHARACTER(32)                               :: invar1,invar2                    ! variable names (if two vars, both used)
CHARACTER(32)                               :: mode                             ! one/two variables from input?
CHARACTER(1)                                :: volc                             ! for ccc: VOLC or NO_VOLC? n=NO_VOLC
INTEGER                                     :: inlev1,inlev2,inlev              ! levels of input data --> inlev inlev* used in loop
INTEGER                                     :: invarn                           ! Counter variable for nvarin
INTEGER                                     :: nvarin=1

! Fields
REAL, dimension(:,:,:,:), allocatable       :: arr,arrin1,arrin2                ! 4D array for input (may be two input files)
REAL, dimension(:,:,:),allocatable          :: arr3out,arr3out2,arr3temp        ! 3D arrays for output (no level left!)

TYPE (attributes)                           :: attlon,attlat,atttime
TYPE (attributes)                           :: attlev,attdat                    ! Containers for attributes of different variables

! Netcdf stuff
CHARACTER(100)                              :: infile, outfile, intype
CHARACTER(128)                              :: invarnames                       ! Name of input variable (used in loop)
CHARACTER(24)                               :: lont,latt,levt,timet             ! Dimension names 
CHARACTER(24),DIMENSION(:),ALLOCATABLE      :: varnames,dimnames                ! Vector with all variable / dimension names
INTEGER,DIMENSION(:),ALLOCATABLE            :: varndims                         ! Vector with number of dims per variable
INTEGER                                     :: ncid, ncido                      ! NetCDF IDs for infile and outfile (ncido)              
INTEGER                                     :: ndims,nvars                      ! Number of dimensions and variables in input file
INTEGER                                     :: nglobatt                         ! Number of global attributes
INTEGER                                     :: svarndims                        ! Number of dimensions in selected variable
INTEGER                                     :: nx=-1,ny=-1,ntime=-1             ! Number of gridpoints (longitude,latitude,time)
INTEGER                                     :: nlev=-1,nmem=-1                  ! Number of levels or members (-1 means no dimension)
INTEGER                                     :: LonDimID=-1,LatDimID=-1          ! Dimension IDs of NetCDF file
INTEGER                                     :: VerDimID=-1,TimeDimID=-1
INTEGER                                     :: LonVarID=-1,LatVarID=-1          ! Variable IDs of NetCDF file
INTEGER                                     :: VerVarID=-1,TimeVarID=-1
INTEGER                                     :: DatVarID=-1
REAL                                        :: missval=-999999                  ! Output missing value
LOGICAL                                     :: packed=.true.                    ! Data packed/not packed
LOGICAL                                     :: axpacked=.false.                 ! Dimension variable packed/not packed

REAL(8),dimension(:),allocatable            :: times,time,lons,lats,levels,mems ! 1-D arrays to hold time,lons,lats, levels (time:in;times:out)
! checker stuff
INTEGER                                     :: year, mn, member                 ! year, month (and member)  of input data
CHARACTER(4)                                :: iyear, imn, ileap                ! character to read arguments
! misc variables

INTEGER                 :: dimn            ! Counter variable for ndims (Number of dimensions)
INTEGER                 :: varn            ! Counter variable for nvars (Number of variables)
INTEGER                 :: natts           ! Number of attributes of a variable
CHARACTER(32)           :: version         ! version number
INTEGER                 :: istat           ! status variable
LOGICAL                 :: exfile=.true.   ! file exists
LOGICAL                 :: leap            ! Leap year switch
LOGICAL                 :: convgp=.false.  ! Convert geopotential -> gph
LOGICAL                 :: convaprx=.false.! Convert aprl kg/m**2s -> mm/timestep
LOGICAL                 :: ens=.false.     ! Level dimension is actually ens member
LOGICAL                 :: dbg=.false.     ! If .true. debugging information will be printed (A LOT!)
! Switches and misc

! netcdf check stuff
INTEGER                                 :: ii, tt, lvl          ! Count variables
INTEGER                                 :: validtmstps,take,pass! Valid tmstps, first valid tmstp,pass
INTEGER                                 :: levelid=-1           ! Select level
INTEGER                                 :: days, ntimes         ! Number of days and timesteps
REAL                                    :: mmax, mmin, mean     ! Get mean 
INTEGER,DIMENSION(:,:),ALLOCATABLE      :: vert2                ! Control variable for vertical axis
LOGICAL                                 :: existfort            ! T/F to check if file with data flaws exists
LOGICAL                                 :: leapy                ! Checks whether it is a leap year
LOGICAL                                 :: checkup              ! Full checkup?
!LOGICAL                                 :: cgrid=.false.
INTEGER,DIMENSION(:,:),ALLOCATABLE      :: corrections          ! Array with all the corrections to be made in a file
INTEGER                                 :: nrows                ! Number of rows in corrections
CHARACTER(100)                          :: outfileraw           ! Save outfilename throughout checking routine
INTEGER,DIMENSION(:),ALLOCATABLE        :: lendims              ! Length of dimensions (feed into nx,ny,...)

!intype="ccc"
version='v0.8.20'
IF(iargc()==1) THEN
    CALL getarg(1,infile)
    SELECT CASE (infile)
        CASE("-h","--help","-H","--Help")
            PRINT*, "This program semi-automatically reads and corrects NetCDF files. &
                   &It was designed to read monthly ccc400 files, but should be able to &
                   &read also yearly 20CR and ERA-interim files. Other files may or may &
                   &not work. CF-convention NetCDF files only! At the moment the only output option is a 3D field. So you &
                   &must select a height level or add/subtract several variables before &
                   &outputting. See details below."
            PRINT*,"NOTE: A 20CR file has the following structure: (/lat,lon,ens_mem,time/)&
                   & and is analyzed differently. The routine will output each member seper&
                   &ately in order to keep the I/O footprint low. However you need at least&
                   & 12 GB RAM to get reasonable performance (and presumably no crashes)!&
                   & The output files are named according to your outfilename (second&
                   & argument, see below), but suffixed with the ensemble number, i.e.:&
                   & z200_1999.nc will give z200_1999.nc.1 ... z200_1999.nc.56"
            PRINT*, "Usage: 7 input arguments are needed:"
            PRINT*, "1: infile --> Path of (full) input file "
            PRINT*, "2: outfile --> Path and/or Name auf output file"
            PRINT*, "3: year --> Year of data, time axis will be overwritten because ccc400&
                  & sometimes provides wrong time axis"
            PRINT*, "4: month --> Month of data, if the input is a yearly file, month=-1"
            PRINT*, "5: leapyears? --> True/False statement (T/F), if true leap years are &
                  & accounted in dataset, otherwise (F) we do not account for leap years &
                  & and they are removed from the dataset"
            PRINT*, "6: Invarnames -->Variable name (maybe check first with ncdump/ncview.) &
                  & A specific level must be selected if necessary (4D-variables), separated&
                  & by a comma. Eg. geopoth,20000. The level should be given in Pa. At the &
                  & moment it is possible to add or subtract two different variables. In &
                  & order to add convective and large-scale precipitation: 'aprc+aprl'. &
                  & This features does not work currently for the same variable (e.g. &
                  & 'geopoth,50000-geopoth,100000'." 
            PRINT*, "7: Intype --> Define type of input data: ccc for ccc400, 20cr for 20cr,&
                  & jra for JRA55, e20c for ERA-20c &
                  & and eraint for ERA-interim. IF you input another dataset, try with eraint, &
                  & as least assumptions are made for this dataset. However this is experimental!&
                  & In order to use the datafails.txt, the member is required after ccc for ccc400,&
                  & i.e. 'ccc001'. In case you access NO_VOLC data, append a 'n' after the ensemble&
                  & member, i.e. 'ccc001n' (Or just use ccc to circumvent the datafails file (not tested)!"
            PRINT*, "Version: ",version
            PRINT*, "Written by M. Rohrer (marco.rohrer@giub.unibe.ch),last changed 20150723"
        CASE("--Version","-V","--version")
            PRINT*, "Version: ",version
            PRINT*, "Written by M. Rohrer (marco.rohrer@giub.unibe.ch),last changed 20150723"
    END SELECT
    STOP
ENDIF


IF(iargc()/=7) THEN
    PRINT*,'Error: Wrong number of Input Arguments given. Correct form is:'
    PRINT*,'Infilename Outfilename Inyear Inmonth Leap?[T/F] Invariable(s) Datatype'
    PRINT*,'Example: ./netcdf.check in.nc out.nc 1999 01 F slp ccc001'
    STOP
ENDIF ! IF Iargc
CALL getarg(5,ileap)   ; READ (ileap,'(L1)') leap

! Read arguments
CALL getarg(1,infile)
CALL getarg(2,outfile)
CALL getarg(3,iyear)   ; READ (iyear,'(I4)') year
CALL getarg(4,imn)     ; READ (imn,'(I2)') mn
CALL getarg(6,invarnames)
CALL getarg(7,intype)
! Read variable names
CALL invar(invarnames,invar1,invar2,inlev1,inlev2,mode)
IF (dbg) THEN
    WRITE(*,*) "Input variables read:"
    PRINT*,"Infile:  ", trim(infile)
    PRINT*,"Outfile: ", trim(outfile)
    PRINT*,"Inyear:  ", trim(iyear)
    PRINT*,"Inmonth: ", trim(imn)
    PRINT*,"Leap?:   ", trim(ileap)
    PRINT*,"Invars:  ", trim(invarnames)
    IF (mode=="onevar") THEN
        PRINT*,"One variable:  ", trim(invar1),inlev1
    ELSE
        PRINT*,"Two variables: ", trim(invar1),inlev1,trim(invar2),inlev2
    ENDIF
    PRINT*,"=============================================================="
ENDIF ! IF dbg

If (mode/="onevar") THEN
    nvarin=2
ENDIF ! IF mode

IF (intype(1:3)=="ccc") THEN
    READ (intype(4:6),'(I3)') member
    IF (len(trim(intype))==7) THEN
        volc="n"
    ELSE
        volc="v"
    ENDIF  
    intype="ccc"
    IF (member==0) THEN
        WRITE(*,*) "WARNING: Type 'ccc' but member not indicated! datafails will not work!"
    ELSE IF (member.LT.1.OR.(member.GT.30.AND.member.NE.101.AND.member.NE.102.AND.member.NE.103)) THEN
        WRITE(*,*) "WARNING: Type 'ccc' but member is not existing! datafails will not work!"
    ENDIF ! PRINT warning message when member not given...
ENDIF ! IF intype

! Input netcdf
! Open NetCDF -> get ncid
istat = NF90_OPEN(infile, nf90_nowrite, ncid)

!istat = NF90_OPEN(infile, nf90_nowrite, ncid)
IF(istat/=nf90_noErr) THEN 
    WRITE(*,*) "Open NetCDF file22:", NF90_STRERROR(istat)
    exfile=.false.
    IF (intype(1:3).NE."ccc") STOP
ENDIF

IF(exfile) THEN
    ! Inquire information -> get number of dimensions, number of variables,
    ! number of global attributes and the ID of the unlimited axis (time)
    istat = NF90_INQUIRE(ncid,ndims,nvars,nglobatt,TimeDimID)
    IF (istat/=nf90_noErr) PRINT*,"Inq info:", NF90_STRERROR(istat)
    
    ALLOCATE(dimnames(ndims),lendims(ndims),varnames(nvars),varndims(nvars))
    IF(dbg) print*,"ndims: ",ndims,"nvars",nvars
    DO invarn=1,nvarin
        IF (invarn==2.AND.invar1==invar2) CYCLE
        IF (invarn==1) invarnames=invar1
        IF (invarn==2) invarnames=invar2
        IF (invarn==1) inlev=inlev1
        IF (invarn==2) inlev=inlev2

        DO dimn=1,ndims
        
            istat = NF90_INQUIRE_DIMENSION(ncid,dimn,dimnames(dimn),lendims(dimn))
            IF (istat/=nf90_noErr) PRINT*,"Inquire Dimension:", NF90_STRERROR(istat)
            ! Select case seems not work as we compare the select variable to another variable (invarnames)
            IF (dimnames(dimn)=="lon" .OR. dimnames(dimn)=="long" .OR.& 
&               dimnames(dimn)=="longitude" .OR. dimnames(dimn)=="g0_lon_3") THEN
                nx=lendims(dimn)
                lont=dimnames(dimn)
                istat = NF90_INQ_DIMID(ncid,dimnames(dimn),LonDimID)
            ELSE IF (dimnames(dimn)=="lat" .OR. dimnames(dimn)=="g0_lat_2" .OR.&
&                    dimnames(dimn)=="latitude" ) THEN
                ny=lendims(dimn)
                latt=dimnames(dimn)
                LatDimID=dimn
                istat = NF90_INQ_DIMID(ncid,dimnames(dimn),LatDimID)
            ELSE IF (dimnames(dimn)=="time" .OR. dimnames(dimn)=="initial_time0_hours") THEN
                ntime=lendims(dimn)
                timet=dimnames(dimn)
                istat = NF90_INQ_DIMID(ncid,dimnames(dimn),TimeDimID)
            ELSE IF (dimnames(dimn)=="lev" .OR. dimnames(dimn)== "level" .OR. dimnames(dimn)=="lv_ISBL1") THEN
                nlev=lendims(dimn)
                levt=dimnames(dimn)
                VerDimID=dimn
                istat = NF90_INQ_DIMID(ncid,dimnames(dimn),VerDimID)
            ELSE IF (dimnames(dimn)=="ensemble_member" .OR. dimnames(dimn)=="ens") THEN
                nlev=lendims(dimn)
                levt=dimnames(dimn)
                VerDimID=dimn
                istat = NF90_INQ_DIMID(ncid,dimnames(dimn),VerDimID)
                IF (istat/=nf90_noErr) PRINT*,"Inquire Dimension:", NF90_STRERROR(istat)
                ens=.true.
            ELSE IF (dimnames(dimn)==invarnames) THEN
                istat = NF90_INQ_VARID(ncid,dimnames(dimn),DatVarID)
            ENDIF  ! IF dimnames(dimn)               
        ENDDO ! DO dimn

        IF(dbg) THEN ; DO dimn=1,ndims ; print*,"Dimension: ",dimnames(dimn),lendims(dimn) ; ENDDO ; ENDIF
        ! Inquire variables, get their names and dimensions
        IF(dbg) print*,"==============================================================="
        DO varn=1,nvars
            istat = NF90_INQUIRE_VARIABLE(ncid,varn,name=varnames(varn),ndims=varndims(varn),natts=natts)
            IF (varnames(varn)=="lon" .OR. varnames(varn)=="long" .OR. &
&               varnames(varn)=="longitude" .OR. varnames(varn)=="g0_lon_3") THEN
                istat = NF90_INQ_VARID(ncid,varnames(varn),LonVarID)
                CALL getattributes(ncid,LonVarID,attlon)
            ELSE IF (varnames(varn)=="lat" .OR. varnames(varn)=="g0_lat_2" .OR. &
&                    varnames(varn)=="latitude") THEN
                istat = NF90_INQ_VARID(ncid,varnames(varn),LatVarID)
                CALL getattributes(ncid,LatVarID,attlat)
            ELSE IF (varnames(varn)=="time" .OR. varnames(varn)=="initial_time0_hours") THEN
                istat = NF90_INQ_VARID(ncid,varnames(varn),TimeVarID)
                CALL getattributes(ncid,TimeVarID,atttime)
            ELSE IF (varnames(varn)=="lev" .OR. varnames(varn)== "level" .OR. varnames(varn)=="lv_ISBL1") THEN
                istat = NF90_INQ_VARID(ncid,varnames(varn),VerVarID)
                !print*,"L305: VerVarID=",VerVarID
                CALL getattributes(ncid,VerVarID,attlev)
            ELSE IF (varnames(varn)=="ensemble_member" .OR. varnames(varn)=="ens") THEN
                istat = NF90_INQ_VARID(ncid,varnames(varn),VerVarID)
                !print*,"L308: VerVarID=",VerVarID
                CALL getattributes(ncid,VerVarID,attlev)
            ELSE IF (varnames(varn)==invarnames) THEN
                istat = NF90_INQ_VARID(ncid,varnames(varn),DatVarID)
                CALL getattributes(ncid,DatVarID,attdat)
                svarndims=varndims(varn)
            ENDIF ! IF varnames matches...IF(dbg) THEN ; DO dimn=1,ndims ; print*,"Dimension: ",dimnames(dimn),lendims(dimn) ; ENDDO ; ENDIF


            IF(dbg) write(*,'(I3,A12,A6,I2,A9,I3,I3,A9,I3,I3,A9,I3,I3,A9,I3,I3,A7,I3)') varn,trim(varnames(varn)),&
                    &"ndims",varndims(varn),"LoDimVar",LonDimID,LonVarID,"LaDimVar",LatDimID,LatVarID,&
                    &"VeDimVar",VerDimID,VerVarID,"TiDimVar",TimeDimID,TimeVarID,"DatVar",DatVarID
        ENDDO ! DO varn
        IF(dbg) print*,"svarndims",svarndims
        IF(DatVarID==-1) THEN ; print*,"Error! Variable not found. Program aborts!" ; CALL abort ; ENDIF
        IF(dbg) write(*,*) '======================================================================='
        IF(dbg) write(*,'(A76,5(A16),3(E12.3))')' ATT Lon:[unit|long_name|std_name|positive|calendar|_FillV|scalef|offset]:',attlon
        IF(dbg) write(*,'(A76,5(A16),3(E12.3))')' ATT Lat:[unit|long_name|std_name|positive|calendar|_FillV|scalef|offset]:',attlat
        IF(dbg) write(*,'(A76,5(A16),3(E12.3))')' ATT Tim:[unit|long_name|std_name|positive|calendar|_FillV|scalef|offset]:',atttime
        IF(dbg) write(*,'(A76,5(A16),3(E12.3))')' ATT Lev:[unit|long_name|std_name|positive|calendar|_FillV|scalef|offset]:',attlev
        IF(dbg) write(*,'(A76,5(A16),3(E12.3))')' ATT Dat:[unit|long_name|std_name|positive|calendar|_FillV|scalef|offset]:',attdat

        ! Check whether conversion required
        SELECT CASE (attdat%units)
            CASE ('m**2s**-2','m**2 s**-2')
                convgp=.true.
                IF(dbg) print*,"Data needs conversion: Geopotential to geopotential height"
            CASE ('kg/m**2s')
                IF(dbg) print*,"Data needs conversion: PREC rate to PREC amount per timestep"
                convaprx=.true.
            CASE DEFAULT
        END SELECT

        IF(dbg) print*,"==================================================================="
        IF(dbg) WRITE(*,'(A13,I3,3(A6,I3),A8,I3)') "VarIDs: Lon: ", LonVarID,"Lat: ",LatVarID,"Lev: ",VerVarID,"Time: ",TimeVarID

        IF (invarn==1) THEN
             !allocate(mems(1))
             DO dimn=1,nvars
                  CALL readaxis(ncid,dimn,LonVarID,LatVarID,VerVarID,TimeVarID,time,lons,lats,levels,mems,&
                                &nx,ny,ntime,nlev,attlon,attlat,attlev,atttime,axpacked)
             ENDDO ! DO dimn
             IF(.NOT.ALLOCATED(mems)) allocate(mems(1))
             IF(.NOT.ALLOCATED(levels)) allocate(levels(1))
        ENDIF ! IF invarn==1 --> Axis needs to be read just once (allocate and overhead...)
        IF(dbg) THEN
            write(*,'(A3,I4,A4,I4,A5,I4,A5,I4,A4,I4)') "nx=",size(lons),"ny=",size(lats),"nme=",size(mems),&
                                                   "nle=",size(levels),"nt=",size(time)
            print*,"lon",lons(1),lons(UBOUND(lons))
            print*,"lat",lats(1),lats(UBOUND(lats))
            if (VerVarID.gt.0) print*,"lev",levels(1),levels(UBOUND(levels))
            print*,"tim",time(1),time(UBOUND(time))
            Print*,"========================================================================"
        ENDIF ! IF dbg
        ! Read data
        IF(dbg) print*,"getdata:",ncid,DatVarID,nx,ny,nlev,ntime,svarndims,dbg,packed,convgp,convaprx
        CALL getdata(ncid,DatVarID,nx,ny,nlev,ntime,svarndims,dbg,packed,convgp,convaprx,arr)
        
        ! Do arithmetics if nvarin==2
        IF (dbg) WRITE(*,*) "Check whether two variables are read in and do arithmetics if necessary"
        IF (nvarin==2.AND.invarn==1.AND.invar1/=invar2) THEN
            IF (dbg) WRITE(*,*) "Allocate arrin1,arrin2 to perform the arithmetics" 
            ALLOCATE(arrin1(nx,ny,nlev,ntime),arrin2(nx,ny,nlev,ntime))
            arrin1=arr
        ELSE IF (nvarin==2.AND.invarn==2.AND.invar1/=invar2) THEN
            IF (dbg) WRITE(*,*) "Do arithmetics on invar1 and invar2:"
            arrin2=arr
            arr=0
            CALL arith(arrin1,arrin2,mode,arr)
            DEALLOCATE(arrin1,arrin2)
            
            IF ((invar1=="aprl".OR.invar2=="aprl").AND.(invar1=="aprc".OR.invar2=="aprc") ) THEN
                invarnames="apr"
            ENDIF ! IF total prec
        ENDIF ! IF nvarin
    ENDDO ! DO invarn
    istat =nf90_close(ncid)
    !!! ===================== Input Done ============================================
    IF(dbg) print*,'=====================Input Done==================================='

ENDIF ! IF exfile


outfileraw=outfile
!!! ====================== NetCDF check =======================================================
!==============================================================================================
IF(dbg) print*,"L0391: array is:",arr(1,1,:,1)
IF(dbg.AND.svarndims.EQ.3) print*,"L0392: array dims: nx=",size(arr,1),"ny=",size(arr,2),"nl=",size(arr,3),"nt=",size(arr,4)

! Read fort.fails --> If existing. 
! Otherwise first produce it with the datafails.txt file --> R script provided
existfort=.FALSE.
CALL createfortfail(year,mn,member,existfort,nrows)
IF (dbg) print*,"L359: existfort: ",existfort,"invarnames: ",invarnames
IF ( existfort .AND. intype=="ccc" .AND. volc/="n" ) THEN
    IF (.NOT.( invarnames=="geopoth".OR.invarnames=="slp" )) THEN
        CALL fortfail(corrections,nrows)
        944 FORMAT (A10,I3,A32,I4,A7,I2)
        PRINT*,invarnames
        IF (nrows.EQ.1) WRITE(*,944) "There is  ", nrows, " correction  to be made in year ", year, " month ", mn
        IF (nrows.GT.1) WRITE(*,944) "There are 1", nrows, " corrections to be made in year ", year, " month ", mn
    ELSE
        nrows=-1
    ENDIF ! if invarnames
ELSE
    nrows=-1
ENDIF ! IF existfort

! Correct lon/lat if data was packed
IF (axpacked) THEN
    !DEALLOCATE(lons,lats)
    CALL dummyaxis(lons,lats,nx,ny,intype,attdat,attlon,attlat,atttime,invarnames,missval)
ENDIF ! IF axpacked

! Correctly prepare missing file and get rid of leap year
IF(dbg) print*,'====== Check whether file exists and correct leap year ==========='
attdat%fillval=-999999
leapy=.false.
IF (((MOD(year,4).lt.1).AND.(MOD(year,100).GE.1)).OR.(MOD(year,400).LT.1)) THEN
    IF (.NOT.leap) leapy=.true.  ! If leap=F, we remove leap years
ENDIF ! IF leap year and we account for it

IF(dbg) print*,"Keep leap day, if present?", leap
IF(dbg) print*,"Leapyear?(leapy)", leapy
IF (.NOT.exfile) then ! IF file is not existing
    print*,"File missing, creating a dummy one for",year*100+mn
    missval=attdat%fillval
    leapy=.false.
    IF (((MOD(year,4).lt.1).AND.(MOD(year,100).GE.1)).OR.(MOD(year,400).LT.1)) THEN
        IF (leap) leapy=.true.  ! If leap=T, we keep the leap years (necessary, as we need to extend from 28 -> 29 days instead crop form 29 to 28 (if exfile)
    ENDIF ! IF leap year and we account for it

    SELECT CASE (mn)
        CASE (1,3,5,7,8,10,12)
            days=31
        CASE (4,6,9,11)
            days=30
        CASE (2)
            days=28
            IF (leapy) days=29
        CASE (-1) ! Wild card to indicate that full year is loaded
            days=365
            IF (leapy) days=366
    END SELECT

    ntimes=days*4
    IF (intype=="ccc") THEN ; nx=192; ny=96; nlev=1 ; ENDIF ! IF intype=="ccc"
    !IF (intype=="20cr".OR.intype=="20CR") THEN ; nx=180; ny=91; nlev=1 ; ENDIF ! IF intype="20cr"
    ALLOCATE(arr(nx,ny,nlev,ntimes),times(ntimes),levels(nlev),mems(nlev)) !Includes dummy allocations
    IF (dbg) WRITE(*,'(A24,I3,A4,I3,A6,I3,A8,I4)') "Create dummy matrix. nx=",nx," ny=",ny," nlev=",nlev,"ntimes=",ntimes
    arr=missval
    times=0
    CALL timeaxis(times,year,mn,ntimes,leapy,leap)
    levelid=1
    CALL dummyaxis(lons,lats,nx,ny,intype,attdat,attlon,attlat,atttime,invarnames,missval)
    IF (dbg) print*,"time",ntimes,times(1),times(UBOUND(times))
    IF (dbg) print*,"lat",nx,lons(1),lons(UBOUND(lons))
    IF (dbg) print*,"lon",ny,lats(1),lats(UBOUND(lats))
ELSE IF (leapy) THEN ! IF file exists and has a leap year (which has been declared NOT to be removed!)
    IF (mn.EQ.2) THEN
        IF(dbg) WRITE(*,'(A28,I4)') 'Remove leap year from year: ', year
        days=28
        ntimes=days*4
        arr=arr(:,:,:,1:ntimes)
        ALLOCATE(times(ntimes))
        CALL timeaxis(time,year,mn,ntimes,leapy,leap)
        times=time(1:112) ! Ugly hack!
    ELSE IF (mn.EQ.-1) THEN
        !IF(leap.AND.() ) 
        IF(dbg) WRITE(*,'(A28,I4)') 'Remove leap year from year: ', year
        days=365
        ntimes=days*4
        arr(:,:,:,237:(ntime-4))=arr(:,:,:,241:ntime) !remove leapday
        arr(:,:,:,(ntime-3):ntime)=attdat%fillval
        ALLOCATE(times(ntimes))
        CALL timeaxis(time,year,mn,ntimes,leapy,leap)
        times=time(1:1460) ! Ugly hack!
    ELSE
        CALL timeaxis(time,year,mn,ntime,leapy,leap)
        ntimes=ntime
        ALLOCATE(times(ntimes))
        times=time(1:ntimes)
    ENDIF ! if month.eq.2   
ELSE ! If file exists and is not declared a leap year (or leap year is removed)
    ntimes=ntime
    ALLOCATE(times(ntimes))
    CALL timeaxis(time,year,mn,ntimes,leapy,leap)
    times=time(1:ntimes)
ENDIF ! Flicke missing files and Schaltjahre
atttime%units="day as %Y%m%d.%f"
IF(dbg) PRINT*,"intype=",intype,"exfile=",exfile,"yyyy=",year,"mn=",mn

! If ccc we do a full check of data
IF (intype=="ccc" .AND. exfile) THEN
    !print*,"==> Full check of data!"
    checkup=.true.
ENDIF ! IF intype and exfile

IF (dbg) WRITE(*,'(A28,L2,A14,L2,A6,I3)') "Start checkup: CCC checkup? ", checkup, " File exists? ", exfile, " nlev=",nlev

nmem=size(mems); IF(nmem<1) nmem=1
IF (dbg) print*,"nlev",nlev,"ens?",ens,nmem,size(mems)
DO lvl=1,nmem ! Loop over all members (if present, don't do that for levels!

   IF (nlev.GT.1 .AND. .NOT.ens) THEN
      IF (dbg) print*, "Select level: ", inlev
      IF (dbg) print*, "Avail levels: ", levels
    
      ! We need to select the correct level
      IF (checkup.AND.invarnames=="geopoth") THEN
         IF (.NOT.allocated(vert2)) ALLOCATE(vert2(nlev,ntimes))
         validtmstps=ntimes  ! Initially assume that all timesteps are valid
         take=1              ! Set take to one, the first time step, will be changed later, if the first timestep is erroneous
         DO ii= 1,nlev
            DO tt= 1,ntimes
               mmax = MAXVAL(arr(:,1:3,ii,tt))
               mmin = MINVAL(arr(:,1:3,ii,tt))
               mean = SUM(arr(:,1:3,ii,tt))/(nx*3)
               IF (dbg.AND.tt==1) WRITE(*,'(A4,I2,4(A8,F13.3))') "lev=",ii,"mean=",mean," mmin=",&
                                                                &mmin," mmax=",mmax," offset=",attdat%offset
                
               if ( (mmax-mean.LT.0.001).OR.(mean-mmin.LT.0.001).OR.(ABS(mmin-mmax).LT.0.001) ) THEN
                  times(tt)=year*40000+mn*400+tt+3 ! Most often the time axis is broken too, so we fix this because later on we need a correct timestamp for selyear
                  times(tt)=times(tt)/4
                  arr(:,:,ii,tt)=missval  
                  IF (ii.eq.1) WRITE(*,'(A29,F13.2,A8,3(F13.3))') 'Data missing for time (geop1)',&
                                                                 & times(tt)," Stats: ",mmin,mean,mmax
                  vert2(ii,tt)=INT(missval)

                  IF (ii.eq.1) validtmstps=validtmstps-1
                  IF ((ii.eq.1).and.(tt.eq.take).and.(take.ne.ntime)) take=take+1  ! If missing time steps are present, we have to take a valid one for overwritting the vertical axis
               endif ! if data missing north pole

               mmax = MAXVAL(arr(:,(ny-2):ny,ii,tt))
               mmin = MINVAL(arr(:,(ny-2):ny,ii,tt))
               mean = SUM(arr(:,(ny-2):ny,ii,tt))/(nx*3)
               !print*,"Timestep",tt,mmax,mean,mmin
               if ( (ABS(mmin-mmax).LT.0.001).AND.(vert2(ii,tt)/=INT(missval) ) ) THEN
                  times(tt)=year*40000+mn*400+tt+3 
                  times(tt)=times(tt)/4
                  IF (ii.eq.1) WRITE(*,'(A29,F13.2,A8,3(F13.3))') 'Data missing for time (geop2)',&
                                                                 & times(tt),"Stats:",mmin,mean,mmax
                  vert2(ii,tt)=INT(missval)
                  arr(:,:,ii,tt)=missval 
                  IF (ii.eq.1) validtmstps=validtmstps-1
                  IF ((ii.eq.1).and.(tt.eq.take).and.(take.ne.ntime)) take=take+1 
                 
               endif ! if data missing south pole

               if ( vert2(ii,tt)/=INT(missval) ) THEN
                  mean=SUM(arr(:,:,ii,tt))/(nx*ny)
                  SELECT CASE (INT(mean))
                     CASE (0:300)
                        vert2(ii,tt)=100000
                     CASE (1000:1500)
                        vert2(ii,tt)=85000
                     CASE (5000:6000)
                        vert2(ii,tt)=50000
                     CASE (6500:7500)
                        vert2(ii,tt)=40000
                     CASE (8500:9499)
                        vert2(ii,tt)=30000
                     CASE (9500:10999)
                        vert2(ii,tt)=25000
                     CASE (11000:12499)
                        vert2(ii,tt)=20000
                     CASE (12500:14000)
                        vert2(ii,tt)=15000
                     CASE (14500:16500)
                        vert2(ii,tt)=10000
                     CASE DEFAULT
                        print*,'Strange data, please check:', times(tt), levels(ii),mean
                  END SELECT ! Select mean
               endif ! IF Data here and else CASE level
            ENDDO !tt
         ENDDO !ii
        
         ! Check whether axis is inverted
         pass=0
         !print*,levels,inlev,vert2(:,1)
         DO tt=1,ntimes
            IF(ALL(INT(vert2(:,tt))==INT(levels))) pass = pass + 1  ! Number of time steps with valid levels
            IF (pass==1) THEN
               DO ii=1,nlev
                  IF (ABS(levels(ii)-inlev)<=0.01) levelid=ii
               ENDDO ! DO ii
            ENDIF ! pass==1
         ENDDO ! DO tt

         IF (dbg) print*,"If pass is >0, no problem:",pass,"valid timesteps=",validtmstps,"levelid=",levelid
         IF((pass==0).AND.(validtmstps>0)) THEN
            WRITE(*,'(A32,I4,A7,I2)') "Vertical axis corrected in year ", year," month ",mn
            IF(dbg) print*,"Take,levels,vert2(:,take)",take!,levels,vert2(:,take)
            levels=vert2(:,take)
            DO ii=1,nlev
               IF(ABS(levels(ii)-inlev)<=0.01) levelid=ii
            ENDDO ! DO ii
         ENDIF !pass=0
         IF(levelid==-1) THEN
            levelid=1
            WRITE(*,*) "Level selection failed! Taking levelid 1"
         ENDIF ! IF levelid==-1

      ELSEIF (checkup.AND.invarnames=="slp") THEN
         levelid=1
         CALL timeaxis(time,year,mn,ntimes,leapy,leap)
         times=time(1:ntimes)

         ! Check whether data is missing
         DO tt=1,ntimes
            mmax = MAXVAL(arr(:,1:3,levelid,tt))
            mmin = MINVAL(arr(:,1:3,levelid,tt))
            mean = SUM(arr(:,1:3,levelid,tt))/(nx*3) ! Only first three ny
            IF (dbg.AND.tt==1) WRITE(*,'(A4,I2,4(A8,F13.3))') "lev=",ii,"mean=",mean," mmin=",&
                                                             &mmin," mmax=",mmax," offset=",attdat%offset

            IF ( (ABS((mmax-mean)/mean).LT.0.0001).OR.(ABS((mean-mmin)/mmin).LT.0.0001)&
               &.OR.(ABS((mmin-mmax)/mmax).LT.0.0001) ) THEN
               WRITE(*,'(A32,F13.2,A8,3(F13.3))') 'Data missing for timestep (slp1)',&
                                                   &times(tt)," Stats: ",mmin,mean,mmax
               arr(:,:,levelid,tt)=missval
            ENDIF ! IF all data 
            IF (ABS(arr(1,1,levelid,tt)-missval).GT.0.0001) THEN
               mmax = MAXVAL(arr(:,(ny-2):ny,levelid,tt))
               mmin = MINVAL(arr(:,(ny-2):ny,levelid,tt))
               mean = SUM(arr(:,(ny-2):ny,levelid,tt))/(nx*3)
               IF (dbg.AND.tt==1) WRITE(*,'(A4,I2,4(A8,F13.3))') "lev=",ii,"mean=",mean," mmin=",&
                                                                 &mmin," mmax=",mmax," offset=",attdat%offset

               IF ( (ABS((mmax-mean)/mean).LT.0.0001).OR.(ABS((mean-mmin)/mmin).LT.0.0001)&
                  &.OR.(ABS((mmin-mmax)/mmax).LT.0.0001) ) THEN
                  WRITE(*,'(A32,F13.2,A8,3(F13.3))') 'Data missing for timestep(slp2)',&
                                                                &times(tt)," Stats: ",mmin,mean,mmax
                  arr(:,:,levelid,tt)=missval
               ENDIF ! IF all data 
            ENDIF ! IF timestep not already set to missval
         ENDDO ! tt in ntimes

      ELSE ! IF checkup.AND.invarnames
         DO ii= 1,nlev
            IF (ABS(levels(ii)-inlev)<=0.01) levelid=ii  ! As we compare an Integer with a REAL, with expect slight differences
         ENDDO
         IF (levelid==-1.AND.svarndims==3) THEN
            levelid=1
         ELSE IF (levelid==-1) THEN
            WRITE(*,'(A27,I5,A18)') "Selected level not found! (",inlev,"). Select level 1."
            levelid=1
         ENDIF ! levelid
         IF(dbg) WRITE(*,'(A12,I6,A41,I2)') "Model level ",inlev," selected. Corresponding to level number ",levelid,"."
         CALL timeaxis(time,year,mn,ntimes,leapy,leap)
         times=time(1:ntimes)
      ENDIF ! if checkup

   ELSEIF (ens) THEN
      levelid=lvl
      outfile=outfileraw
      CALL outfilemodder(outfile,lvl)
   ELSEIF (nlev.LE.1 .AND. .NOT.ens) THEN
      levelid=1
   ENDIF ! IF nlev>1

   IF (dbg) PRINT*, "=>Select and repare axis done"
   IF (dbg) PRINT*, "=>Selected level number:",levelid
   IF (dbg.AND.svarndims==4) WRITE(*,'(A22,I3,A5,I3,A6,I2,A5,I4)') "In array has dims: nx=",size(arr,1),",&
							&ny=",size(arr,2),",nlev=",size(arr,3),"nt=",size(arr,4)
   IF (dbg.AND.svarndims==3) WRITE(*,'(A22,I3,A5,I3,A5,I4)') "In array has dims: nx=",size(arr,1),",ny=",&
                                &size(arr,2),",nt=",size(arr,4)
   IF (dbg) WRITE(*,'(A27,I3,A5,I3,A5,I4)') "Allocate output arrays. nx=",nx," ,ny=",ny," nt=",ntimes
   IF (dbg) PRINT*,"ntimes",size(times)
   IF (.NOT.allocated(arr3out)) allocate(arr3out(nx,ny,ntimes))
   IF (.NOT.allocated(arr3out2)) allocate(arr3out2(nx,ny,ntimes))
   IF (.NOT.allocated(arr3temp)) allocate(arr3temp(nx,ny,ntimes))
   !print*,size(arr,1),size(arr,2),size(arr,3),size(arr,4)
   !print*,size(arr3out,1),size(arr3out,2),size(arr3out,3)
   arr3out(:,:,:)=arr(:,:,levelid,1:ntimes)

   ! corrgrid
   ! Check whether NA values are abundant
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   IF( checkup.AND.(invarnames=="geopoth".OR.invarnames=="slp") ) THEN
      CALL corrgrid(arr3out,missval,invarnames)
   ELSEIF ( checkup.AND.(invarnames/="geopoth".OR.invarnames/="slp") )THEN
        SELECT CASE(nrows)
            CASE(-1)
                IF(dbg) PRINT*,"Checkup: Nothing to be done"
            CASE(1:)
                CALL corrgridnongph(arr3out,corrections,missval)
        END SELECT
   ENDIF ! corrgird
   !dbg=.TRUE.
   IF(dbg) print*,"dim arrin:",size(arr),size(arr,1),size(arr,2),size(arr,3),size(arr,4)
   IF(dbg) print*,"dim arrout",size(arr3out),size(arr3out,1),size(arr3out,2),size(arr3out,3)


   ! Possibly change outvarname
   IF(invarnames.EQ."z") invarnames="Z"

   !!! =============================================================================
   !!! ===================== Output Start===========================================
   IF(dbg) print*,'=====================Output Start================================'
   !!! =============================================================================
   !!! ===================== We only write 3D arrays at the moment... ==============

   !CALL nf90_save(nx,ny,ntime,lons,lats,times,arrout,attlon,attlat,attlev,atttime,attdat,outfile,invarnames)
   istat=0
   ! Create new NetCDF file
   istat = nf90_create(trim(outfile), NF90_CLOBBER, ncido) !NF90_HDF5, ncido)
   IF(istat /= nf90_NoErr) THEN ; WRITE(*,*) 'Create file: ', NF90_STRERROR(istat); ENDIF !CALL ABORT ; ENDIF
   IF (dbg) WRITE(*,*) "Create NetCDF file for variable ", trim(invarnames)," and save it as: ",trim(outfile) 
   ! Define Dimension
   istat = nf90_def_dim(ncido,"lon",nx,LonDimID)     ; IF(istat/=NF90_NoErr) print*,'Create Lon', NF90_STRERROR(istat)
   istat = nf90_def_dim(ncido,"lat",ny,LatDimID)     ; IF(istat/=NF90_NoErr) print*,'Create Lat', NF90_STRERROR(istat)
   !istat = nf90_def_dim(ncido,"lev",1,VerDimID)	 ; IF(istat/=NF90_NoErr) print*,'Create Lev', NF90_STRERROR(istat)
   istat = nf90_def_dim(ncido,"time",nf90_unlimited,TimeDimID); IF(istat/=NF90_NoErr) print*,'Create Time', NF90_STRERROR(istat)

   ! Define Variables
   istat = NF90_DEF_VAR(ncido,"lon",NF90_DOUBLE,(/ LonDimID /),LonVarID)
   IF(istat/=NF90_NoErr) print*,'Define lon', NF90_STRERROR(istat)
   istat = NF90_PUT_ATT(ncido,LonVarID, "standard_name",attlon%standard_name)
   istat = NF90_PUT_ATT(ncido,LonVarID, "long_name",attlon%long_name)
   istat = NF90_PUT_ATT(ncido,LonVarID, "units",attlon%units)


   istat = NF90_DEF_VAR(ncido,"lat",NF90_DOUBLE,(/ LatDimID /),LatVarID)
   IF(istat/=NF90_NoErr) print*,'Define lat', NF90_STRERROR(istat)
   istat = NF90_PUT_ATT(ncido,LatVarID, "standard_name",attlat%standard_name)
   istat = NF90_PUT_ATT(ncido,LatVarID, "long_name",attlat%long_name)
   istat = NF90_PUT_ATT(ncido,LatVarID, "units",attlat%units)

   !istat = NF90_DEF_VAR(ncido,"lev",NF90_DOUBLE,(/ VerDimID /),VerVarID)
   !IF(istat/=NF90_NoErr) print*,'Define lev', NF90_STRERROR(istat)
   !istat = NF90_PUT_ATT(ncido,VerVarID, "standard_name",attlev%standard_name)
   !istat = NF90_PUT_ATT(ncido,VerVarID, "long_name",attlev%long_name)
   !istat = NF90_PUT_ATT(ncido,VerVarID, "units",attlev%units)
   !istat = NF90_PUT_ATT(ncido,VerVarID, "positive",attlev%positive)

   istat = NF90_DEF_VAR(ncido,"time",NF90_DOUBLE,(/ TimeDimID /),timeVarID)
   IF(istat/=NF90_NoErr) print*,'Define time ', NF90_STRERROR(istat)
   istat = NF90_PUT_ATT(ncido,TimeVarID, "standard_name",atttime%standard_name)
   istat = NF90_PUT_ATT(ncido,TimeVarID, "long_name",atttime%long_name)
   istat = NF90_PUT_ATT(ncido,TimeVarID, "units",atttime%units)
   istat = NF90_PUT_ATT(ncido,TimeVarID, "calendar",atttime%calendar)

   !istat = NF90_DEF_VAR(ncido,invarnames,NF90_DOUBLE,(/LonDimID,LatDimID,VerDimID,TimeDimID /),DatVarID)
   istat = NF90_DEF_VAR(ncido,invarnames,NF90_DOUBLE,(/LonDimID,LatDimID,TimeDimID /), DatVarID)
   IF(istat/=NF90_NoErr) print*,'Define time', NF90_STRERROR(istat)
   istat = NF90_PUT_ATT(ncido,DatVarID, "standard_name",attdat%standard_name)
   istat = NF90_PUT_ATT(ncido,DatVarID, "long_name",attdat%long_name)
   istat = NF90_PUT_ATT(ncido,DatVarID, "units",attdat%units)
   istat = NF90_PUT_ATT(ncido,DatVarID, "missing_value",attdat%fillval)

   istat = NF90_PUT_ATT(ncido,NF90_GLOBAL, 'title', 'NetCDF checker ')
   istat = NF90_PUT_ATT(ncido,NF90_GLOBAL, 'institution','University of Bern: Group Broennimann')

   istat=NF90_ENDDEF(ncido)
   IF(istat/=NF90_NoErr) print*,'Error while closing definition mode: ', NF90_STRERROR(istat)

   istat = NF90_PUT_VAR(ncido,LonVarID,lons)
   IF(istat/=NF90_NoErr) print*,'Put lon ', NF90_STRERROR(istat)
   istat = NF90_PUT_VAR(ncido,LatVarID,lats)
   IF(istat/=NF90_NoErr) print*,'Put lat ', NF90_STRERROR(istat)
   !istat = NF90_PUT_VAR(ncido,VerVarID,levels)
   !IF(istat/=NF90_NoErr) print*,'Put lev ', NF90_STRERROR(istat)
   istat = NF90_PUT_VAR(ncido,TimeVarID,times)
   IF(istat/=NF90_NoErr) print*,'Put time ', NF90_STRERROR(istat)
   !PRINT*,arr3out(:,:,1)
   istat = NF90_PUT_VAR(ncido,DatVarID,arr3out)
   IF(istat/=NF90_NoErr) print*,'Put data: ', NF90_STRERROR(istat)

   IF(dbg) print*,"NetCDF creation: ", NF90_STRERROR(istat)

   ! Finish NetCDF
   istat =NF90_CLOSE(ncido)
   IF(dbg) print*,"LVL=",lvl,"Nmem=",nmem
ENDDO ! DO lvl (loop of we have several members in same file

! Deallocate variables
DEALLOCATE(lons,lats,levels,times,mems,arr,arr3out,arr3out2,arr3temp)
IF(ALLOCATED(time))        DEALLOCATE(time)
IF(ALLOCATED(arrin1))      DEALLOCATE(arrin1)
IF(ALLOCATED(arrin2))      DEALLOCATE(arrin2)
IF(ALLOCATED(varndims))    DEALLOCATE(varndims)
IF(ALLOCATED(dimnames))    DEALLOCATE(dimnames)
IF(ALLOCATED(varnames))    DEALLOCATE(varnames)
IF(ALLOCATED(vert2))       DEALLOCATE(vert2)
IF(ALLOCATED(corrections)) DEALLOCATE(corrections)
IF(ALLOCATED(lendims))     DEALLOCATE(lendims)

! ================================================================================
!=================================================================================

! ================================================================================
!=================================================================================
CONTAINS

SUBROUTINE getdata(ncid,DatVarID,nx,ny,nlev,ntime,svarndims,dbg,packed,convgp,convaprx,arr)
    USE netcdf
    IMPLICIT NONE
    INTEGER,INTENT(IN)                                  :: ncid,DatVarID,nx,ny,nlev,ntime,svarndims
    LOGICAL,INTENT(IN)                                  :: dbg
    LOGICAL,INTENT(INOUT)                               :: packed,convgp,convaprx
    REAL, DIMENSION(:,:,:,:),ALLOCATABLE,INTENT(OUT)    :: arr
    REAL, DIMENSION(:,:,:,:),ALLOCATABLE                :: arrscale
    REAL, DIMENSION(:,:,:),ALLOCATABLE                  :: arrtest
    ! Read data
IF(dbg) print*,"@GetData:lendims:",lendims,svarndims
IF(dbg) print*,"@GetData:nx,ny,nlev,ntime:",nx,ny,nlev,ntime

! Notes: istat changed to arr instead of arrscale --> arr = arrscale * ...
IF (svarndims==3) THEN
    IF (.NOT. ALLOCATED(arr)) ALLOCATE(arr(nx,ny,1,ntime))
    IF (.NOT. ALLOCATED(arrscale)) ALLOCATE(arrscale(nx,ny,1,ntime))
    IF (.NOT. ALLOCATED(arrtest)) ALLOCATE(arrtest(nx,ny,ntime))

    istat = NF90_GET_VAR(ncid,DatVarID,arrtest)
    IF (istat/=NF90_NoErr) PRINT*,"Get Data:", NF90_STRERROR(istat)
    IF(dbg .AND. istat==NF90_NoErr) PRINT*,'Allocation and GetVar successful, got ',nx*ny*ntime/128/1024, 'MB data'
    arr(:,:,1,:)=arrtest
ELSE IF (svarndims==4) THEN
    IF (.NOT. ALLOCATED(arr)) ALLOCATE(arr(nx,ny,nlev,ntime))
    IF (.NOT. ALLOCATED(arrscale)) ALLOCATE(arrscale(nx,ny,nlev,ntime))
    istat = NF90_GET_VAR(ncid,DatVarID,arr)
    IF (istat/=NF90_NoErr) PRINT*,"Get Data:", NF90_STRERROR(istat)
    IF(dbg .AND. istat==NF90_NoErr) PRINT*,'Allocation and GetVar successful, got ',nx*ny*nlev*ntime/128/1024, 'MB data'
ELSE
    WRITE(*,*) "Unknown number of dimensions, no allocation possible!"
    CALL ABORT
ENDIF ! IF svarndims

! Unpack and convert data
IF ( Compare_Float(attdat%scalef,-1.) .AND. Compare_Float(attdat%offset,-1.) ) packed=.false.
IF (dbg) PRINT*,"Data packed? ", packed
IF (packed) THEN
    IF (dbg) print*,"Unpacking NetCDF"
    IF (convgp) THEN
        IF (dbg) print*,"Conversion gp -> gph"
        arr= ( arr * attdat%scalef + attdat%offset ) / 9.80655
        attdat%units="m"
    ELSEIF (convaprx) THEN
        IF (dbg) print*,"Conversion kgm**2/s --> mm/(6h)"
        arr= ( arr * attdat%scalef + attdat%offset ) * 3600 * 6
        attdat%units="mm/(6h)"
        attdat%long_name="Total precipitation"
        attdat%standard_name="Total precipitation"
    ELSE
        arr= arr * attdat%scalef + attdat%offset
    ENDIF ! IF conversions

ELSE
    IF (convgp) THEN
        IF (dbg) print*,"Conversion gp -> gph"
        arr = arr / 9.80655
        attdat%units="m"
    ELSE IF (convaprx) THEN
        IF (dbg) print*,"Conversion kgm**2/s --> mm/(6h)"
        arr = arr * 3600 * 6
        attdat%units="mm/(6h)"
        attdat%long_name="Total precipitation"
        attdat%standard_name="Total precipitation"
    ELSE
        arr = arr
    ENDIF ! IF conversions
ENDIF ! IF packed
!print*,"Hallo",size(arr,1),size(arr,2),size(arr,3),size(arr,4)

IF(ALLOCATED(arrtest)) DEALLOCATE(arrtest)
IF(ALLOCATED(arrscale)) DEALLOCATE(arrscale)
END SUBROUTINE getdata

! =====================================================================================

SUBROUTINE getattributes(ncid,VarID,attraw)
    USE netcdf
    IMPLICIT NONE
    INTEGER,INTENT(IN)                              :: ncid,VarID
    TYPE (attributes), INTENT(INOUT),OPTIONAL       :: attraw
    
    IF (VarID.lt.0) THEN
       attraw%units="NULL"
    ELSE
       istat = NF90_GET_ATT(ncid,VarID,"units",attraw%units)
       if (istat/=NF90_NoErr) attraw%units="NULL"
       istat = NF90_GET_ATT(ncid,VarID,"long_name",attraw%long_name)
       if (istat/=NF90_NoErr) attraw%long_name="NULL"
       istat = NF90_GET_ATT(ncid,VarID,"standard_name",attraw%standard_name)
       if (istat/=NF90_NoErr) attraw%standard_name="NULL"
       istat = NF90_GET_ATT(ncid,VarID,"positive",attraw%positive)
       if (istat/=NF90_NoErr) attraw%positive="NULL"
       istat = NF90_GET_ATT(ncid,VarID,"calendar",attraw%calendar)
       if (istat/=NF90_NoErr) attraw%calendar="NULL"
       istat = NF90_GET_ATT(ncid,VarID,"_FillValue",attraw%FillVal)
       if (istat/=NF90_NoErr) attraw%FillVal=-1
       istat = NF90_GET_ATT(ncid,VarID,"scale_factor",attraw%scalef)
       if (istat/=NF90_NoErr) attraw%scalef=-1
       istat = NF90_GET_ATT(ncid,VarID,"add_offset",attraw%offset)
       if (istat/=NF90_NoErr) attraw%offset=-1
   ENDIF
END SUBROUTINE getattributes

! =====================================================================================

SUBROUTINE readaxis(ncid,dimn,LonDimID,LatDimID,VerDimID,TimeDimID,times,lons,lats,levels,mems,&
                    &nx,ny,ntime,nlev,attlon,attlat,attlev,atttime,axpacked)
    USE netcdf
    IMPLICIT NONE
    INTEGER,INTENT(IN)                                      :: ncid,dimn
    INTEGER,INTENT(IN)                                      :: LonDimID, LatDimID, VerDimID, TimeDimID
    INTEGER,INTENT(IN)                                      :: nx,ny,ntime,nlev
    INTEGER                                                 :: istat
    REAL(8),DIMENSION(:),ALLOCATABLE,INTENT(INOUT),OPTIONAL :: times,lons,lats,levels,mems
    TYPE (attributes), INTENT(INOUT),OPTIONAL               :: attlon,attlat,attlev,atttime
    LOGICAL                                                 :: packed
    LOGICAL,INTENT(INOUT)                                   :: axpacked
    
    packed=.TRUE.
    IF (dimn==LonDimID ) THEN
        allocate(lons(nx))
        istat = NF90_GET_VAR(ncid,dimn,lons)!; print*, NF90_STRERROR(istat)
        IF ( Compare_Float(attlon%scalef,-1.) .AND. Compare_Float(attlon%offset,-1.) ) packed=.FALSE.
        IF (packed) THEN
            lons = lons * attlon%scalef + attlon%offset 
            axpacked=.TRUE.
        ENDIF ! IF packed

    ELSE IF (dimn==LatDimID) THEN
        allocate(lats(ny))
        istat = NF90_GET_VAR(ncid,dimn,lats)
        IF ( Compare_Float(attlat%scalef,-1.) .AND. Compare_Float(attlat%offset,-1.) ) packed=.false.
        IF (packed) THEN
            lats = lats * attlat%scalef + attlat%offset
            axpacked=.TRUE.
        ENDIF ! IF packed

    ELSE IF (dimn==VerDimID) THEN !VerDimID) THEN
        IF (ens) THEN
            ALLOCATE(mems(nlev))
            ALLOCATE(levels(1)) ; levels=-1
            istat = NF90_GET_VAR(ncid,dimn,mems)
        ELSE
            ALLOCATE(levels(nlev))
            ALLOCATE(mems(1)) ; mems=-1
            istat = NF90_GET_VAR(ncid,dimn,levels)

            IF ( Compare_Float(attlev%scalef,-1.) .AND. Compare_Float(attlev%offset,-1.) ) packed=.false.
            IF (packed) THEN
                levels = levels * attlev%scalef + attlev%offset
                axpacked=.TRUE.
            ENDIF
        ENDIF ! If level or member
    ELSE IF (dimn==TimeDimID) THEN !TimeDimID) THEN
        ALLOCATE(times(ntime))
        istat = NF90_GET_VAR(ncid,dimn,times)
        IF ( Compare_Float(atttime%scalef,-1.) .AND. Compare_Float(atttime%offset,-1.) ) packed=.false.
        IF (packed) THEN
            times = times * atttime%scalef + atttime%offset
            axpacked=.TRUE.
        ENDIF

    ENDIF
END SUBROUTINE readaxis

! =====================================================================================

SUBROUTINE dummyaxis(lons,lats,nx,ny,intype,attdat,attlon,attlat,atttime,invarnames,missval)
    IMPLICIT NONE
    ! Create dummy axis for non-existant files. This consists of lon,lat and the attributes
    CHARACTER(32),INTENT(IN)                        :: intype,invarnames
    INTEGER,INTENT(IN)                              :: nx, ny
    REAL(8),DIMENSION(:),ALLOCATABLE,INTENT(INOUT)  :: lons,lats
    REAL,INTENT(IN)                                 :: missval
    TYPE (ATTRIBUTES), INTENT(INOUT)                :: attdat,attlon,attlat,atttime

    INTEGER                                         :: ii
    REAL                                            :: dx,dy
  
    IF (ALLOCATED(lons) .EQV. .FALSE. ) ALLOCATE(lons(nx))
    IF (ALLOCATED(lats) .EQV. .FALSE. ) ALLOCATE(lats(ny))
    IF (intype=="ccc") THEN

        IF (ny/=96) PRINT*,"WARNING! Number of latitudes differs from 96 for ccc:", ny

        lats=(/88.5721685140073, 86.7225309546681, 84.8619702920424, &
     82.9989416428375, 81.1349768376774, 79.2705590348597, 77.4058880820788, &
     75.541061452879, 73.6761323132091, 71.8111321142745, 69.9460806469834, &
     68.0809909856513, 66.2158721139987, 64.3507304088721, 62.4855705220364, &
     60.6203959268265, 58.7552092693799, 56.8900126013571, 55.0248075383117, &
     53.1595953700197, 51.2943771389511, 49.429153697123, 47.5639257479787, &
     45.6986938777018, 43.8334585789513, 41.9682202690754, 40.1029793042494, &
     38.2377359905648, 36.3724905928122, 34.507243341501, 32.6419944385177, &
     30.7767440617232, 28.9114923687178, 27.0462394999448, 25.1809855812706, &
     23.3157307261409, 21.4504750373982, 19.5852186088223, 17.7199615264474, &
     15.8547038696949, 13.9894457123567, 12.1241871234558, 10.2589281680064, &
     8.39366890769239, 6.52840940147998, 4.66314970617789, 2.79788987695673, &
     0.932629967838004, -0.932629967838004, -2.79788987695673, &
     -4.66314970617789, -6.52840940147998, -8.39366890769239, &
     -10.2589281680064, -12.1241871234558, -13.9894457123567, &
     -15.8547038696949, -17.7199615264474, -19.5852186088223, &
     -21.4504750373982, -23.3157307261409, -25.1809855812706, &
     -27.0462394999448, -28.9114923687178, -30.7767440617232, &
     -32.6419944385177, -34.507243341501, -36.3724905928122, &
     -38.2377359905648, -40.1029793042494, -41.9682202690754, &
     -43.8334585789513, -45.6986938777018, -47.5639257479787, &
     -49.429153697123, -51.2943771389511, -53.1595953700197, &
     -55.0248075383117, -56.8900126013571, -58.7552092693799, &
     -60.6203959268265, -62.4855705220364, -64.3507304088721, &
     -66.2158721139987, -68.0809909856513, -69.9460806469834, &
     -71.8111321142745, -73.6761323132091, -75.541061452879, &
     -77.4058880820788, -79.2705590348597, -81.1349768376774,&
     -82.9989416428375, -84.8619702920424, -86.7225309546681, -88.5721685140073 /)

        lons(1)=0; dy=1.875 ; dx=-1
        DO ii=2,nx
            lons(ii)=lons(ii-1)+dy
        END DO ! DO lat
        attlon%long_name="longitude"    ; attlon%standard_name="longitude"
        attlat%long_name="latitude"     ; attlat%standard_name="latitude"
        atttime%long_name="time"        ; atttime%standard_name="time"
        attlon%units="degrees_east"     ; attlat%units="degrees_north"
        atttime%units="day as %Y%m%d.%f"; atttime%calendar="standard"
        attdat%fillval=missval
        
        SELECT CASE (invarnames)
            CASE ("geopoth")
                attdat%long_name="geopotential height"
                attdat%standard_name="geopotential heigth"
                attdat%units="m"
            CASE ("aprl","aprs","aprc","apr")
                attdat%long_name="precipitation"
                attdat%standard_name="precipitation"
                attdat%units="mm/timestep"
            CASE ("slp")
                attdat%long_name="mean sea level pressure"
                attdat%standard_name="mean sea level pressure"
                attdat%units="Pa"
            CASE DEFAULT
                PRINT*,"Create dummy file: Invarnames: Variable not found"
                attdat%long_name="NULL"
                attdat%standard_name="NULL"
                attdat%units="NULL"
        END SELECT ! invarnames
    ENDIF ! IF intype=="CCC"
    IF (dbg) PRINT*,"Dummy axes successfully created"
END SUBROUTINE dummyaxis

! =====================================================================================

SUBROUTINE invar(invarnames,invar1,invar2,inlev1,inlev2,mode)
    IMPLICIT NONE

    CHARACTER(32),INTENT(IN)        :: invarnames
    CHARACTER(32),INTENT(OUT)       :: invar1, invar2, mode
    INTEGER, INTENT(OUT)            :: inlev1, inlev2
    CHARACTER(32)                   :: inl1, inl2
    INTEGER                         :: strlen,ii
    CHARACTER(32)                   :: var1, var2
    
    IF (dbg) PRINT*,"====================== Read in variables ====================="
    strlen=len(trim(invarnames))
    IF (dbg) PRINT*,strlen,trim(invarnames)

    mode="onevar"
    ii=SCAN(invarnames,"+")
    IF (ii>0) THEN
        mode="plus"
        var1=invarnames(1:ii-1); var2=invarnames(ii+1:strlen)
    ELSE
        ii=SCAN(invarnames,"-")
        IF (ii>0) mode="minus"
        var1=invarnames(1:ii-1); var2=invarnames(ii+1:strlen)
    ENDIF ! IF ii>0
    IF (dbg.AND.mode/="onevar") WRITE(*,'(A21,A16,A16)') "Two variables found: ",var1,var2

    IF (mode=="onevar") THEN
        IF (dbg) WRITE(*,'(A21,A16)') "One variable found: ",invarnames
        ii=SCAN(invarnames,",")
        IF(ii>0) THEN
            IF (dbg) PRINT*,"Height level given"
            invar1=invarnames(1:ii-1); invar2="NULL"
            inl1=invarnames(ii+1:strlen); inl2="-1"
        ELSE
            invar1=trim(invarnames); invar2="NULL"
            inl1="-1" ; inl2="-1"
        ENDIF ! IF ii>0
    ELSE
        IF (dbg) PRINT*,"Two variables found"
        ii=SCAN(var1,",")
        IF (ii>0) THEN
            print*,"height level given for var1"
            invar1=var1(1:ii-1); inl1=var1(ii+1:len(trim(var1)))
        ELSE
            invar1=var1; inl1="-1"
        ENDIF ! IF height level found for var1
        ii=SCAN(var2,",")
        IF (ii>0) THEN
            print*,"height level given for var2"
            invar2=var2(1:ii-1); inl2=var2(ii+1:len(trim(var2)))
        ELSE
            invar2=var2; inl2="-1"
        ENDIF ! IF height level found for var1
    
    ENDIF
    READ (inl1,'(I6)') inlev1
    READ (inl2,'(I6)') inlev2
END SUBROUTINE invar

! =====================================================================================

SUBROUTINE timeaxis(times,year,mn,ntimes,leapy,leap)
    IMPLICIT NONE

    REAL(8),DIMENSION(:),INTENT(INOUT)              :: times
    INTEGER, INTENT(IN)                             :: year,mn,ntimes
    LOGICAL, INTENT(INOUT)                             :: leapy,leap
    INTEGER                                         :: jj


    leapy=.false.
    IF (leap) THEN
        IF (((MOD(year,4).lt.1).AND.(MOD(year,100).GE.1)).OR.(MOD(year,400).LT.1)) THEN
            leapy=.true.
        ENDIF
    ENDIF ! leap
    DO jj=1,ntimes
        IF ( (INT(times(jj))-year*10000.gt.9999).or.(INT(times(jj))-year*10000.lt.0) ) then
            IF (dbg.AND.jj==1) print*,' Time axis corrected for ', year
        ENDIF ! IF time axis wrong
        
        IF (mn > 0 ) THEN
            times(jj)=year*40000+mn*400+jj+3
            times(jj)=times(jj)/4  ! I don't know, why it's not working directly.. but it works (precision?)
        ELSE IF (mn == -1 .AND. leapy) THEN
            select case (jj)
               case(:124) ! January
                times(jj)=year*40000+400+jj+3
               case(125:240) ! February
                times(jj)=year*40000+800+jj-121
               case(241:364) ! March
                times(jj)=year*40000+1200+jj-237
               case(365:484) ! April
                times(jj)=year*40000+1600+jj-361
               case(485:608) ! May
                times(jj)=year*40000+2000+jj-481
               case(609:728) ! June
                times(jj)=year*40000+2400+jj-605
               case(729:852) ! July
                times(jj)=year*40000+2800+jj-725
               case(853:976) ! August
                times(jj)=year*40000+3200+jj-849
               case(977:1096) ! September
                times(jj)=year*40000+3600+jj-973
               case(1097:1220) ! October
                times(jj)=year*40000+4000+jj-1093
               case(1221:1340) ! November
                times(jj)=year*40000+4400+jj-1217
               case(1341:1464) ! December
                times(jj)=year*40000+4800+jj-1337
            end select
            times(jj)=times(jj)/4
        ELSE IF (mn == -1 .AND. .NOT.leapy) THEN
            select case (jj)
               case(:124) ! January
                times(jj)=year*40000+400+jj+3
               case(125:236) ! February
                times(jj)=year*40000+800+jj-121
               case(237:360) ! March
                times(jj)=year*40000+1200+jj-233
               case(361:480) ! April
                times(jj)=year*40000+1600+jj-357
               case(481:604) ! May
                times(jj)=year*40000+2000+jj-477
               case(605:724) ! June 
                times(jj)=year*40000+2400+jj-601
               case(725:848) ! July
                times(jj)=year*40000+2800+jj-721
               case(849:972) ! August
                times(jj)=year*40000+3200+jj-845
               case(973:1092) ! September
                times(jj)=year*40000+3600+jj-969
               case(1093:1216) ! October
                times(jj)=year*40000+4000+jj-1089
               case(1217:1336) ! November
                times(jj)=year*40000+4400+jj-1213
               case(1337:1460) ! December
                times(jj)=year*40000+4800+jj-1333
            end select
            times(jj)=times(jj)/4
        ELSE
            WRITE(*,'(A26,I3)') "Invalid monthly input. mn="
        ENDIF ! Date wrong
    ENDDO ! do jj=1,nTime

END SUBROUTINE timeaxis

! =====================================================================================

SUBROUTINE corrgrid(arr3out,missval,invarnames)
    IMPLICIT NONE

    REAL(4),DIMENSION(:,:,:),INTENT(INOUT)              :: arr3out
    REAL,INTENT(IN)                                     :: missval
    CHARACTER(100),INTENT(IN)                           :: invarnames
    INTEGER                                             :: ntimes
    INTEGER                                             :: nx,ny ! #lon,#lat grid points 
    INTEGER                                             :: tt,jj,kk ! Count variables
    REAL(4),DIMENSION(:,:,:),ALLOCATABLE                :: arr3corr ! Helper array
    REAL(4)                                             :: mmin,mmax ! Calculate min/max of timestep
    REAL(4)                                             :: difference1,difference2 ! Difference between adjacent latitudes
    REAL(4)                                             :: maxdiffnp,maxdiffsp ! Maximum difference at North/South Pole
    INTEGER                                             :: offset1,offset2
    LOGICAL                                             :: shift
    INTEGER                                             :: nywa ! ny without antarctica
    REAL(4),DIMENSION(:),ALLOCATABLE                    :: zonmean, dzonmeandlat
    REAL(4),DIMENSION(:),ALLOCATABLE                    :: zonvar, relzvard
    REAL(4),DIMENSION(:,:),ALLOCATABLE                  :: dzdlat
    INTEGER                                             :: xzonloc,xvarloc
    INTEGER                                             :: xdzonloc
    REAL                                                :: xdzon
    REAL                                                :: zvarthres ! Measures ratio of variance difference between to latitude bands
    LOGICAL,DIMENSION(:,:),ALLOCATABLE                  :: mask

    nx=size(arr3out,1); ny=size(arr3out,2) ; ntimes=size(arr3out,3)
    SELECT CASE (invarnames)
        CASE ("geopoth")
            maxdiffnp=180 ; maxdiffsp=500  ; zvarthres=1.3 ! 
        CASE ("slp")
            maxdiffnp=2000; maxdiffsp=6000 ; zvarthres=1.3 ! 40 hPa at South Pole
    END SELECT ! CASE invarnames
    ! CASE st,u,v,relhum,v10,u10,temp2,..., cannot be adapted because we have realistic values at the bottom/top (shift cannot be detectd...)

    ALLOCATE(arr3corr(nx,ny,ntimes),mask(nx,ny))
    arr3corr=0

    DO tt=1,ntimes
        mask = INT(arr3out(:,:,tt)).EQ.INT(missval)
        mmax = maxval(arr3out(:,:,tt))
        mmin = minval(arr3out(:,:,tt))
        offset1=0
       
        DO kk=1,nx
            IF (dbg) write(*,'(3(F11.1),I3)') arr3out(kk,1:3,tt),kk
            difference1 = abs(arr3out(kk,3,tt)-arr3out(kk,2,tt))
            difference2 = abs(arr3out(kk,2,tt)-arr3out(kk,1,tt))
            IF ( (difference2.GE.(5*difference1)).AND.(difference2.GE.maxdiffnp) ) THEN
                offset1=kk
            ELSE
                EXIT ! exit loop
            ENDIF ! diff1.ge.5*diff2
        ENDDO ! do kk
        IF (offset1.GT.0) THEN
            ! Move gridpoint back to expected position
            arr3corr((nx-offset1+1):nx,1:(ny-1),tt) =arr3out(1:offset1,2:ny,tt)
            arr3corr((nx-offset1+1):nx,ny,tt)       =arr3out(1:offset1,ny,tt)
            arr3corr(1:(nx-offset1),:,tt)           =arr3out((offset1+1):nx,:,tt)
            arr3out(:,:,tt)                         =arr3corr(:,:,tt)
        ENDIF ! if offset.lt.0
    
        offset2=0
        DO kk=nx,1,-1 ! For each longitude starting from highest
            IF (dbg) write(*,'(3(F11.1),I3)') arr3out(kk,(ny-2):ny,tt),kk
            difference1 = abs(arr3out(kk,(ny-2),tt)-arr3out(kk,(ny-1),tt))
            difference2 = abs(arr3out(kk,(ny-1),tt)-arr3out(kk,(ny-0),tt))
            IF ( (difference2.GE.(5*difference1)).AND.(difference2.GE.maxdiffsp) ) THEN ! Difference2 must be larger, because of problems in Antarctica...
                offset2=offset2+1
            ELSE
                EXIT ! exit loop
            ENDIF ! id diff2.ge.5*diff1
        ENDDO ! do kk=nLong,1,-1
    
        IF (offset2.GT.0) THEN
            ! Move gridpoints back to expected position
            arr3corr(1:offset2,2:ny,tt)             =arr3out((nx-offset2+1):nx,1:(ny-1),tt)
            arr3corr(1:offset2,1,tt)                =arr3out((nx-offset2+1):nx,1,tt)
            arr3corr((offset2+1):nx,:,tt)           =arr3out(1:(nx-offset2),:,tt)
            arr3out(:,:,tt)                         =arr3corr(:,:,tt)
        ENDIF ! if offset2
    
        ! Check for unrealistic jumps within the gridded field
        IF(tt.EQ.1) ALLOCATE(zonmean(ny),dzonmeandlat(ny),zonvar(ny),relzvard(ny),dzdlat(nx,ny))
        zonmean=0 ; dzonmeandlat=0 ; zonvar=0 ; relzvard=0 ; shift=.false.
        nywa=INT(0.9*ny) ! Exclude antarctica (~72-90S)from the analysis, as we get strange signal there which obfuscicate the results
        IF ((offset1.GT.0).OR.(offset2.GT.0)) THEN
            DO jj=1,nywa ! do Latitude (without (Ant)arctica), only longitudes which were not shifted before!
                DO kk=(1+offset2),(nx-1-offset1) ! do Longitude (exclude lons which were offsetted, because do not have the jumps
                    dzdlat(kk,jj)=arr3out(kk,jj,tt)-arr3out(kk,jj+1,tt)
                ENDDO !kk
                zonmean(jj)=SUM(dzdlat(:,jj))/float(nx) ! No abs for Variance first
                DO kk=(1+offset2),(nx-1-offset1) ! do Longitude
                    zonvar(jj)=sum(abs(dzdlat(:,jj)-zonmean(jj)))/float(nx)
                ENDDO !KK (zonvar)
                zonmean(jj)=sum(abs(dzdlat(:,jj)))/float(nx) ! Just want the absolute value ?!?!

                IF (jj.GT.1) THEN
                    dzonmeandlat(jj)=abs( zonmean(jj)-zonmean(jj-1))
                    IF(jj.GT.2) relzvard(jj-1)=zonvar(jj-1)/( zonvar(jj-2)+zonvar(jj) ) ! Calculate relative zonal variance 
                    ! The absolute zonal variance varies greatly (lower at eq, higher at midlats, so we use relative values
                ENDIF ! if jj
            ENDDO ! do jj=1,nywa
    
            relzvard(1:3)=0 ! Exclude Arctic ?
            xzonloc = maxloc(dzonmeandlat(:),dim=1) ! maxloc/minloc return a 1D array, but we just want on integer
            xvarloc = maxloc(relzvard(:),dim=1)
            xdzonloc = maxloc( zonmean((xzonloc-1):(xzonloc+1)),dim=1 ) +(xzonloc-2) ! Get the correct index, thus move output by maxzonmean-2
            xdzon = sum(dzonmeandlat(xdzonloc:xdzonloc+1))/2
            IF (dbg) print*,"maxzonloc|maxdzon|maxdzonloc|maxvarloc",xzonloc,xdzon,xdzonloc,xvarloc
 
            IF(maxval(relzvard).GE.zvarthres.AND.(xvarloc.LT.90)) THEN
                if (offset1.GT.0) THEN
                    IF(dbg)print*,"Shifting=true by offset1:",offset1
                    arr3corr(:,:,tt)                          =missval
                    arr3corr(:,1:xvarloc,tt)                  =arr3out(:,1:xvarloc,tt)
                    arr3corr(1:offset1,(xvarloc+2):ny,tt)     =arr3out((nx-offset1+1):nx,(xvarloc+1):(ny-1),tt)
                    arr3corr((offset1+1):nx,(xvarloc+1):ny,tt)=arr3out(1:(nx-offset1),(xvarloc+1):ny,tt)
                    arr3corr(1:offset1,(xvarloc+1),tt)        =arr3out((nx-offset1+1):nx,xvarloc,tt)
                    arr3corr((nx-offset1+1):nx,xvarloc,tt)  =sum(arr3corr((nx-offset1+1):nx,[(xvarloc-1),(xvarloc+1)],tt),dim=2)/2
                    arr3out(:,:,tt)                           =arr3corr(:,:,tt)
                    shift=.true.
                endif ! offset1
        
                if (offset2.GT.0) THEN
                    print*,'Not yet supported, set timestep to missing values'
                    shift=.true.
                    arr3out(:,:,tt)=missval
                endif ! offset2
            ENDIF ! zonvar too large (maxval(relzvard)
        ENDIF ! IF offset1 or offset2 > 0
    
        ! Checking whether shifting helped...
        IF(dbg.AND.shift) print*,'Shift,Offset1,Offset2,',shift,offset1,offset2!offset1=0; offset2=0
        IF ((shift).OR.((offset1.EQ.0).AND.(offset2.EQ.0))) THEN
            IF(dbg.AND.shift) print*,'    Check whether shifting helped...'
            IF(dbg.AND.shift) print*,"nywa",nywa,tt,"offset1|2",offset1,offset2,shift,"nx",nx
       
            DO jj=1,nywa ! do Latitude, again without the critical parts at the poles, only longitudes which were not shifted before
                DO kk=(1+offset2),(nx-1-offset1) ! do Longitude
                    dzdlat(kk,jj)=arr3out(kk,jj,tt)-arr3out(kk,jj+1,tt)
                ENDDO !kk
                zonmean(jj)=sum(dzdlat(:,jj))/float(nx)
                DO kk=(1+offset2),(nx-1-offset1) ! do Longitude
                    zonvar(jj)=sum(abs(dzdlat(:,jj)-zonmean(jj)))/float(nx)
                ENDDO !KK (zonvar)
                zonmean(jj)=sum(abs(dzdlat(:,jj)))/float(nx)

                IF (jj.GT.1) THEN
                    dzonmeandlat(jj)=abs( zonmean(jj)-zonmean(jj-1))
                    IF(jj.GT.2) relzvard(jj-1)=zonvar(jj-1)/( zonvar(jj-2)+zonvar(jj) )
                ENDIF ! if jj.gt.1)
                relzvard(1:3)=0
                xvarloc=maxloc(relzvard,dim=1) ! Exclude Antarctica
            ENDDO ! JJ
        ENDIF ! Shift control

        200 FORMAT(A18,I3,A58,F12.2,F7.3,I3)
        201 FORMAT(A58,F12.2,F7.3,I3)
        202 FORMAT(A58,F12.2,F7.3,I3)
        203 FORMAT(A25,I3,A30,F12.2,F7.3,I3)
        !204 FORMAT(A11,F12.2,F7.3)

        IF ( (maxval(relzvard).GE.zvarthres).AND.(shift) ) THEN !2
            WRITE(*,200)'Grid offsetted by ',offset2-offset1,' to the right and shifted,&
                         &unsuccessful recovery, set to NA', times(tt),maxval(relzvard),xvarloc
            arr3out(:,:,tt)=missval
        ELSEIF ( (maxval(relzvard).LT.zvarthres).AND.(shift) ) THEN !4
            WRITE(*,201)'Grid shifted, successful recovery, data usable!', times(tt), maxval(relzvard), xvarloc
        ELSEIF ( (maxval(relzvard).GE.zvarthres).AND.(.NOT.shift).AND.( (xvarloc.LE.nywa).OR.(xvarloc.GE.INT(0.05*ny)) ) ) THEN !4
            WRITE(*,202)'Unrecoverable shift detected (offset=0), set to NA at       ', times(tt), maxval(relzvard), xvarloc
            arr3out(:,:,tt)=missval
        ELSEIF ( (maxval(relzvard).LT.zvarthres).AND.( (offset1.GT.0).OR.(offset2.GT.0) ) ) THEN !3
            WRITE(*,203)'Grid offset corrected by ',offset1-offset2,' to the left, everything ok at',&
                        &times(tt), maxval(relzvard),xvarloc
        ELSEIF ( (maxval(relzvard).LT.zvarthres).AND.(.NOT.shift).AND.(offset1.EQ.0).AND.(offset2.EQ.0) ) THEN !2
            !WRITE(*,204),'Timestep ok',times(tt),maxval(relzvard)
        ELSEIF ( (maxval(relzvard).GE.zvarthres).AND.(.NOT.shift).AND.((xvarloc.GT.nywa).or.(xvarloc.LT.INT(0.05*ny)) ) ) THEN !2
            !WRITE(*,204),'Timestep ok',times(tt),maxval(relzvard)
        ENDIF ! if shift
    
    ENDDO ! do tt

    DEALLOCATE(arr3corr,zonmean,dzonmeandlat,zonvar,relzvard,dzdlat,mask)

END SUBROUTINE corrgrid

! =====================================================================================

SUBROUTINE createfortfail(year,mn,member,existfort,nn2)
    IMPLICIT NONE

    INTEGER,DIMENSION(:,:),ALLOCATABLE :: indata, outdata
    INTEGER                            :: nn, istat, ii
    INTEGER,INTENT(IN)                 :: year, mn, member
    INTEGER,INTENT(OUT)                :: nn2
    LOGICAL,INTENT(INOUT)              :: existfort

    ! Read fort.fails
    INQUIRE(FILE="datafails.txt", EXIST=existfort)
    IF(dbg) print*,"L1338: existfort: ",existfort
    IF (existfort) THEN
        OPEN(100,file="datafails.txt",status="old",ACTION="read")
        nn=0
        DO 
            READ(100,*,IOSTAT=istat)
            IF ( istat == -1 ) EXIT
            nn=nn+1
        ENDDO ! Get number of lines in file
        ALLOCATE(indata(nn-1,6),outdata(nn-1,6))
        REWIND(100) ! Rewind file to start
        indata=0; outdata=0
        READ(100,*)
        DO ii=2,nn
            READ(100,*,IOSTAT=istat) indata(ii-1,:)
        ENDDO
  
        nn2=0 ! ==nrows
        DO ii=1,nn-1
            IF (indata(ii,1)==member) THEN
                IF (indata(ii,2)==year.OR.indata(ii,2)==0) THEN
                    IF (indata(ii,3)==mn.OR.indata(ii,3)==0) THEN
                        nn2=nn2+1
                        outdata(nn2,:)=(indata(ii,:))
                    ENDIF ! if month matching (OR wildcard)
                ENDIF ! if year matching (OR wildcard)
            ENDIF ! if member matching
        ENDDO ! Do for each entry in datafails.txt
        IF (dbg) print*,"L1366: nn2: ", nn2, "member: ",member,"year:",year,"mn:",mn
        IF (nn2.GT.0) THEN
            OPEN(101,file="fort.fails",ACTION="write")
            DO ii=1,nn2
                WRITE(101,'(3(I8))') outdata(ii,4:6)
            ENDDO ! Do for each column
            CLOSE(101)
            existfort=.TRUE.
        ELSE
            existfort=.FALSE.
        ENDIF ! if nn2.GT.0 (nrows in fort.fails)
        CLOSE(100)
        DEALLOCATE(indata,outdata)
    ENDIF ! If existfort

END SUBROUTINE createfortfail 

! =====================================================================================

SUBROUTINE fortfail(corrections,nrows)
    
    INTEGER,INTENT(OUT)                                     :: nrows
    INTEGER                                                 :: ii               ! Count variable
    INTEGER                                                 :: istat
    INTEGER,DIMENSION(:,:),ALLOCATABLE,INTENT(OUT)          :: corrections
    OPEN(80,FILE="fort.fails",STATUS="old",ACTION="read")
    
    nrows=0
    DO ! Get the number of rows
        READ(80,*, IOSTAT=istat)
        IF ( istat == -1 ) EXIT
        nrows=nrows+1
    ENDDO ! Get nrows
    REWIND(80) ! Rewind file to beginning to read in 

    ALLOCATE(corrections(nrows,3))
    DO ii = 1,nrows
        READ(80,*) corrections(ii,:)
    END DO ! Read file into array
    !944 FORMAT (A10,I3,A32,I4,A7,I2)
    !IF (nrows.GT.0) WRITE(*,944) "There are ", nrows, " corrections to be made in year ", year, " month ", mn

END SUBROUTINE fortfail

! =====================================================================================

SUBROUTINE corrgridnongph(arr3out,corrections,missval)

    REAL(4),DIMENSION(:,:,:),INTENT(INOUT)                  :: arr3out
    INTEGER,DIMENSION(:,:),INTENT(IN)                       :: corrections
    REAL,INTENT(IN)                                         :: missval
    INTEGER                                                 :: nrows,nt,nx,ny
    INTEGER                                                 :: ii
    INTEGER                                                 :: dxa,dxb,reason       ! Selected timesteps to correct data | ID of error
    INTEGER                                                 :: shift                ! Calculated shift out of error number
    REAL(4),DIMENSION(:,:,:),ALLOCATABLE                    :: arr3

    nrows=size(corrections,dim=1)
    nx=size(arr3out,dim=1); ny=size(arr3out,dim=2); nt=size(arr3out,dim=3)
    PRINT*,nx,ny,nt
    ALLOCATE(arr3(nx,ny,nt))
    
    DO ii=1,nrows ! Make all corrections

    IF (corrections(ii,1)==0) then
        dxa=1 ! 0 is the wildcard for the whole month
        dxb=nt
    ELSE
        dxa=MAX(1,corrections(ii,1)-1) ! To be sure that we erase every timestep affected, we also erase the timestep
        dxb=MIN(nt,corrections(ii,2)+1) ! before and after known erroneous timesteps (dxa,dxb cannot be <1,>nt)
    ENDIF ! Correction month is 0 or [1-12]
    reason=corrections(ii,3)

    
    SELECT CASE (reason)
        CASE(1)
            951 FORMAT (A23,I3,A4,I3)
            WRITE(*,951) "Data missing from step ", dxa," to ", dxb
            arr3out(:,:,dxa:dxb)=missval
        CASE(200:299)
            952 FORMAT (A18,I3,A4,I3,A7)
            WRITE(*,952) "Grid shifted from ", dxa, " to ", dxb, ": Erase"
            shift=reason-200
            arr3((nx-shift+1):nx,1:(ny-1),dxa:dxb)      =missval!arr3out(1:shift,2:ny,dxa:dxb)
            arr3((nx-shift+1):nx,ny,dxa:dxb)            =missval!arr3out(1:shift,ny,dxa:dxb)
            arr3(1:(nx-shift),:,dxa:dxb)                =missval!arr3out((shift+1):nx,:,dxa:dxb)
            arr3out(:,:,dxa:dxb)                        =arr3(:,:,dxa:dxb)
            arr3=missval
        CASE(300:399)
            953 FORMAT (A18,I3,A4,I3,A7)
            WRITE(*,953) "Grid shifted from ", dxa, " to ", dxb, ": Erase"
            shift=reason-300
            arr3(1:shift,2:ny,dxa:dxb)                  =missval!arr3out((nx-shift+1):nx,1:(ny-1),dxa:dxb)
            arr3(1:shift,1,dxa:dxb)                     =missval!arr3out((nx-shift+1):nx,1,dxa:dxb)
            arr3((shift+1):nx,:,dxa:dxb)                =missval!arr3out(1:(nx-shift),:,dxa:dxb)
            arr3out(:,:,dxa:dxb)                        =missval!arr3(:,:,dxa:dxb)
        CASE(30000:49999)
            954 FORMAT (A50,I3,A4,I3)
            WRITE(*,954) "Grid shifted and recovery needed - Set to NA from ", dxa," to ", dxb ! Not worth the effort...
            arr3out(:,:,dxa:dxb)=missval
        CASE(5)
            955 FORMAT (A37,I3,A4,I3)
            WRITE(*,955) "Unrecoverable shift from ", dxa, " to ", dxb
            arr3out(:,:,dxa:dxb)=missval
        CASE(6)
            PRINT*,"Vertical axis inverted - Nothing to be done"
            ! Nothing to be done
        CASE(7)
            PRINT*,"Adjust time axis"
        CASE(8)
            PRINT*,"Geopotential data missing - Nothing to be done"
            ! Nothing to be done
        CASE(9)
            PRINT*,"File missing"
        CASE(10)
            PRINT*,"Only file header - Nothing to be done"
            arr3out(:,:,:)=missval
        CASE DEFAULT
            945 FORMAT (A20,I5,A12,I4,I2)
            WRITE(*,945) "Unknown error code (",reason,"), for file ",year,mn
    END SELECT ! Select case depending on error code

    ENDDO ! do iy - number of corrections

    DEALLOCATE(arr3)

END SUBROUTINE ! corrgridnonqph

! =====================================================================================

SUBROUTINE outfilemodder(outfile,lvl)

    IMPLICIT NONE
    
    CHARACTER(100),INTENT(INOUT)                :: outfile
    INTEGER,INTENT(IN)                          :: lvl
    INTEGER                                     :: ii
    CHARACTER(100)                              :: out1,out2,clvl

    ii=SCAN(outfile,".",.TRUE.)
    out1=outfile(1:ii); out2=outfile(ii:len(trim(outfile)))
    WRITE(clvl,'(I5)') lvl
    outfile=trim(out1)//trim(adjustl(clvl))
    outfile=trim(outfile)//trim(out2)

    WRITE(*,*) "Prepare to write ensemble number ",trim(clvl)

END SUBROUTINE outfilemodder

! =====================================================================================

SUBROUTINE arith(arrin1,arrin2,mode,arr)
    REAL, DIMENSION(:,:,:,:),INTENT(IN)            :: arrin1,arrin2
    REAL, DIMENSION(:,:,:,:),INTENT(INOUT)         :: arr
    CHARACTER(32),INTENT(IN)                       :: mode

    IF (mode=="plus") THEN
        arr=arrin1+arrin2
        IF (dbg) WRITE(*,*) "Mode plus: arrin1+arrin2" 
    ELSE IF (mode=="minus") THEN
        arr=arrin1-arrin2
        IF (dbg) WRITE(*,*) "Mode minus: arrin1-arrin2"
    ELSE
        WRITE(*,*) "Arithmetic mode not supported/found. Arr=0"
        arr=0 
    ENDIF

END SUBROUTINE arith

! =====================================================================================

FUNCTION Compare_Float( x, y, ulp ) RESULT( Compare )
   REAL,              INTENT( IN )  :: x
   REAL,              INTENT( IN )  :: y
   INTEGER, OPTIONAL, INTENT( IN )  :: ulp
   LOGICAL :: Compare
   REAL :: Rel

   Rel = 1.0
   IF ( PRESENT( ulp ) ) THEN
       Rel = REAL( ABS(ulp) )
   END IF
   Compare = ABS( x - y ) < ( Rel * SPACING( MAX(ABS(x),ABS(y)) ) )

!   ULP:         Unit of data precision. The acronym stands for "unit in
!                the last place," the smallest possible increment or decrement
!                that can be made using a machine's floating point arithmetic.
!                A 0.5 ulp maximum error is the best you could hope for, since
!                this corresponds to always rounding to the nearest representable
!                floating-point number. Value must be positive - if a negative
!                negative value is supplied, the absolute value is used.
!                If not specified, the default value is 1.
!
!   The test performed is
!
!     ABS( x - y ) < ( ULP * SPACING( MAX(ABS(x),ABS(y)) ) )
!
!   If the result is .TRUE., the numbers are considered equal.
!
!   The intrinsic function SPACING(x) returns the absolute spacing of numbers
!   near the value of x,
!
!                  {     EXPONENT(x)-DIGITS(x)
!                  {  2.0                        for x /= 0
!     SPACING(x) = {
!                  {
!                  {  TINY(x)                    for x == 0
!
!   The ULP optional argument scales the comparison.
!
!   James Van Buskirk and James Giles suggested this method for floating
!   point comparisons in the comp.lang.fortran newsgroup. 
END FUNCTION Compare_Float 


END PROGRAM netcdfchecker
