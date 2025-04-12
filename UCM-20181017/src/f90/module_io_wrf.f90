module module_io_wrf

 use module_params, only: wrf_infile, WRF_POINT, start_year
 use netcdf
 USE module_vars, only : GRDP_WRF
!=============
 contains
!=============


SUBROUTINE READ_INTERP_WRF ( THETA_WRF, QV_WRF, U_WRF, V_WRF, &
                             SWDOWN_WRF, GLW_WRF,             &
                             COSZ_WRF, OMG_WRF, DECLIN_WRF      )
  implicit none
    
    real(8), dimension(:), allocatable, INTENT(OUT)   :: THETA_WRF, QV_WRF, U_WRF, V_WRF 
    real(8), dimension(:), allocatable, INTENT(OUT)   :: SWDOWN_WRF, GLW_WRF
    real(8), dimension(:), allocatable, INTENT(OUT)   :: COSZ_WRF, OMG_WRF, DECLIN_WRF 

!TYPE wrflu
!INTEGER          :: LU_INDEX
!INTEGER          :: UTYPE_URB
!REAL(8)          :: FRC_URB
!REAL(8)          :: ALBEDO
!REAL(8)          :: EMISS
!REAL(8)          :: XLAT
!REAL(8)          :: XLONG
!REAL(8)          :: HGT
!END TYPE wrflu
!TYPE(wrflu)      ::  GRDP_WRF

    real(8), dimension(:), allocatable   :: HEIGHT_WRF   
    character (len = 19), dimension(:), allocatable  :: TIME_WRF


    real(8), dimension(:,:,:), allocatable :: THETA, QV, U, V 
    real(8), dimension(:,:), allocatable   :: SWDOWN, GLW, &
                                              COSZ_URB, OMG_URB 
    real(8), dimension(:), allocatable     :: HEIGHT   
    character (len = 19), dimension(:), allocatable  :: TIME

    INTEGER :: I1, I2, I3, IS, IH, IT
    REAL(8), dimension(:), allocatable  ::     &
                   LU_INDEX, FRC_URB, UTYPE_URB, ALBEDO,       &
                   EMISS, XLAT, XLONG, HGT, DECLIN_URB

!-----------------------------------------------------------!




    CALL READ_WRF ( THETA, QV, U, V,             &
                      SWDOWN, GLW, HEIGHT, TIME, &
                      LU_INDEX, FRC_URB,         &
                      UTYPE_URB, ALBEDO, EMISS,  &
                      XLAT, XLONG, HGT,          &
                      COSZ_URB, OMG_URB, DECLIN_URB )



    IT = 25 - 7  ! start time

    I1 = Ubound(THETA, 1)   ! STATION
    I2 = Ubound(THETA, 2)   ! HEIGHT
    I3 = Ubound(THETA, 3) - IT ! TIME
!print*, i1, i2, i3
    ALLOCATE( THETA_WRF(I3) )
    ALLOCATE( QV_WRF(I3) )
    ALLOCATE( U_WRF(I3) )
    ALLOCATE( V_WRF(I3) )
    ALLOCATE( SWDOWN_WRF(I3) )
    ALLOCATE( GLW_WRF(I3) )
    ALLOCATE( TIME_WRF(I3) )
    ALLOCATE( COSZ_WRF(I3) )
    ALLOCATE( OMG_WRF(I3) )
    ALLOCATE( DECLIN_WRF(I3) )
   


    IH = 5          ! height 5
    IS = WRF_POINT  ! station
    IT = 25 - 7     ! start time

     THETA_WRF  =  THETA      ( IS, IH, IT:IT+I3 )
     QV_WRF     =  QV         ( IS, IH, IT:IT+I3 )
     U_WRF      =  U          ( IS, IH, IT:IT+I3 )
     V_WRF      =  V          ( IS, IH, IT:IT+I3 )
     SWDOWN_WRF =  SWDOWN     ( IS,     IT:IT+I3 )
     GLW_WRF    =  GLW        ( IS,     IT:IT+I3 )
     COSZ_WRF   =  COSZ_URB   ( IS,     IT:IT+I3 )
     OMG_WRF    =  OMG_URB    ( IS,     IT:IT+I3 )
     DECLIN_WRF =  DECLIN_URB (         IT:IT+I3 )
!TIME_WRF = TIME
    
     GRDP_WRF%LU_INDEX  = INT(LU_INDEX( IS ))
     GRDP_WRF%UTYPE_URB = INT(UTYPE_URB( IS ))
     GRDP_WRF%FRC_URB   = FRC_URB( IS )
     GRDP_WRF%ALBEDO    = ALBEDO( IS )
     GRDP_WRF%EMISS     = EMISS( IS )
     GRDP_WRF%XLAT      = XLAT( IS )
     GRDP_WRF%XLONG     = XLONG( IS )
     GRDP_WRF%HGT       = HGT( IS )
     GRDP_WRF%HEIGHT_TOP = HEIGHT( IH )     
END SUBROUTINE 
!------------------------------------------------------


















SUBROUTINE INTERP_WRF_DATA (ISTEP, N,                        & !-IN
                            DT, DT_OBS,                      & !-IN
                            THETA_WRF, QV_WRF, U_WRF, V_WRF, & !-IN
                            SWDOWN_WRF, GLW_WRF,             & !-IN
                            COSZ_WRF, OMG_WRF, DECLIN_WRF,   & !-IN
                            THETA_T, QV_T, U_T, V_T,         & !-OUT
                            SW_B, LW_B,                      & !-OUT
                            COSZ, OMG, DECLIN )                !-OUT

! This subroutine to intepolate the wrf input values to time step
!
  IMPLICIT NONE
    INTEGER, INTENT(IN)                 :: ISTEP, N
    real(8), dimension(N), INTENT(IN)   :: THETA_WRF, QV_WRF, U_WRF, V_WRF 
    real(8), dimension(N), INTENT(IN)   :: SWDOWN_WRF, GLW_WRF  
    real(8), dimension(N), INTENT(IN)   :: COSZ_WRF, OMG_WRF, DECLIN_WRF  
    REAL(8), INTENT(IN)                 :: DT, DT_OBS

    REAL(8), INTENT(OUT)                :: THETA_T, QV_T, U_T, V_T,  &
                                           SW_B, LW_B,               &
                                           COSZ, OMG, DECLIN

! Local variable
    INTEGER                      :: TSTEP
    real(8)                      :: dtt, time, PI
!----------------------------------
    PI = 4.D0 * DATAN(1.D0)
   

    time     =  dfloat(istep)*dt/dt_obs + 1.d0      ! time in real (h)
    TSTEP    =  int(time) 
    dtt = time - dfloat(TSTEP) 
    IF(TSTEP==n) STOP 'ERROR: number of run step is invalid'
    THETA_T = THETA_WRF(TSTEP)  + (THETA_WRF(TSTEP+1)  - THETA_WRF(TSTEP))   * dtt  
    QV_T    = QV_WRF(TSTEP)     + (QV_WRF(TSTEP+1)     - QV_WRF(TSTEP))      * dtt
    U_T     = U_WRF(TSTEP)      + (U_WRF(TSTEP+1)      - U_WRF(TSTEP))       * dtt
    V_T     = V_WRF(TSTEP)      + (V_WRF(TSTEP+1)      - V_WRF(TSTEP))       * dtt
    SW_B    = SWDOWN_WRF(TSTEP) + (SWDOWN_WRF(TSTEP+1) - SWDOWN_WRF(TSTEP))  * dtt
    LW_B    = GLW_WRF(TSTEP)    + (GLW_WRF(TSTEP+1)    - GLW_WRF(TSTEP))     * dtt
    COSZ    = COSZ_WRF(TSTEP)   + (COSZ_WRF(TSTEP+1)   - COSZ_WRF(TSTEP))     * dtt
    DECLIN  = DECLIN_WRF(TSTEP) + (DECLIN_WRF(TSTEP+1) - DECLIN_WRF(TSTEP))  * dtt

    IF(OMG_WRF(TSTEP)*OMG_WRF(TSTEP+1)<0.D0) THEN
    OMG     = OMG_WRF(TSTEP)    + (15.D0*PI/180.D0) * DTT
    ELSE
    OMG     = OMG_WRF(TSTEP)    + (OMG_WRF(TSTEP+1)    - OMG_WRF(TSTEP))     * DTT
    END IF
!if(mod(istep,900)==0) write(*,'(f5.1,i4,3f7.2)'), TIME, TSTEP, omg, omg_WRF(TSTEP), omg_WRF(TSTEP+1)


END SUBROUTINE


























SUBROUTINE READ_WRF ( THETA, QV, U, V,           &
                      SWDOWN, GLW, HEIGHT, TIME, &
                      LU_INDEX, FRC_URB,         &
                      UTYPE_URB, ALBEDO, EMISS,  &
                      XLAT, XLONG, HGT,          &
                      COSZ_URB, OMG_URB, DECLIN_URB )


  implicit none
    real(8), dimension(:,:,:), allocatable, intent(out) :: THETA, QV, U, V 
    real(8), dimension(:,:), allocatable, intent(out)   :: SWDOWN, GLW, &
                                                           COSZ_URB, OMG_URB

    real(8), dimension(:), allocatable, intent(out)     :: HEIGHT   
    character (len = 19),  dimension(:), allocatable, intent(out)  :: TIME
    REAL(8), dimension(:), allocatable, INTENT(OUT)     :: LU_INDEX,  &
                                                           FRC_URB,   &
                                                           UTYPE_URB, &
                                                           ALBEDO,    &
                                                           EMISS,     &
                                                           XLAT,      &
                                                           XLONG,     &
                                                           HGT, DECLIN_URB
    INTEGER   :: I, NF, IOS
    character (len = 200)    :: A1, A2, A3 
    INTEGER  :: id, error, p
    
 
    INTEGER      :: ncid_in
    INTEGER      :: irec, iv, iat
    INTEGER      :: numdims
    INTEGER      :: numatts
    INTEGER      :: nDimensions
    INTEGER      :: nVariables
    INTEGER      :: nAttributes
    INTEGER      :: unlimitedDimId
    INTEGER      :: formatNum 
    INTEGER, dimension(:), allocatable  :: lendim
    INTEGER, dimension(:), allocatable  :: odimid
    INTEGER, dimension(:), allocatable  :: IVAR_ID
    INTEGER, dimension(:), allocatable  :: start, count
    INTEGER, dimension(:), allocatable  :: dimids

    character (len = 100), dimension(:), allocatable  :: dimname, var_attname

    real(8), dimension(:,:,:), allocatable :: VAR3D
    real(8), dimension(:,:), allocatable   :: VAR2D


    character (len = 100) :: VAR_NAME
    character (len = 100) :: WRF_INPUT_FILE
    character (len = 4)   :: RUN_YEAR    
    CHARACTER (LEN = 100), DIMENSION(:), ALLOCATABLE :: VARNAMES
    INTEGER :: i1, i2, i3, i4 , nvar, VARID, ID1, ID2, ID3, ID4
!--------------------------------------




    OPEN(12,FILE="WRFVAR.TBL")
    NVAR = 0
    DO I = 1, 1000
      READ(12,*,IOSTAT=IOS) A1, A2
      IF(IOS /= 0) EXIT
      IF(I == 1000) THEN
         PRINT*, "ERROR: MAXIMUM NUMBER OF RECORDS EXCEEDED..."
         STOP
      END IF
    NVAR = NVAR + 1
   END DO
   CLOSE(12)
   ALLOCATE(VARNAMES(NVAR))


   OPEN(12,FILE="WRFVAR.TBL")
   DO I = 1, NVAR
      READ(12,*,IOSTAT=IOS) VARNAMES(I), A2
      IF(IOS /= 0) EXIT
   END DO
   CLOSE(12)
  
!PRINT*, VARNAMES
       write(run_year,      '(i4)'  ) start_year
       WRF_INPUT_FILE = "WRF/CU_"//RUN_YEAR//".nc"
       PRINT*, "   INPUT WRF FILE: ",  trim(wrf_infile) !TRIM(WRF_INPUT_FILE)
! OPEN FILE
       call check( nf90_open(trim(wrf_infile),nf90_nowrite, ncid_in) )
! INQUIRE DIMENSION, VARIABLE, ATTRIBUTES, ETC. INDEXES
       call check( nf90_inquire ( ncid_in, nDimensions, nVariables, nAttributes,  &
                                  unlimitedDimId, formatNum ))
! print*, nDimensions, nVariables, nAttributes, unlimitedDimId
       allocate(dimname(nDimensions))
       allocate(lendim(nDimensions))
! INQUIRE DIMENSION NAME AND LENGTH
       do id=1, nDimensions
        call check(nf90_inquire_dimension(ncid_in,id, dimname(id), lendim(id)))
!print*, dimname(id), lendim(id)
       end do

! READ VARIABLES
    do iv = 1, NVAR ! NVAR IS NUMBER OF VARIABLES THAT READED IN 12
! FIRST, INQUIRE VARIABLES INDEX
       CALL CHECK( NF90_INQ_VARID (NCID_IN, VARNAMES(IV), VARID) )
! CALL CHECK( NF90_INQ_DIMID (NCID_IN, VARNAMES(IV), DIMIDS
! SECOND, INQUIRE VARIABLES DIMENSION NUMBER
        call check( nf90_inquire_variable(ncid_in, VARID, A1, i1, i2))
! THIRD, AFTER GET THE DIMENSION SIZES, ALLOCATE THE DIM ID
        allocate(dimids(i2))
        call check( nf90_inquire_variable(ncid_in, VARID, A1, i1, i2, dimids))
! GET DIMENSION ID

! GET VARIABLES
! IN THE CASE OF 3 DIMENSION
         IF(I2==3) THEN
          ID1 = LENDIM(DIMIDS(1))
          ID2 = LENDIM(DIMIDS(2))
          ID3 = LENDIM(DIMIDS(3))
        
! FOR THETA
          IF( VARNAMES(IV)=="THETA" ) then
            ALLOCATE( THETA( ID1, ID2, ID3 ) )
            CALL CHECK ( NF90_GET_VAR(NCID_IN, VARID, THETA) )
!PRINT*, THETA
          end if
! FOR QV
          IF( VARNAMES(IV)=="QV" ) then
            ALLOCATE( QV( ID1, ID2, ID3 ) )
            CALL CHECK ( NF90_GET_VAR(NCID_IN, VARID, QV) )
          end if
        
! FOR U
          IF( VARNAMES(IV)=="U" ) then
            ALLOCATE( U( ID1, ID2, ID3 ) )
            CALL CHECK ( NF90_GET_VAR(NCID_IN, VARID, U) )
          end if
! FOR V
          IF( VARNAMES(IV)=="V" ) then
            ALLOCATE( V( ID1, ID2, ID3 ) )
            CALL CHECK ( NF90_GET_VAR(NCID_IN, VARID, V) )
          end if

         END IF  ! END IF DIMSIZE IS 3







! GET 2 DIMENSIONS VARIABLES: TIME,POINT
         IF(I2==2) THEN

           ID1 = LENDIM(DIMIDS(1))
           ID2 = LENDIM(DIMIDS(2))
! FOR SWDOWN
           IF( VARNAMES(IV)=="SWDOWN" ) then
             ALLOCATE( SWDOWN( ID1, ID2 ) )
             CALL CHECK ( NF90_GET_VAR(NCID_IN, VARID, SWDOWN) )
           end if
             
! FOR GLW
           IF( VARNAMES(IV)=="GLW" ) then
             ALLOCATE( GLW( ID1, ID2 ) )
             CALL CHECK ( NF90_GET_VAR(NCID_IN, VARID, GLW) )
           end if

! FOR COSZ_URB
           IF( VARNAMES(IV)=="COSZ_URB" ) then
             ALLOCATE( COSZ_URB( ID1, ID2 ) )
             CALL CHECK ( NF90_GET_VAR(NCID_IN, VARID, COSZ_URB) )
           end if

! FOR OMG_URB
           IF( VARNAMES(IV)=="OMG_URB" ) then
             ALLOCATE( OMG_URB( ID1, ID2 ) )
             CALL CHECK ( NF90_GET_VAR(NCID_IN, VARID, OMG_URB) )
           end if

! FOR TIMES
           IF( VARNAMES(IV)=="Times" ) then
             ALLOCATE( TIME( ID2 ) )
             CALL CHECK ( NF90_GET_VAR(NCID_IN, VARID, TIME) )
!print*, time
           end if

         END IF







! GET 2 DIMENSIONS VARIABLES: HEIGHT
         IF(I2==1) THEN
           ID1 = LENDIM(DIMIDS(1))
! FOR HEIGHT
           IF( VARNAMES(IV)=="Height" ) then
             ALLOCATE( HEIGHT( ID1) )
             CALL CHECK ( NF90_GET_VAR(NCID_IN, VARID, HEIGHT) )
           end if
         END IF

! FOR LU_INDEX
            IF( VARNAMES(IV)=="LU_INDEX" ) then
            ALLOCATE( LU_INDEX( ID1) )
            CALL CHECK ( NF90_GET_VAR(NCID_IN, VARID, LU_INDEX) )
            END IF

! FOR FRC_URB
            IF( VARNAMES(IV)=="FRC_URB" ) then
            ALLOCATE( FRC_URB( ID1) )
            CALL CHECK ( NF90_GET_VAR(NCID_IN, VARID, FRC_URB) )
            END IF

! FOR UTYPE_URB
            IF( VARNAMES(IV)=="UTYPE_URB" ) then
            ALLOCATE( UTYPE_URB( ID1) )
            CALL CHECK ( NF90_GET_VAR(NCID_IN, VARID, UTYPE_URB) )
            END IF

! FOR ALBEDO
            IF( VARNAMES(IV)=="ALBEDO" ) then
            ALLOCATE( ALBEDO( ID1) )
            CALL CHECK ( NF90_GET_VAR(NCID_IN, VARID, ALBEDO) )
            END IF

! FOR EMISS
            IF( VARNAMES(IV)=="EMISS" ) then
            ALLOCATE( EMISS( ID1) )
            CALL CHECK ( NF90_GET_VAR(NCID_IN, VARID, EMISS) )
            END IF

! FOR XLAT
            IF( VARNAMES(IV)=="XLAT" ) then
            ALLOCATE( XLAT( ID1) )
            CALL CHECK ( NF90_GET_VAR(NCID_IN, VARID, XLAT) )
            END IF

! FOR XLONG
            IF( VARNAMES(IV)=="XLONG" ) then
            ALLOCATE( XLONG( ID1) )
            CALL CHECK ( NF90_GET_VAR(NCID_IN, VARID, XLONG) )
            END IF

! FOR HGT
            IF( VARNAMES(IV)=="HGT" ) then
            ALLOCATE( HGT( ID1) )
            CALL CHECK ( NF90_GET_VAR(NCID_IN, VARID, HGT) )
            END IF

! FOR DECLIN
            IF( VARNAMES(IV)=="DECLIN_URB" ) then
            ALLOCATE( DECLIN_URB( ID1) )
            CALL CHECK ( NF90_GET_VAR(NCID_IN, VARID, DECLIN_URB) )
            END IF


  


       deallocate(dimids)
     end do
     deallocate(dimname)
     deallocate(lendim)
  
!print*, size(Hgt)
        print*, "---END READ NETCDF FILE"
END SUBROUTINE READ_WRF
























!=================================================!
 subroutine check(status)
    INTEGER, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check  
!=================================================!









end module module_io_wrf
