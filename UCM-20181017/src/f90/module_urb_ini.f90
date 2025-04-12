module module_urb_ini
implicit none
   INTEGER                         :: ICATE
   REAL(8), ALLOCATABLE, DIMENSION(:) :: ZR_TBL
   REAL(8), ALLOCATABLE, DIMENSION(:) :: ROOF_WIDTH_TBL
   REAL(8), ALLOCATABLE, DIMENSION(:) :: ROAD_WITH_TBL
   REAL(8), ALLOCATABLE, DIMENSION(:) :: AH_TBL
   REAL(8), ALLOCATABLE, DIMENSION(:) :: FRC_URB_TBL
   REAL(8), ALLOCATABLE, DIMENSION(:) :: CAPR_TBL
   REAL(8), ALLOCATABLE, DIMENSION(:) :: CAPB_TBL
   REAL(8), ALLOCATABLE, DIMENSION(:) :: CAPG_TBL
   REAL(8), ALLOCATABLE, DIMENSION(:) :: AKSR_TBL
   REAL(8), ALLOCATABLE, DIMENSION(:) :: AKSB_TBL
   REAL(8), ALLOCATABLE, DIMENSION(:) :: AKSG_TBL
   REAL(8), ALLOCATABLE, DIMENSION(:) :: ALBR_TBL
   REAL(8), ALLOCATABLE, DIMENSION(:) :: ALBB_TBL
   REAL(8), ALLOCATABLE, DIMENSION(:) :: ALBG_TBL
   REAL(8), ALLOCATABLE, DIMENSION(:) :: Z0B_TBL
   REAL(8), ALLOCATABLE, DIMENSION(:) :: Z0G_TBL
   REAL(8), ALLOCATABLE, DIMENSION(:) :: Z0R_TBL
   REAL(8), ALLOCATABLE, DIMENSION(:) :: TRLEND_TBL
   REAL(8), ALLOCATABLE, DIMENSION(:) :: TBLEND_TBL
   REAL(8), ALLOCATABLE, DIMENSION(:) :: TGLEND_TBL
   REAL(8), ALLOCATABLE, DIMENSION(:) :: DDZR_TBL
   REAL(8), ALLOCATABLE, DIMENSION(:) :: DDZB_TBL
   REAL(8), ALLOCATABLE, DIMENSION(:) :: DDZG_TBL
   REAL(8), DIMENSION(24)             :: AHDIUPRF_TBL
   INTEGER         :: THERMAL_INSOL_ROOF_TBL
   INTEGER         :: THERMAL_INSOL_WALL_TBL
   INTEGER         :: RLNU_TBL
   INTEGER         :: BLNU_TBL
   INTEGER         :: GLNU_TBL
   INTEGER         :: BOUNDR_TBL
   INTEGER         :: BOUNDB_TBL
   INTEGER         :: BOUNDG_TBL
   INTEGER         :: AHOPTION_TBL
   INTEGER         :: FRCURBOPT_TBL

contains


!======================================================================
  subroutine read_urbparm ()
!======================================================================

  implicit none
  integer  :: IOSTATUS, allocate_status 
  integer  :: nuvar
  character(len=300)  :: string   
  character(len=50) :: name
  integer :: indx

   open ( 11, file='URBPARM.TBL', access='SEQUENTIAL',        &
          status='OLD', action='READ', position='REWIND',     &
          iostat=IOSTATUS )

   if( IOSTATUS > 0 ) print*, '   ERROR: open URBPARM.TBL'


  nuvar = 0
  READLOOP: do 
     
      read(11,'(A300)', iostat=iostatus) string
      if (iostatus /= 0) exit READLOOP
      
      if (string(1:1) == "#") cycle READLOOP
      if (trim(string) == "") cycle READLOOP
      indx = index(string,":")
      if (indx <= 0) cycle READLOOP
      name = trim(adjustl(string(1:indx-1)))
      
      
      IF (name == "Number of urban categories") then
         read(string(indx+1:),*) icate
!write(991,*) icate
!IF (.not. ALLOCATED(ZR_TBL)) then

         ALLOCATE( ZR_TBL(ICATE), stat=allocate_status )
           if(allocate_status /= 0) print*, 'Error allocating ZR_TBL in urban_param_init'
         ALLOCATE( ROOF_WIDTH_TBL(ICATE), stat=allocate_status )
           if(allocate_status /= 0) print*, 'Error allocating ROOF_WIDTH_TBL in urban_param_init'
         ALLOCATE( ROAD_WITH_TBL(ICATE), stat=allocate_status )
           if(allocate_status /= 0) print*, 'Error allocating ROAD_WITH_TBL in urban_param_init'
         ALLOCATE( AH_TBL(ICATE), stat=allocate_status )
           if(allocate_status /= 0) print*, 'Error allocating AH_TBL in urban_param_init'
         ALLOCATE( FRC_URB_TBL(ICATE), stat=allocate_status )
           if(allocate_status /= 0) print*, 'Error allocating FRC_URB_TBL in urban_param_init'
         ALLOCATE( CAPR_TBL(ICATE), stat=allocate_status )
           if(allocate_status /= 0) print*, 'Error allocating CAPR_TBL in urban_param_init'
         ALLOCATE( CAPB_TBL(ICATE), stat=allocate_status )
           if(allocate_status /= 0) print*, 'Error allocating CAPB_TBL in urban_param_init'
         ALLOCATE( CAPG_TBL(ICATE), stat=allocate_status )
           if(allocate_status /= 0) print*, 'Error allocating CAPG_TBL in urban_param_init'
         ALLOCATE( AKSR_TBL(ICATE), stat=allocate_status )
           if(allocate_status /= 0) print*, 'Error allocating AKSR_TBL in urban_param_init'
         ALLOCATE( AKSB_TBL(ICATE), stat=allocate_status )
           if(allocate_status /= 0) print*, 'Error allocating AKSB_TBL in urban_param_init'
         ALLOCATE( AKSG_TBL(ICATE), stat=allocate_status )
           if(allocate_status /= 0) print*, 'Error allocating AKSG_TBL in urban_param_init'
         ALLOCATE( ALBR_TBL(ICATE), stat=allocate_status )
           if(allocate_status /= 0) print*, 'Error allocating ALBR_TBL in urban_param_init'
         ALLOCATE( ALBB_TBL(ICATE), stat=allocate_status )
           if(allocate_status /= 0) print*, 'Error allocating ALBB_TBL in urban_param_init'
         ALLOCATE( ALBG_TBL(ICATE), stat=allocate_status )
           if(allocate_status /= 0) print*, 'Error allocating ALBG_TBL in urban_param_init'
         ALLOCATE( Z0B_TBL(ICATE), stat=allocate_status )
           if(allocate_status /= 0) print*, 'Error allocating Z0B_TBL in urban_param_init'
         ALLOCATE( Z0G_TBL(ICATE), stat=allocate_status )
           if(allocate_status /= 0) print*, 'Error allocating Z0G_TBL in urban_param_init'
         ALLOCATE( Z0R_TBL(ICATE), stat=allocate_status )
           if(allocate_status /= 0) print*, 'Error allocating Z0R_TBL in urban_param_init'
         ALLOCATE( TRLEND_TBL(ICATE), stat=allocate_status )
           if(allocate_status /= 0) print*, 'Error allocating TRLEND_TBL in urban_param_init'
         ALLOCATE( TBLEND_TBL(ICATE), stat=allocate_status )
           if(allocate_status /= 0) print*, 'Error allocating TBLEND_TBL in urban_param_init'
         ALLOCATE( TGLEND_TBL(ICATE), stat=allocate_status )
           if(allocate_status /= 0) print*, 'Error allocating TGLEND_TBL in urban_param_init'

         ALLOCATE( DDZR_TBL(ICATE), stat=allocate_status )
           if(allocate_status /= 0) print*, 'Error allocating DDZR_TBL in urban_param_init'
         ALLOCATE( DDZB_TBL(ICATE), stat=allocate_status )
           if(allocate_status /= 0) print*, 'Error allocating DDZB_TBL in urban_param_init'
         ALLOCATE( DDZG_TBL(ICATE), stat=allocate_status )
           if(allocate_status /= 0) print*, 'Error allocating DDZG_TBL in urban_param_init'


!ENDIF
      else if (name == "BUILDING_HEIGHT") then   
          read(string(indx+1:),*,iostat=iostatus) ZR_TBL(1:icate)
!write(991,'(a40,3f7.2)')  "BUILDING_HEIGHT :    ", ZR_TBL(:)
          if (iostatus /= 0) then
              print*, "ERROR when READING: BUILDING_HEIGHT " 
              stop
          end if
      else if (name == "ROOF_WIDTH") then   
          read(string(indx+1:),*,iostat=iostatus) ROOF_WIDTH_TBL(1:icate)
!write(991,'(a40,3f7.2)')  "ROOF_WIDTH :    ", ROOF_WIDTH_TBL(:)
          if (iostatus /= 0) then
              print*, "ERROR when READING: ROOF_WIDTH " 
              stop
          end if
      else if (name == "ROAD_WIDTH") then   
          read(string(indx+1:),*,iostat=iostatus) ROAD_WITH_TBL(1:icate)
!write(991,'(a40,3f7.2)')  "ROAD_WIDTH :    ", ROAD_WITH_TBL(:)
          if (iostatus /= 0) then
              print*, "ERROR when READING: ROAD_WIDTH " 
              stop
          end if
      else if (name == "AH") then   
          read(string(indx+1:),*,iostat=iostatus) AH_TBL(1:icate)
!write(991,'(a40,3f7.2)')  "AH :    ", AH_TBL(:)
          if (iostatus /= 0) then
              print*, "ERROR when READING: AH " 
              stop
          end if
      else if (name == "FRC_URB") then   
          read(string(indx+1:),*,iostat=iostatus) FRC_URB_TBL(1:icate)
!write(991,'(a40,3f7.2)')  "FRC_URB :    ", FRC_URB_TBL(:)
          if (iostatus /= 0) then
              print*, "ERROR when READING: FRC_URB " 
              stop
          end if
      else if (name == "HEAT_CAPACITY_ROOF") then   
          read(string(indx+1:),*,iostat=iostatus) CAPR_TBL(1:icate)
!write(991,'(a40,3f9.0)')  "HEAT_CAPACITY_ROOF :    ", CAPR_TBL(:)
          if (iostatus /= 0) then
              print*, "ERROR when READING: HEAT_CAPACITY_ROOF " 
              stop
          end if
      else if (name == "HEAT_CAPACITY_WALL") then   
          read(string(indx+1:),*,iostat=iostatus) CAPB_TBL(1:icate)
!write(991,'(a40,3f9.0)')  "HEAT_CAPACITY_WALL :    ", CAPB_TBL(:)
          if (iostatus /= 0) then
              print*, "ERROR when READING: HEAT_CAPACITY_WALL " 
              stop
          end if
      else if (name == "HEAT_CAPACITY_GROUND") then   
          read(string(indx+1:),*,iostat=iostatus) CAPG_TBL(1:icate)
!write(991,'(a40,3f9.0)')  "HEAT_CAPACITY_GROUND :    ", CAPG_TBL(:)
          if (iostatus /= 0) then
              print*, "ERROR when READING: HEAT_CAPACITY_GROUND " 
              stop
          end if
      else if (name == "THERMAL_CONDUCTIVITY_ROOF") then   
          read(string(indx+1:),*,iostat=iostatus) AKSR_TBL(1:icate)
!write(991,'(a40,3f7.2)')  "THERMAL_CONDUCTIVITY_ROOF :    ", AKSR_TBL(:)
          if (iostatus /= 0) then
              print*, "ERROR when READING: THERMAL_CONDUCTIVITY_ROOF " 
              stop
          end if
      else if (name == "THERMAL_CONDUCTIVITY_WALL") then   
          read(string(indx+1:),*,iostat=iostatus) AKSB_TBL(1:icate)
!write(991,'(a40,3f7.2)')  "THERMAL_CONDUCTIVITY_WALL :    ", AKSB_TBL(:)
          if (iostatus /= 0) then
              print*, "ERROR when READING: THERMAL_CONDUCTIVITY_WALL " 
              stop
          end if
      else if (name == "THERMAL_CONDUCTIVITY_GROUND") then   
          read(string(indx+1:),*,iostat=iostatus) AKSG_TBL(1:icate)
!write(991,'(a40,3f7.2)')  "THERMAL_CONDUCTIVITY_GROUND :    ", AKSG_TBL(:)
          if (iostatus /= 0) then
              print*, "ERROR when READING: THERMAL_CONDUCTIVITY_GROUND " 
              stop
          end if
      else if (name == "ALBEDO_ROOF") then   
          read(string(indx+1:),*,iostat=iostatus) ALBR_TBL(1:icate)
!write(991,'(a40,3f7.2)')  "ALBEDO_ROOF :    ", ALBR_TBL(:)
          if (iostatus /= 0) then
              print*, "ERROR when READING: ALBEDO_ROOF " 
              stop
          end if
      else if (name == "ALBEDO_WALL") then   
          read(string(indx+1:),*,iostat=iostatus) ALBB_TBL(1:icate)
!write(991,'(a40,3f7.2)')  "ALBEDO_WALL :    ", ALBB_TBL(:)
          if (iostatus /= 0) then
              print*, "ERROR when READING: ALBEDO_WALL " 
              stop
          end if
      else if (name == "ALBEDO_GROUND") then   
          read(string(indx+1:),*,iostat=iostatus) ALBG_TBL(1:icate)
!write(991,'(a40,3f7.2)')  "ALBEDO_GROUND :    ", ALBG_TBL(:)
          if (iostatus /= 0) then
              print*, "ERROR when READING: ALBEDO_GROUND " 
              stop
          end if
      else if (name == "Z0_WALL") then   
          read(string(indx+1:),*,iostat=iostatus) Z0B_TBL(1:icate)
!write(991,'(a40,3f7.2)')  "Z0_WALL :    ", Z0B_TBL(:)
          if (iostatus /= 0) then
              print*, "ERROR when READING: Z0_WALL " 
              stop
          end if
      else if (name == "Z0_GROUND") then   
          read(string(indx+1:),*,iostat=iostatus) Z0G_TBL(1:icate)
!write(991,'(a40,3f7.2)')  "Z0_GROUND :    ", Z0G_TBL(:)
          if (iostatus /= 0) then
              print*, "ERROR when READING: Z0_GROUND " 
              stop
          end if
      else if (name == "Z0_ROOF") then   
          read(string(indx+1:),*,iostat=iostatus) Z0R_TBL(1:icate)
!write(991,'(a40,3f7.2)')  "Z0_ROOF :    ", Z0R_TBL(:)
          if (iostatus /= 0) then
              print*, "ERROR when READING: Z0_ROOF " 
              stop
          end if
      else if (name == "T_IN_ROOF") then   
          read(string(indx+1:),*,iostat=iostatus) TRLEND_TBL(1:icate)
!write(991,'(a40,3f7.2)')  "T_IN_ROOF :    ", TRLEND_TBL(:)
          if (iostatus /= 0) then
              print*, "ERROR when READING: T_IN_ROOF " 
              stop
          end if
      else if (name == "T_IN_WALL") then   
          read(string(indx+1:),*,iostat=iostatus) TBLEND_TBL(1:icate)
!write(991,'(a40,3f7.2)')  "T_IN_WALL :    ", TBLEND_TBL(:)
          if (iostatus /= 0) then
              print*, "ERROR when READING: T_IN_WALL " 
              stop
          end if
      else if (name == "T_IN_GROUND") then   
          read(string(indx+1:),*,iostat=iostatus) TGLEND_TBL(1:icate)
!write(991,'(a40,3f7.2)')  "T_IN_GROUND :    ", TGLEND_TBL(:)
          if (iostatus /= 0) then
              print*, "ERROR when READING: T_IN_GROUND " 
              stop
          end if
 

      else if (name == "ROOF_THICKNESS") then   
          read(string(indx+1:),*,iostat=iostatus) DDZR_TBL(1:icate)
!write(991,'(a30,3f7.2)')  "ROOF_THICKNESS :    ", DDZR_TBL(:)
          if (iostatus /= 0) then
              print*, "ERROR when READING: ROOF_THICKNESS " 
              stop
          end if
      else if (name == "WALL_THICKNESS") then   
          read(string(indx+1:),*,iostat=iostatus) DDZB_TBL(1:icate)
!write(991,'(a30,3f7.2)')  "WALL_THICKNESS :    ", DDZB_TBL(:)
          if (iostatus /= 0) then
              print*, "ERROR when READING: WALL_THICKNESS " 
              stop
          end if
      else if (name == "GROUND_THICKNESS") then   
          read(string(indx+1:),*,iostat=iostatus) DDZG_TBL(1:icate)
!write(991,'(a30,3f7.2)')  "GROUND_THICKNESS :    ", DDZG_TBL(:)
          if (iostatus /= 0) then
              print*, "ERROR when READING: GROUND_THICKNESS " 
              stop
          end if
      else if (name == "AH_DIURNAL_PROFILE") then   
          read(string(indx+1:),*,iostat=iostatus) AHDIUPRF_TBL(1:24)
!write(991,'(a30,24f4.1)')  "AH_DIURNAL_PROFILE :    ", AHDIUPRF_TBL(:)
          if (iostatus /= 0) then
              print*, "ERROR when READING: AH_DIURNAL_PROFILE " 
              stop
          end if




!---------------------------------------

        else if (name == "THERMAL_INSOL_ROOF") then   
          read(string(indx+1:),*,iostat=iostatus) THERMAL_INSOL_ROOF_TBL
!write(991,'(a30,i4)')  "THERMAL_INSOL_ROOF :    ", THERMAL_INSOL_ROOF_TBL
          if (iostatus /= 0) then
              print*, "ERROR when READING: THERMAL_INSOL_ROOF " 
              stop
          end if
      else if (name == "THERMAL_INSOL_WALL") then   
          read(string(indx+1:),*,iostat=iostatus) THERMAL_INSOL_WALL_TBL
!write(991,'(a30,i4)')  "THERMAL_INSOL_WALL :    ", THERMAL_INSOL_WALL_TBL
          if (iostatus /= 0) then
              print*, "ERROR when READING: THERMAL_INSOL_WALL " 
              stop
          end if
      else if (name == "ROOF_LAYERS_NU") then   
          read(string(indx+1:),*,iostat=iostatus) RLNU_TBL
!write(991,'(a30,i4)')  "ROOF_LAYERS_NU :    ", RLNU_TBL
          if (iostatus /= 0) then
              print*, "ERROR when READING: ROOF_LAYERS_NU " 
              stop
          end if
      else if (name == "WALL_LAYERS_NU") then   
          read(string(indx+1:),*,iostat=iostatus) BLNU_TBL
!write(991,'(a30,i4)')  "WALL_LAYERS_NU :    ", BLNU_TBL
          if (iostatus /= 0) then
              print*, "ERROR when READING: WALL_LAYERS_NU " 
              stop
          end if
      else if (name == "GROUND_LAYERS_NU") then   
          read(string(indx+1:),*,iostat=iostatus) GLNU_TBL
!write(991,'(a30,i4)')  "GROUND_LAYERS_NU :    ", GLNU_TBL
          if (iostatus /= 0) then
              print*, "ERROR when READING: GROUND_LAYERS_NU " 
              stop
          end if
      else if (name == "BOUND_CONDITION_ROOF") then   
          read(string(indx+1:),*,iostat=iostatus) BOUNDR_TBL
!write(991,'(a30,i4)')  "BOUND_CONDITION_ROOF :    ", BOUNDR_TBL
          if (iostatus /= 0) then
              print*, "ERROR when READING: BOUND_CONDITION_ROOF " 
              stop
          end if
      else if (name == "BOUND_CONDITION_WALL") then   
          read(string(indx+1:),*,iostat=iostatus) BOUNDB_TBL
!write(991,'(a30,i4)')  "BOUND_CONDITION_WALL :    ", BOUNDB_TBL
          if (iostatus /= 0) then
              print*, "ERROR when READING: BOUND_CONDITION_WALL " 
              stop
          end if
      else if (name == "BOUND_CONDITION_GROUND") then   
          read(string(indx+1:),*,iostat=iostatus) BOUNDG_TBL
!write(991,'(a30,i4)')  "BOUND_CONDITION_GROUND :    ", BOUNDG_TBL
          if (iostatus /= 0) then
              print*, "ERROR when READING: BOUND_CONDITION_GROUND " 
              stop
          end if
      else if (name == "AH_OPTION") then   
          read(string(indx+1:),*,iostat=iostatus) AHOPTION_TBL
!write(991,'(a30,i4)')  "AH_OPTION :    ", AHOPTION_TBL
          if (iostatus /= 0) then
              print*, "ERROR when READING: AH_OPTION " 
              stop
          end if
     
      else if (name == "FRC_URB_OPT") then   
          read(string(indx+1:),*,iostat=iostatus) FRCURBOPT_TBL
!write(991,'(a30,i4)')  "FRC_URB_OPT :    ", FRCURBOPT_TBL
          if (iostatus /= 0) then
              print*, "ERROR when READING: FRC_URB_OPT " 
              stop
          end if



     ENDIF


        nuvar = nuvar + 1

     ENDDO READLOOP


     string = "end"

!print*, nuvar

       if( nuvar /= 35 ) then 
           print*, '   ERROR: READING URBPARM.TBL'
           stop
       end if

  close(11)


 return
 end subroutine read_urbparm
!============================================================






!===========================================================!
subroutine zero_clear_tbl ( )


    ZR_TBL(:) = 0.d0
    ROOF_WIDTH_TBL(:) = 0.d0
    ROAD_WITH_TBL(:) = 0.d0
    AH_TBL(:) = 0.d0
    FRC_URB_TBL(:) = 0.d0
    CAPR_TBL(:) = 0.d0
    CAPB_TBL(:) = 0.d0
    CAPG_TBL(:) = 0.d0
    AKSR_TBL(:) = 0.d0
    AKSB_TBL(:) = 0.d0
    AKSG_TBL(:) = 0.d0
    ALBR_TBL(:) = 0.d0
    ALBB_TBL(:) = 0.d0
    ALBG_TBL(:) = 0.d0
    Z0B_TBL(:) = 0.d0
    Z0G_TBL(:) = 0.d0
    Z0R_TBL(:) = 0.d0
    TRLEND_TBL(:) = 0.d0
    TBLEND_TBL(:) = 0.d0
    TGLEND_TBL(:) = 0.d0
    DDZR_TBL(:) = 0.d0
    DDZB_TBL(:) = 0.d0
    DDZG_TBL(:) = 0.d0
    THERMAL_INSOL_ROOF_TBL = 0
    THERMAL_INSOL_WALL_TBL = 0
    RLNU_TBL = 0
    BLNU_TBL = 0
    GLNU_TBL = 0
    BOUNDR_TBL = 0
    BOUNDB_TBL = 0
    BOUNDG_TBL = 0
    AHOPTION_TBL = 0
    AHDIUPRF_TBL = 0
    FRCURBOPT_TBL = 0

return
end subroutine zero_clear_tbl



end module 
!=====================================









