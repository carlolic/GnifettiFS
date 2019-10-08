FUNCTION Get_SurfaceEnthalpy(  Model, Node, InputArray) RESULT(SurfaceEnthalpy)

  USE DefUtils
  
  IMPLICIT NONE
  
  !----------------------------------------------------------------------------
  ! External variables
  !----------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  INTEGER :: Node
  REAL(KIND=dp) :: InputArray(2)
  REAL(KIND=dp) :: SurfaceEnthalpy
  
  !----------------------------------------------------------------------------
  ! Internal variables
  !----------------------------------------------------------------------------
  REAL(KIND=dp) ::  T_surf, Time, Tdev_yy, Tpos_yy  
  REAL(KIND=dp) ::  MeltwaterHeatFlux, MeltFactor, H_surf
  REAL(KIND=dp) ::  HeatCap_A, HeatCap_B, T_ref  
  
  INTEGER :: it, jt, iw, jw, lines_t, lines_w, io_error_t, io_error_w, dt_found, dw_found
  INTEGER, ALLOCATABLE :: yy_t(:), yy_w(:) 
  REAL(KIND=dp), ALLOCATABLE ::  dT(:), dW(:)

  LOGICAL :: GotIt
  TYPE(ValueList_t), POINTER :: BC
  

  ! Read variables  
  T_surf = InputArray(1)                   ! kelvin
  Time = InputArray(2)                     ! years  

  
  ! Read temperature variations  
  lines_t=0
  OPEN(unit=20,file='DATA/Tdev_yy.dat',status='old',action='read', iostat=io_error_t)
  DO 
    READ (20,*,END=10)
    lines_t = lines_t+1
  END DO
  10 CLOSE(unit=20) 
  
  ALLOCATE( yy_t(lines_t), dT(lines_t) )
  OPEN(unit=20,file='DATA/Tdev_yy.dat',status='old',action='read', iostat=io_error_t)
  DO it=1,lines_t 
    READ (20,*) yy_t(it), dT(it)
  END DO
  CLOSE(unit=20) 
    
  Tdev_yy = 0.0
  dt_found = 0
  DO jt=1,SIZE(yy_t)
    IF( INT(Time-1)==yy_t(jt) ) THEN
      IF( jt>0 ) Tdev_yy = dT(jt)
      dt_found = 1
    END IF  
  END DO
  IF( dt_found==0 ) PRINT *, 'No Tdev_yy found at given time. Tdev_yy=0 set as default'
  
  
  ! Read year-integrated positive temperatures
  lines_w=0
  OPEN(unit=30,file='DATA/Tpos_yy.dat',status='old',action='read', iostat=io_error_w)
  DO 
    READ (30,*,END=15)
    lines_w = lines_w+1
  END DO
  15 CLOSE(unit=30) 
  
  ALLOCATE( yy_w(lines_w), dW(lines_w) )
  OPEN(unit=30,file='DATA/Tpos_yy.dat',status='old',action='read', iostat=io_error_w)
  DO iw=1,lines_w 
    READ (30,*) yy_w(iw), dW(iw)
  END DO
  CLOSE(unit=30) 
    
  Tpos_yy = 0.0
  dw_found = 0
  DO jw=1,SIZE(yy_w)
    IF( INT(Time-1)==yy_w(jw) ) THEN
      IF( jw>0 ) Tpos_yy = dW(jw)
      dw_found = 1
    END IF  
  END DO
  IF( dw_found==0 ) PRINT *, 'No Tpos_yy found at given time. Tpos_yy=0 set as default'
  
  
  
  ! Read constants
  BC => GetBC()
  IF ( .NOT. ASSOCIATED(BC) ) THEN
     CALL FATAL('Get_SurfaceEnthalpy','No boundary condition found')
  END IF
  
  MeltFactor = GetConstReal( BC, 'Melt Factor', GotIt)
  IF (.NOT. GotIt) THEN
     CALL FATAL('Get_SurfaceEnthalpy', 'Melt Factor not found in BC section')
  END IF

  HeatCap_A = ListGetConstReal( Model % Constants, 'Enthalpy Heat Capacity A', GotIt )
  IF (.NOT. GotIt) THEN 
     CALL FATAL('Get_SurfaceEnthalpy', 'Enthalpy Heat Capacity A not found in constant section')
  END IF
  
  HeatCap_B = ListGetConstReal( Model % Constants, 'Enthalpy Heat Capacity B', GotIt )
  IF (.NOT. GotIt) THEN 
     CALL FATAL('Get_SurfaceEnthalpy', 'Enthalpy Heat Capacity B not found in constant section')
  END IF

  T_ref = ListGetConstReal( Model % Constants, 'T_ref_enthalpy', GotIt )
  IF (.NOT. GotIt) THEN 
     CALL FATAL('Get_SurfaceEnthalpy', 'T_ref_enthalpy not found in constant section')
  END IF
    
  
  MeltwaterHeatFlux = 0.0
    
  H_surf = HeatCap_A/2.*( (T_surf+Tdev_yy)**2 - T_ref**2 ) + HeatCap_B*( T_surf+Tdev_yy - T_ref )
 
  MeltwaterHeatFlux = MeltFactor * Tpos_yy
  
  SurfaceEnthalpy = H_surf + MeltwaterHeatFlux

  PRINT *, Time, Tdev_yy, Tpos_yy
  
  RETURN

END FUNCTION Get_SurfaceEnthalpy
