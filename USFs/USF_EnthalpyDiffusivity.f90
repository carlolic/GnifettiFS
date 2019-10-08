FUNCTION Get_EnthalpyDiffusivity(  Model, Node, InputArray) RESULT(EnthalpyDiffusivity)

  USE DefUtils
  
  IMPLICIT NONE
  
  !----------------------------------------------------------------------------
  ! External variables
  !----------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  INTEGER :: Node
  REAL(KIND=dp) :: InputArray(2)
  REAL(KIND=dp) :: EnthalpyDiffusivity 
  
  !----------------------------------------------------------------------------
  ! Internal variables
  !----------------------------------------------------------------------------
  REAL(KIND=dp) ::  Temperature, RelativeDensity, HeatCapacity, HeatConductivity, IceDensity, & 
                      HeatCap_A, HeatCap_B, HeatCon_A, HeatCon_B, HeatCon_C,  &
                      K_dens, K_ice, T_ptr, K_ptr
  LOGICAL :: GotIt
  TYPE(ValueList_t), POINTER :: Material
    
  
  ! Read variables
    
  Temperature = InputArray(1)  
  RelativeDensity = InputArray(2)
  
  ! Read constants
  
  Material => GetMaterial()
  IceDensity = GetConstReal( Material,'Enthalpy Density', GotIt )
  IF (.NOT. GotIt) THEN 
     CALL FATAL('Get_EnthalpyDiffusivity', 'Enthalpy Density not found in material section')
  END IF
  
  HeatCap_A = ListGetConstReal( Model % Constants, 'Enthalpy Heat Capacity A', GotIt )
  IF (.NOT. GotIt) THEN 
     CALL FATAL('Get_EnthalpyDiffusivity', 'Enthalpy Heat Capacity A not found in constant section')
  END IF
  
  HeatCap_B = ListGetConstReal( Model % Constants, 'Enthalpy Heat Capacity B', GotIt )
  IF (.NOT. GotIt) THEN 
     CALL FATAL('Get_EnthalpyDiffusivity', 'Enthalpy Heat Capacity B not found in constant section')
  END IF
  
  HeatCon_A = ListGetConstReal( Model % Constants, 'Enthalpy Heat Conductivity A', GotIt )
  IF (.NOT. GotIt) THEN 
     CALL FATAL('Get_EnthalpyDiffusivity', 'Enthalpy Heat Conductivity A not found in constant section')
  END IF
  
  HeatCon_B = ListGetConstReal( Model % Constants, 'Enthalpy Heat Conductivity B', GotIt )
  IF (.NOT. GotIt) THEN 
     CALL FATAL('Get_EnthalpyDiffusivity', 'Enthalpy Heat Conductivity B not found in constant section')
  END IF
  
  HeatCon_C = ListGetConstReal( Model % Constants, 'Enthalpy Heat Conductivity C', GotIt )
  IF (.NOT. GotIt) THEN 
     CALL FATAL('Get_EnthalpyDiffusivity', 'Enthalpy Heat Conductivity C not found in constant section')
  END IF
  
  T_ptr = ListGetConstReal( Model % Constants, 'T_triple', GotIt )
  IF (.NOT. GotIt) THEN 
     CALL FATAL('Get_EnthalpyDiffusivity', 'T_triple not found in constant section')
  END IF
  
  ! Calculations
  
  HeatCapacity = HeatCap_A*(273.16+Temperature) + HeatCap_B
 
  K_dens = HeatCon_A*(RelativeDensity*IceDensity)**2 + HeatCon_B*RelativeDensity*IceDensity + HeatCon_C
  
  K_ice = 9.828*EXP(-5.7e-3*(273.16+Temperature))
  
  K_ptr = 9.828*EXP(-5.7e-3*T_ptr)
  
  HeatConductivity = K_ice / K_ptr * K_dens
  
  EnthalpyDiffusivity = HeatConductivity / HeatCapacity  
  
  RETURN

END FUNCTION Get_EnthalpyDiffusivity
