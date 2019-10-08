FUNCTION Fluidity_B( model, n, InputArray) RESULT(FluidityParameter)

  USE DefUtils

  IMPLICIT NONE
  
  !----------------------------------------------------------------------------
  ! External variables
  !----------------------------------------------------------------------------
  TYPE(Model_t) :: model         
  INTEGER :: n, i, s
  REAL(KIND=dp) :: InputArray(1)
  REAL(KIND=dp) :: FluidityParameter, Temperature
  REAL(KIND=dp) :: T_pat(12), A_pat(12), M(11), C(11)
  

  Temperature = InputArray(1)
 
  T_pat = [0,-2,-5,-10,-15,-20,-25,-30,-35,-40,-45,-50]
  A_pat = [2.4e-24,1.7e-24,9.3e-25,3.5e-25,2.1e-25,1.2e-25,6.8e-26,3.7e-26,2.0e-26,1.0e-26,5.2e-27,2.6e-27]
  
  M = [3.5e-25,2.567e-25,1.16e-25,2.8e-26,1.8e-26,1.04e-26,6.2e-27,3.4e-27,2e-27,9.6e-28,5.2e-28]
  C = [2.4e-24,2.213e-24,1.51e-24,6.3e-25,4.8e-25,3.28e-25,2.23e-25,1.39e-25,9e-26,4.84e-26,2.86e-26]

  s = 0
  DO i=1,SIZE(T_pat)-1
    IF( Temperature <= T_pat(i) .AND. Temperature > T_pat(i+1) ) THEN
      s = i
    ELSE IF( Temperature == -50) THEN 
      s = 11    
    END IF
  END DO
  
  FluidityParameter = M(s)*Temperature + c(s)

  IF( Temperature > 0) FluidityParameter = A_pat(1)
  IF( Temperature < -50) FluidityParameter = A_pat(12)
  
  FluidityParameter = 2*FluidityParameter*3600*24*365.25*1.0e18

  !PRINT *, Temperature, FluidityParameter
   
  RETURN

END FUNCTION 

 
