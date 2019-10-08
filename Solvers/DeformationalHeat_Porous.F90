!/*****************************************************************************/
! *
! *  Elmer/Ice, a glaciological add-on to Elmer
! *  http://elmerice.elmerfem.org
! *
! * 
! *  This program is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU General Public License
! *  as published by the Free Software Foundation; either version 2
! *  of the License, or (at your option) any later version.
! * 
! *  This program is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! *  GNU General Public License for more details.
! *
! *  You should have received a copy of the GNU General Public License
! *  along with this program (in file fem/GPL-2); if not, write to the 
! *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
! *  Boston, MA 02110-1301, USA.
! *
! *****************************************************************************/
! ******************************************************************************
! *
! *  Authors: 
! *  Email:   
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 
! * 
! *****************************************************************************
!> DOXYGEN INFORMATION TO BE ADDED
SUBROUTINE DeformationalHeatSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
!
!
!  ARGUMENTS:
!
!  TYPE(Model_t) :: Model,  
!     INPUT: All model information (mesh, materials, BCs, etc...)
!
!  TYPE(Solver_t) :: Solver
!     INPUT: Linear & nonlinear equation solver options
!
!  REAL(KIND=dp) :: dt,
!     INPUT: Timestep size for time dependent simulations
!
!  LOGICAL :: TransientSimulation
!     INPUT: Steady state or transient simulation
!
!******************************************************************************
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Element_t),POINTER :: Element
  TYPE(ValueList_t), POINTER :: SolverParams, Material, Constants
  TYPE(Variable_t), POINTER :: PointerToVariable, RelDensitySol, StrainRateSol, DeviatoricStressSol
  TYPE(Solver_t), POINTER :: PointerToSolver

  LOGICAL :: AllocationsDone = .FALSE., Found, UnFoundFatal=.TRUE.

  INTEGER :: i, j,n, m, t, istat,k
  INTEGER, POINTER :: Permutation(:), RelDensityPerm(:), StrainRatePerm(:), DeviatoricStressPerm(:), NodeIndexes(:)

  REAL(KIND=dp), POINTER :: VariableValues(:), RelDensitySolution(:), StrainRateSolution(:), DeviatoricStressSolution(:) 
  REAL(KIND=dp) :: Norm
  Integer :: STDOFs, RD_NSDOFs, SR_NSDOFs, DS_NSDOFs,dim

  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), LOAD(:), FORCE(:), IceDensity(:), RelDensity(:), &
           StrainRate(:,:), DeviatoricStress(:,:)
  
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName

  SAVE STIFF, LOAD, FORCE,  AllocationsDone, IceDensity, RelDensity, StrainRate, DeviatoricStress
!------------------------------------------------------------------------------
  SolverName = 'Deformational Heat Solver'
  PointerToVariable => Solver % Variable
  Permutation  => PointerToVariable % Perm
  VariableValues => PointerToVariable % Values
  STDOFs = PointerToVariable % DOFs

  dim = CoordinateSystemDimension()

  RelDensitySol => VariableGet( Solver % Mesh % Variables, 'Relative Density', UnFoundFatal=UnFoundFatal)
  RelDensityPerm => RelDensitySol % Perm
  RD_NSDOFs = RelDensitySol % DOFs
  RelDensitySolution => RelDensitySol % Values
  
  StrainRateSol => VariableGet( Solver % Mesh % Variables, 'StrainRate', UnFoundFatal=UnFoundFatal)
  StrainRatePerm => StrainRateSol % Perm
  SR_NSDOFs = StrainRateSol % DOFs
  StrainRateSolution => StrainRateSol % Values

  DeviatoricStressSol => VariableGet( Solver % Mesh % Variables, 'DeviatoricStress', UnFoundFatal=UnFoundFatal)
  DeviatoricStressPerm => DeviatoricStressSol % Perm
  DS_NSDOFs = DeviatoricStressSol % DOFs
  DeviatoricStressSolution => DeviatoricStressSol % Values
  
  !Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------
  IF ( (.NOT. AllocationsDone) .OR. Solver % Mesh % Changed  ) THEN
     N = Solver % Mesh % MaxElementNodes ! just big enough for elemental arrays
     M = Model % Mesh % NumberOfNodes
     IF (AllocationsDone) DEALLOCATE(FORCE, LOAD, STIFF, IceDensity, RelDensity, StrainRate, DeviatoricStress)

     ALLOCATE( FORCE(2*STDOFs*N), LOAD(2*STDOFs*N), STIFF(2*STDOFs*N,2*STDOFs*N), IceDensity(N), RelDensity(N), &
          StrainRate(6,N), DeviatoricStress(6,N), STAT=istat )
     IF ( istat /= 0 ) THEN
        CALL Fatal( 'HessianSolve', 'Memory allocation error.' )
     END IF     

     AllocationsDone = .TRUE.
  END IF

  !Initialize the system and do the assembly:
  !------------------------------------------
  CALL DefaultInitialize()

  ! bulk assembly
  DO t=1,Solver % NumberOfActiveElements
     Element => GetActiveElement(t)
     n = GetElementNOFNodes()
     NodeIndexes => Element % NodeIndexes
    
     Material => GetMaterial()
     IceDensity(1:n) = GetReal( Material,'Enthalpy Density', Found )
     IF (.NOT.Found) CALL FATAL(SolverName,'Could not find  >Enthalpy Density<')
    
     RelDensity = 0.0d0
     Do i=1,n
        j = RD_NSDOFs*StrainRatePerm(NodeIndexes(i))
        RelDensity(i) =  RelDensitySolution(j)
     End do
     
     ! This possibly only works for dim=3, please check! (Carlo)
     StrainRate = 0.0d0
     Do i=1,n
        j = SR_NSDOFs*StrainRatePerm(NodeIndexes(i))
        Do k=1,dim*2
          StrainRate(k,i) =  StrainRateSolution( j-(dim*2)+k )
        End do
     End do
     
     ! This possibly only works for dim=3, please check! (Carlo)
     DeviatoricStress = 0.0d0
     Do i=1,n
        j = DS_NSDOFs*DeviatoricStressPerm(NodeIndexes(i))
        Do k=1,dim*2
          DeviatoricStress(k,i) =  DeviatoricStressSolution( j-(dim*2)+k )
        End do
     End do
     
     CALL LocalMatrix(  STIFF, FORCE, Element, n, IceDensity, RelDensity, StrainRate, DeviatoricStress )
     CALL DefaultUpdateEquations( STIFF, FORCE )
  END DO
  

  CALL DefaultFinishAssembly()

  Norm = DefaultSolve()

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  STIFF, FORCE, Element,  n, IceDensity, RelDensity, StrainRate, DeviatoricStress )
!------------------------------------------------------------------------------
    USE MaterialModels

    REAL(KIND=dp) :: STIFF(:,:), FORCE(:), IceDensity(:), RelDensity(:), StrainRate(:,:), DeviatoricStress(:,:)
    INTEGER :: n, dim
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3), DetJ
    REAL(KIND=dp) :: IceDens, DeformationalHeat, t1, t2, t3, t4, t5, t6                    
    LOGICAL :: Stat
    INTEGER :: t, p, q, i
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    FORCE = 0.0d0
        
    dim = CoordinateSystemDimension()

    IP = GaussPoints( Element )
    DO t=1,IP % n
       
       stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t),  detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )
       
       IceDens = SUM( IceDensity(1:n) * Basis(1:n) )

       ! This only works for dim=3 !! (Carlo)
       t1 = SUM( StrainRate(1,1:n) * DeviatoricStress(1,1:n) / RelDensity(1:n) * Basis(1:n) )
       t2 = SUM( StrainRate(2,1:n) * DeviatoricStress(2,1:n) / RelDensity(1:n) * Basis(1:n) )
       t3 = SUM( StrainRate(3,1:n) * DeviatoricStress(3,1:n) / RelDensity(1:n) * Basis(1:n) ) 
       t4 = SUM( 2 * StrainRate(4,1:n) * DeviatoricStress(4,1:n) / RelDensity(1:n) * Basis(1:n) )
       t5 = SUM( 2 * StrainRate(5,1:n) * DeviatoricStress(5,1:n) / RelDensity(1:n) * Basis(1:n) )
       t6 = SUM( 2 * StrainRate(6,1:n) * DeviatoricStress(6,1:n) / RelDensity(1:n) * Basis(1:n) )       
       
       DeformationalHeat = (t1 + t2 + t3 + t4 + t5 + t6) / IceDens

       DO p=1,n
          DO q=1,n
              STIFF(p,q) = STIFF(p,q) + IP % S(t) * detJ * Basis(p) * Basis(q)
          END DO
       END DO
       DO p=1,n
          Force(p) = Force(p) + IP % S(t) * detJ * DeformationalHeat * Basis(p)
       END DO
    END DO       
    
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
END SUBROUTINE DeformationalHeatSolver
!------------------------------------------------------------------------------

