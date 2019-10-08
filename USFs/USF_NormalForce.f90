FUNCTION CalcNormalForce( model, n, InputArray) RESULT(normal_force)

  USE DefUtils

  IMPLICIT NONE
  
  !----------------------------------------------------------------------------
  ! External variables
  !----------------------------------------------------------------------------
  TYPE(Model_t) :: model         
  INTEGER :: n        
  REAL(KIND=dp) :: InputArray(3)
  REAL(KIND=dp) :: normal_force  
  
  !----------------------------------------------------------------------------
  ! Internal variables
  !----------------------------------------------------------------------------
  REAL(KIND=dp) :: depth, reldens, coord2, Gravity, IceDensity, CrevDepth, StressGrad
  INTEGER :: DIM, other_body_id, nboundary, nparent,&
       BoundaryElementNode, ParentElementNode
  Logical :: GotIt
  TYPE(ValueList_t), POINTER :: BC
  TYPE(Element_t), POINTER :: BoundaryElement, ParentElement
  TYPE(Nodes_t) :: ElementNodes
  
  
  ! -----------------------
  ! Get element information
  ! -----------------------
  BoundaryElement => CurrentModel % CurrentElement
  IF ( .NOT. ASSOCIATED(BoundaryElement) ) THEN
     CALL FATAL('CalcNormalForce','No boundary element found')
  END IF

  BC => GetBC()
  IF ( .NOT. ASSOCIATED(BC) ) THEN
     CALL FATAL('CalcNormalForce','No boundary condition found')
  END IF
  
  other_body_id = BoundaryElement % BoundaryInfo % outbody
  ! only one body in simulation
  IF (other_body_id < 1) THEN
     ParentElement => BoundaryElement % BoundaryInfo % Right
     IF ( .NOT. ASSOCIATED(ParentElement) )&
          ParentElement => BoundaryElement % BoundaryInfo % Left
  ! we are dealing with a body-body boundary
  ! and assume that the normal is pointing outwards
  ELSE
     ParentElement => BoundaryElement % BoundaryInfo % Right
     IF (ParentElement % BodyId == other_body_id)&
          ParentElement => BoundaryElement % BoundaryInfo % Left
  END IF
  
  ! just to be on the save side, check again
  IF ( .NOT. ASSOCIATED(ParentElement) ) THEN
     WRITE(Message, *)&
          'No parent element found for boundary element no. ', n
     CALL FATAL('CalcNormalForce',Message)
  END IF
  
  ! get the corresponding node in the elements
  nboundary = GetElementNOFNodes(BoundaryElement)
  DO BoundaryElementNode=1,nboundary
     IF ( n == BoundaryElement % NodeIndexes(BoundaryElementNode) ) EXIT
  END DO
  nparent = GetElementNOFNodes(ParentElement)
  DO ParentElementNode=1,nboundary
    IF ( n == ParentElement % NodeIndexes(ParentElementNode) ) EXIT
  END DO

  ! get element nodes
  CALL GetElementNodes(ElementNodes, BoundaryElement)  
  
  
  ! ----------------------------------------------------
  ! Get parameters from constants section and BC section
  ! ----------------------------------------------------
  DIM = CoordinateSystemDimension()
  Gravity = ListGetConstReal( Model % Constants, 'Gravity',GotIt)
  IF (.NOT. GotIt) THEN 
     CALL FATAL('CalcNormalForce', 'Gravity not found in constant section')
  END IF
  
  IceDensity = ListGetConstReal( Model % Constants, 'Ice Density',GotIt)
  IF (.NOT. GotIt) THEN 
     CALL FATAL('CalcNormalForce', 'Ice Density not found in constant section')
  END IF
  
  CrevDepth = GetConstReal( BC,'Crevasse Depth',GotIt)
  IF (.NOT. GotIt) THEN
     CALL FATAL('CalcNormalForce', 'Crevasse Depth not found in BC section')
  END IF
  
  StressGrad = GetConstReal( BC,'Stress Gradient',GotIt)
  IF (.NOT. GotIt) THEN
     CALL FATAL('CalcNormalForce', 'Stress Gradient not found in BC section')
  END IF
  
  !---------------------
  ! Compute normal force
  !---------------------    
  depth = InputArray(1)
  reldens = InputArray(2)
  coord2 = InputArray(3)
  
  IF (coord2 < 86270) THEN
    CrevDepth = 0.0
    StressGrad = 0.0    
  END IF
  
  IF (coord2 > 86500) THEN
    CrevDepth = 0.0
  END IF
  
  IF (depth > CrevDepth) THEN
     normal_force = StressGrad*(depth-CrevDepth)
  ELSE
     normal_force = 0.0
  END IF
  
  RETURN

END FUNCTION CalcNormalForce
