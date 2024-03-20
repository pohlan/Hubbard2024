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
! *  Authors: Olivier Gagliardini
! *  Email:   
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 10 April 2015 
! * 
! *****************************************************************************
!> Solver for creating a mask to attribute the basin connected to each moulin 
!> The mask number is the node number to which the Moulin belongs
!> A node is part of the basin of its closest moulin
SUBROUTINE MaskBasin( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
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
  TYPE(ValueList_t), POINTER :: Material, SolverParams, Constants
  TYPE(Variable_t), POINTER :: PointerToVariable, ZsSol, AMSol, EiiSol, TimeVar
  TYPE(Nodes_t), SAVE :: Nodes

  INTEGER :: i, j, n, t, Nn, Nm, EiiDOFs, ClosestM
  INTEGER, POINTER :: Permutation(:), ZsPerm(:), ActiveMoulinPerm(:), EiiPerm(:)

  REAL(KIND=dp), POINTER :: VariableValues(:), ZsVal(:), ActiveMoulin(:), Eii(:)
  REAL(KIND=dp), ALLOCATABLE :: xm(:), ym(:)
  INTEGER, ALLOCATABLE :: MoulinNumber(:) 
  REAL(KIND=dp) :: x, y, MinD, dist, Eth 

  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName = 'MaskBasin', &
        StrainRateVariableName, ZsName
  LOGICAL :: FirstTime = .True., Found = .False.  
       
  SAVE SolverName, Nm, MoulinNumber, xm, ym, StrainRateVariableName, ZsName, Eth
  !------------------------------------------------------------------------------

  PointerToVariable => Solver % Variable
  Permutation  => PointerToVariable % Perm
  VariableValues => PointerToVariable % Values

  CALL INFO(SolverName, 'Computing MaskMoulin from moulin position', level=3)

  IF (FirstTime) THEN
     FirstTime = .False. 

     Constants => GetConstants()
     StrainRateVariableName = GetString( Constants, 'StrainRate Variable Name', Found )
     IF (.NOT.Found) StrainRateVariableName = 'StrainRate'

     ZsName = GetString( Constants, 'Zs Variable Name', Found )
     IF (.NOT.Found) ZsName = 'Zs'

  !--------------------------------------------------------------
  ! Determine the number of Moulin = number of 101 Boundary elts
  !--------------------------------------------------------------
     Nm = 0
     DO t=1,Solver % Mesh % NumberOfBoundaryElements
        Element => GetBoundaryElement(t)
        IF ( GetElementFamily() == 1 ) THEN 
           Nm = Nm + 1
        END IF
     END DO

  !--------------------------------------------------------------
  ! Get the Moulin coordinates
  !--------------------------------------------------------------
     ALLOCATE(xm(Nm), ym(Nm), MoulinNumber(Nm))
     Nm = 0
     DO t=1,Solver % Mesh % NumberOfBoundaryElements
        Element => GetBoundaryElement(t)
        IF ( GetElementFamily() == 1 ) THEN 
           Nm = Nm + 1
           CALL GetElementNodes( Nodes )
           xm(Nm) = Nodes % x(1) 
           ym(Nm) = Nodes % y(1) 
           MoulinNumber(Nm) = Element % NodeIndexes(1)
        END IF
     END DO
     write(*,*)'found ',Nm,' moulins'

  END IF ! End first Time

  !--------------------------------------------------------------
  ! Read the Zs for each elements (added 04/05/2017)
  !--------------------------------------------------------------

  ZsSol => VariableGet( Solver % Mesh % Variables, ZsName, UnfoundFatal = .TRUE. )
  ZsVal => ZsSol % Values
  ZsPerm => ZsSol % Perm

  ! This variable is < 0 if a moulin is inactive (initial condition)
  ! A moulin is activated (> 0) if Exx > E_th (=0.05) 
  AMSol => VariableGet( Solver % Mesh % Variables, 'Active Moulin', UnfoundFatal = .TRUE. )
  ActiveMoulin => AMSol % Values
  ActiveMoulinPerm => AMSol % Perm


  ! Need strain rate to activate moulins
  EiiSol => VariableGet( Solver % Mesh % Variables, StrainRateVariableName, UnfoundFatal = .TRUE. )
  EiiPerm => EiiSol % Perm
  EiiDOFs =  EiiSol % DOFs
  Eii => EiiSol % Values


  !--------------------------------------------------------------
  ! Loop over all moulins to see which ones are activated       
  ! Save at which date it is activated
  !--------------------------------------------------------------
  Timevar => VariableGet( Model % Variables,'Time')

  DO t = 1, Nm
     j = MoulinNumber(t)
     Element => GetActiveElement(j)
     Material => GetMaterial()
     IF (.NOT.ASSOCIATED(Material)) THEN
        WRITE (Message,'(A)') 'No Material found '
        CALL FATAL(SolverName,Message)
     END IF
     Eth  = GetCReal(Material, 'Crevasse Criteria', Found )
     IF (.NOT.Found) THEN
        WRITE(Message,'(A)')'Keyword >Crevasse Criteria< not found in Material section'
        CALL FATAL(SolverName, Message)
     END IF
     ! only look at unactivated moulins - 
     IF (ActiveMoulin(ActiveMoulinPerm(j)) < 0.0) THEN
        IF ((Eii(EiiDOFs*(EiiPerm(j)-1)+1) > Eth) & 
                       &.OR.(ZsVal(ZsPerm(j))<1000.0)) THEN  
   !!!!!  ATTENTION CHANGED 1600 by 1000 !!!!!!!!  
           ActiveMoulin(ActiveMoulinPerm(j)) = TimeVar % Values(1)
           WRITE (Message,'(A,I6)') 'New Moulin Actived: ', j 
           CALL INFO( SolverName, Message, Level=3 )
        END IF
     END IF
  END DO


  !--------------------------------------------------------------
  ! Loop over all elements to find in which basin they belong       
  !--------------------------------------------------------------
  VariableValues = 0
  DO t = 1, Solver % NumberOfActiveElements
     Element => GetActiveElement(t)
     n = GetElementNOFNodes()
     CALL GetElementNodes( Nodes )
     
     DO i = 1, n
        Nn = Permutation(Element % NodeIndexes(i))
        IF (Nn==0) CYCLE
        ! If we have already visited that node
        ! IF (ABS(VariableValues(Nn)- TimeVar % Values(1)) < 1.0e-10) CYCLE

        !--------------------------------------------------------------
        ! Find the closest ACTIVATED moulin of that node
        ! Take into account the altitude of the node w/ the moulin: Zs(m)>Zs(n)
        ! (added 04/05/2017)
        !--------------------------------------------------------------
        ! Note 14/09 : it was n instead of i, but should not change that much
        ! results????
        x = Nodes % x(i)
        y = Nodes % y(i)

        MinD = 1.0e20     
        ClosestM = 0
        DO j = 1, Nm
           ! Look only at activated moulins at elevation lower than Zs for that node  
           IF (ActiveMoulin(ActiveMoulinPerm(MoulinNumber(j))) < 0.0) CYCLE
           IF (ZsVal(ZsPerm(MoulinNumber(j)))>ZsVal(ZsPerm(Element%NodeIndexes(i)))) CYCLE
           dist = SQRT((x-xm(j))**2+(y-ym(j))**2)  
           IF (dist < MinD) THEN
              MinD = dist
              ClosestM = j
           END IF
        END DO
        IF (ClosestM > 0) THEN ! To account for nodes at lower than lowest moulin 
           VariableValues(Nn) = MoulinNumber(ClosestM)
           ! For nodes (except moulin) store the date this basin become active   
           IF (.Not.ANY(Element % NodeIndexes(i)==MoulinNumber)) THEN
              ActiveMoulin(ActiveMoulinPerm(Element%NodeIndexes(i))) = &
                        ActiveMoulin(ActiveMoulinPerm(MoulinNumber(ClosestM)))
             !write(*,*)'node',Nn, ClosestM, MoulinNumber(ClosestM), ActiveMoulin(ActiveMoulinPerm(MoulinNumber(ClosestM)))
           ELSE
             !write(*,*)'node Moulin',Nn,Element%NodeIndexes(i)
           END IF                
        END IF
     END DO
  END DO
  
  IF ( ParEnv % PEs>1 ) CALL ParallelSumVector( Solver % Matrix, VariableValues, 1 )
 
  CALL INFO( SolverName , 'Done')
END SUBROUTINE MaskBasin
