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
! *  Modified August 2016 to account for new value for rm and rs
! * 
! *****************************************************************************
!> Compute the Melt for each Moulin by integration over its basin of the melt  
SUBROUTINE ComputeQm( Model,Solver,dt,TransientSimulation )
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
  TYPE(Variable_t), POINTER :: PointerToVariable, ZsSol, MaskSol, & 
                               EiiSol, AMSol, TimeVar
  TYPE(Nodes_t), SAVE :: Nodes

  INTEGER :: i, j, p, n, t, Nn, Nm, istat, EiiDOFs
  INTEGER, POINTER :: Permutation(:), ZsPerm(:), MaskPerm(:), &
                      EiiPerm(:)
  REAL(KIND=dp) :: rm, rs, sm, tspr, taut, ddt, year, day, melt

  REAL(KIND=dp), POINTER :: VariableValues(:), ZsVal(:), Mask(:), &
                            Eii(:)
  REAL(KIND=dp), ALLOCATABLE :: rzst(:), Zs(:), Basis(:), dBasisdx(:,:) 

  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName = 'ComputeQm', ZsName, &
           StrainRateVariableName
  LOGICAL :: AllocationsDone = .FALSE., Found=.FALSE., FirstTime=.TRUE.

  REAL(KIND=dp) :: detJ
  LOGICAL :: Stat
  TYPE(GaussIntegrationPoints_t) :: IP
       
  SAVE SolverName, AllocationsDone, Zs, rzst, Basis, dBasisdx, &
       FirstTime, StrainRateVariableName, ZsName
  !------------------------------------------------------------------------------

  PointerToVariable => Solver % Variable
  Permutation  => PointerToVariable % Perm
  VariableValues => PointerToVariable % Values

  CALL INFO(SolverName, 'Computing Qm for each Moulin', level=3)

  !--------------------------------------------------------------
  !Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------
  IF ( (.NOT. AllocationsDone) .OR. Solver % Mesh % Changed  ) THEN
     n = Solver % Mesh % MaxElementNodes ! just big enough for elemental arrays
     IF (AllocationsDone) DEALLOCATE(rzst, Zs, Basis, dBasisdx)
     ALLOCATE( rzst(n), Zs(n), Basis(n), dBasisdx(n,3), STAT=istat )
     IF ( istat /= 0 ) THEN
        CALL Fatal( SolverName, 'Memory allocation error.' )
     END IF
     AllocationsDone = .TRUE.
     CALL INFO( SolverName, 'Memory allocation done.',Level=1 )
  END IF
  
  IF (FirstTime) THEN
     FirstTime = .FALSE.
     Constants => GetConstants()
     StrainRateVariableName = GetString( Constants, 'StrainRate Variable Name', Found )
     IF (.NOT.Found) StrainRateVariableName = 'StrainRate'

     ZsName = GetString( Constants, 'Zs Variable Name', Found )
     IF (.NOT.Found) ZsName = 'Zs'
  END IF

  !--------------------------------------------------------------
  ! Read the MaskMoulin variable 
  ! This variable contains for each node the node number of the 
  ! moulin it belongs to the basin...
  ! !!! It is real value and we will use it as node index...
  !--------------------------------------------------------------
  MaskSol => VariableGet( Solver % Mesh % Variables, 'MaskMoulin', UnfoundFatal = .TRUE. )
  Mask => MaskSol % Values
  MaskPerm => MaskSol % Perm

  EiiSol => VariableGet( Solver % Mesh % Variables, StrainRateVariableName, UnfoundFatal = .TRUE. )
  EiiPerm => EiiSol % Perm
  EiiDOFs =  EiiSol % DOFs
  Eii => EiiSol % Values

  ZsSol => VariableGet( Solver % Mesh % Variables, ZsName, UnfoundFatal = .TRUE. )
  ZsVal => ZsSol % Values
  ZsPerm => ZsSol % Perm

  !--------------------------------------------------------------
  ! Compute for each element the integrated melt 
  !--------------------------------------------------------------

! Real time import and day of the year
  Timevar => VariableGet( Model % Variables,'Time')
  year = TimeVar % Values(1) 
  day = ( year - AINT(year))*365.0_dp

  VariableValues = 0.0_dp
  DO t = 1, Solver % NumberOfActiveElements
     Element => GetActiveElement(t)
     n = GetElementNOFNodes()
     CALL GetElementNodes( Nodes )
     
     IF (t==1) THEN
        Material => GetMaterial()
        IF (.NOT.ASSOCIATED(Material)) THEN
           WRITE (Message,'(A)') 'No Material found '
           CALL FATAL(SolverName,Message)
        END IF
        rm  = GetCReal(Material, 'rm melt Parameter', Found )
        IF (.NOT.Found) THEN
           WRITE(Message,'(A)')'Keyword >rm melt Parameter< not found in Material section'
           CALL FATAL('BasalMelt', Message)
        END IF
        rs  = GetCReal(Material, 'rs melt Parameter', Found )
        IF (.NOT.Found) THEN
           WRITE(Message,'(A)')'Keyword >rs melt Parameter< not found in Material section'
           CALL FATAL('BasalMelt', Message)
        END IF
        tspr  = GetCReal(Material, 'tspr melt Parameter', Found )
        IF (.NOT.Found) THEN
           WRITE(Message,'(A)')'Keyword >tspr melt Parameter< not found in Material section'
           CALL FATAL('BasalMelt', Message)
        END IF
        taut  = GetCReal(Material, 'taut melt Parameter', Found )
        IF (.NOT.Found) THEN
           WRITE(Message,'(A)')'Keyword >taut melt Parameter< not found in Material section'
           CALL FATAL('BasalMelt', Message)
        END IF
        ddt  = GetCReal(Material, 'ddt melt Parameter', Found )
        IF (.NOT.Found) THEN
           WRITE(Message,'(A)')'Keyword >ddt melt Parameter< not found in Material section'
           CALL FATAL('BasalMelt', Message)
        END IF
        sm  = GetCReal(Material, 'sm melt Parameter', Found )
        IF (.NOT.Found) THEN
           WRITE(Message,'(A)')'Keyword >sm melt Parameter< not found in Material section'
           CALL FATAL(SolverName, Message)
        END IF
        WRITE(*,*)'year, day, sm', year, day, sm
     END IF

     Zs(1:n) = ZsVal(ZsPerm(Element % NodeIndexes(1:n)))
     
     ! Values are in mm/d 
     ! melt in m/a : x 1.e-3 * 365
     ! melt is m of water  
     DO i = 1, n
        rzst(i) = 365.0e-3*max(0.0, 0.5*(rm + rs*sm)*(tanh((day-tspr)/ddt)-tanh((day-taut)/ddt)) - rs*Zs(i))
     END DO

     IP = GaussPoints( Element )
     melt = 0.0_dp
     DO p = 1,IP % n
        stat = ElementInfo( Element, Nodes, IP % U(p), IP % V(p), &
                   IP % W(p),  detJ, Basis, dBasisdx )
        melt = SUM(rzst(1:n)*Basis(1:n)) * detJ * IP % S(p)      

        DO i = 1, n
           ! Return the Moulin number for this node
           Nn = Permutation(NINT(Mask(MaskPerm(Element % NodeIndexes(i)))))
           IF (Nn==0) CYCLE
           VariableValues(Nn) = VariableValues(Nn) + melt * Basis(i) 
        END DO
     END DO 
  END DO
  WRITE (Message,'(A,D14.8)') 'Sum Qm: ', SUM(VariableValues)
  CALL INFO( SolverName, Message, Level=3 )
  
  IF ( ParEnv % PEs>1 ) CALL ParallelSumVector( Solver % Matrix, VariableValues, 1 )
 
  CALL INFO( SolverName , 'Done')
END SUBROUTINE ComputeQm 