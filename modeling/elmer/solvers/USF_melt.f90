!!!
! New version modified August 2016 
! rm and rs in the sif file 
!!!


FUNCTION Ablation (Model, nodenumber, Zs) RESULT(Source)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   REAL (KIND=dp) :: Zs, Source
   INTEGER :: nodenumber
   
   Source = -max(7.0*(1.0-Zs/1000.0),0.0)
END FUNCTION Ablation

FUNCTION EvolveSMparameter (Model, nodenumber, t) RESULT(sm)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   REAL(KIND=dp) :: t, sm 
   INTEGER :: nodenumber
   TYPE(Variable_t), POINTER :: TimeVar
   REAL(KIND=dp) :: sm0, dsm0

   Timevar => VariableGet( Model % Variables,'Time')
   t = TimeVar % Values(1)

   sm0 = 300.0_dp
   dsm0 = 0.0_dp !9.0_dp
   ! increase sm each year by dsm0
   sm = sm0 + dsm0*AINT(t)
END FUNCTION EvolveSMparameter 

FUNCTION EvolveSMparameterPic (Model, nodenumber, t) RESULT(sm)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   REAL(KIND=dp) :: t, sm 
   INTEGER :: nodenumber
   TYPE(Variable_t), POINTER :: TimeVar
   REAL(KIND=dp) :: sm0, dsm0

   Timevar => VariableGet( Model % Variables,'Time')
   t = TimeVar % Values(1)

   sm0 = 300.0_dp
   dsm0 = 0.0 !37.5_dp 
   ! increase sm each year by dsm0 every 5 years
   sm = sm0
   IF ((t.GE.5.0).AND.(t.LT.6.0)) sm = sm0 + dsm0  
   IF ((t.GE.10.0).AND.(t.LT.11.0)) sm = sm0 + dsm0*2  
   IF ((t.GE.15.0).AND.(t.LT.16.0)) sm = sm0 + dsm0*3   
   IF ((t.GE.20.0).AND.(t.LT.21.0)) sm = sm0 + dsm0*4 
   IF ((t.GE.25.0).AND.(t.LT.26.0)) sm = sm0 + dsm0*5  
   IF ((t.GE.30.0).AND.(t.LT.31.0)) sm = sm0 + dsm0*6  
   IF ((t.GE.35.0).AND.(t.LT.36.0)) sm = sm0 + dsm0*7  
   IF ((t.GE.40.0).AND.(t.LT.41.0)) sm = sm0 + dsm0*8  
END FUNCTION EvolveSMparameterPic 

! Compute the yearly ablation from the melt function of Hewitt 2013
! r(zs,t) = max{0 ; (rm + rs*sm)/2*(tanh((t-tspr)/dt)-tanh((t-taut)/dt)) - rs*zs}
FUNCTION AblationMelt (Model, nodenumber, Zs) RESULT(Source)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   REAL(KIND=dp) :: Zs, Source
   INTEGER :: nodenumber
   TYPE(Variable_t), POINTER :: TimeVar
   TYPE(ValueList_t), POINTER :: Material
   REAL(KIND=dp) :: t, tNext   
   REAL(KIND=dp) :: rm, rs, sm, tspr, taut, dt, rzst, zss
   INTEGER, PARAMETER :: nh = 31
   REAL(KIND=dp) :: abla(nh), ablazs(nh), accu
   REAL(KIND=dp), TARGET, ALLOCATABLE :: abla2zs(:)
   REAL(KIND=dp), DIMENSION(:), POINTER :: abla2zs_p
   INTEGER :: ss, dd 
   LOGICAL :: FirstTime=.TRUE., Found

   SAVE abla, ablazs, abla2zs, FirstTime, tNext

   IF (FirstTime) THEN
      FirstTime = .FALSE.
      ALLOCATE(abla2zs(nh))
      Timevar => VariableGet( Model % Variables,'Time')
      t = TimeVar % Values(1)
      ! we want tNext Smaller than t
      tNext = t - 1.0 
   END IF

   Material => GetMaterial()
   IF (.NOT.ASSOCIATED(Material)) THEN
      WRITE (Message,'(A)') 'No Material found '
      CALL FATAL('AblationMelt',Message)
   END IF
   sm  = GetCReal(Material, 'sm melt Parameter', Found )
   IF (.NOT.Found) THEN
      WRITE(Message,'(A)')'Keyword >sm melt Parameter< not found in Material section'
      CALL FATAL('AblationMelt', Message)
   END IF
   rm  = GetCReal(Material, 'rm melt Parameter', Found )
   IF (.NOT.Found) THEN
      WRITE(Message,'(A)')'Keyword >rm melt Parameter< not found in Material section'
      CALL FATAL('AblationMelt', Message)
   END IF
   rs  = GetCReal(Material, 'rs melt Parameter', Found )
   IF (.NOT.Found) THEN
      WRITE(Message,'(A)')'Keyword >rs melt Parameter< not found in Material section'
      CALL FATAL('AblationMelt', Message)
   END IF
   tspr  = GetCReal(Material, 'tspr melt Parameter', Found )
   IF (.NOT.Found) THEN
      WRITE(Message,'(A)')'Keyword >tspr melt Parameter< not found in Material section'
      CALL FATAL('AblationMelt', Message)
   END IF
   taut  = GetCReal(Material, 'taut melt Parameter', Found )
   IF (.NOT.Found) THEN
      WRITE(Message,'(A)')'Keyword >taut melt Parameter< not found in Material section'
      CALL FATAL('AblationMelt', Message)
   END IF
   dt  = GetCReal(Material, 'dt melt Parameter', Found )
   IF (.NOT.Found) THEN
      WRITE(Message,'(A)')'Keyword >dt melt Parameter< not found in Material section'
      CALL FATAL('AblationMelt', Message)
   END IF

   accu = 0.5 ! m/a

! Real time import
   Timevar => VariableGet( Model % Variables,'Time')
   t = TimeVar % Values(1)

! ablation is constant over a year. 
   IF ( t > tNext ) THEN
      tNext = t + 1.0_dp - Model % Solver % dt 
! Compute the ablation as the annual mean of melt as a function of Zs
! Computed each 100 m and will then be interpolated
! abla(1) for zs = 0.0 and abla(21) for zs = 3000 m
      DO ss = 1,nh
        zss = (ss-1.0)*100.0_dp
        ablazs(ss) = zss 
        abla(ss) = 0.0_dp
! r(zs,1) = r(zs,365) and dx = 1, then        
        DO dd = 2, 364
           rzst = max(0.0, 0.5*(rm + rs*sm)*(tanh((dd-tspr)/dt)-tanh((dd-taut)/dt)) - rs*zss)
           abla(ss) = abla(ss) + rzst 
        END DO
        !mm/a -> m/a
        abla(ss) = abla(ss)*1.0e-3
      END DO
   CALL CubicSpline(nh,ablazs,abla,abla2zs)
   END IF
   
   abla2zs_p => abla2zs
   ! negative because it is here the surface ablation
   ! add the accu assumed not to be zs or t dependant
   ! Source = accu - rhow/rhoi x runoff to be in meter of ice
   Source = accu - 1.099*max(InterpolateCurve(ablazs,abla,Zs,abla2zs_p),0.0)
END FUNCTION AblationMelt


! Compute the basal melt assumed to be the summ of the surface melt + (friction
! + geothermal) 
! r(zs,t) = max{0 ; (rm + rs*sm)/2*(tanh((t-tspr)/dt)-tanh((t-taut)/dt)) - rs*zs}
FUNCTION BasalMelt (Model, nodenumber, Input) RESULT(Source)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   REAL(KIND=dp) :: Input(2)
   REAL(KIND=dp) :: Zs, SSAbasalMelt, Source
   INTEGER :: nodenumber
   TYPE(Variable_t), POINTER :: TimeVar
   TYPE(ValueList_t), POINTER :: Material
   REAL(KIND=dp) :: t, day 
   REAL (KIND=dp) :: rm, rs, sm, tspr, taut, dt, rzst, zss
   INTEGER :: ss, dd 
   LOGICAL :: FirstTime=.TRUE., Found

   Material => GetMaterial()
   IF (.NOT.ASSOCIATED(Material)) THEN
      WRITE (Message,'(A)') 'No Material found '
      CALL FATAL('BasalMelt',Message)
   END IF
   sm  = GetCReal(Material, 'sm melt Parameter', Found )
   IF (.NOT.Found) THEN
      WRITE(Message,'(A)')'Keyword >sm melt Parameter< not found in Material section'
      CALL FATAL('BasalMelt', Message)
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
   dt  = GetCReal(Material, 'dt melt Parameter', Found )
   IF (.NOT.Found) THEN
      WRITE(Message,'(A)')'Keyword >dt melt Parameter< not found in Material section'
      CALL FATAL('BasalMelt', Message)
   END IF


   Zs = Input(1)
   SSAbasalMelt = Input(2)

! Real time import
   Timevar => VariableGet( Model % Variables,'Time')
   t = TimeVar % Values(1)
! time in day of the year
   day = ( t - AINT(t))*365
   
! Values are in mm/d 
! melt in m/a : x 1.e-3 * 365
   Source = 365.0e-3*max(0.0, 0.5*(rm + rs*sm)*(tanh((day-tspr)/dt)-tanh((day-taut)/dt)) - rs*Zs)
   Source = Source + SSAbasalMelt

END FUNCTION BasalMelt

! Compute the mean ablation from the melt function of Hewitt 2013
! r(zs,t) = max{0 ; (rm + rs*sm)/2*(tanh((t-tspr)/dt)-tanh((t-taut)/dt)) - rs*zs}
! It is converted in mean basal melt (m/a) + friction + geothermal heat flux
FUNCTION MeanBasalMelt (Model, nodenumber, Input) RESULT(Source)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   REAL(KIND=dp) :: Input(2)
   REAL(KIND=dp) :: Zs, SSABasalMelt, Source
   INTEGER :: nodenumber
   TYPE(Variable_t), POINTER :: TimeVar
   TYPE(ValueList_t), POINTER :: Material
   REAL(KIND=dp) :: t, tNext  
   REAL(KIND=dp) :: rm, rs, sm, tspr, taut, dt, rzst, zss
   INTEGER, PARAMETER :: nh = 31
   REAL(KIND=dp) :: abla(nh), ablazs(nh)
   REAL(KIND=dp), TARGET, ALLOCATABLE :: abla2zs(:)
   REAL(KIND=dp), DIMENSION(:), POINTER :: abla2zs_p
   INTEGER :: ss, dd 
   LOGICAL :: FirstTime=.TRUE., Found

   SAVE abla, ablazs, abla2zs, FirstTime, tNext

   IF (FirstTime) THEN
      FirstTime = .FALSE.
      ALLOCATE(abla2zs(nh))
      Timevar => VariableGet( Model % Variables,'Time')
      t = TimeVar % Values(1)
      ! we want tNext Smaller than t
      tNext = t - 1.0 
   END IF

   Zs = Input(1)
   SSABasalMelt = Input(2)

   Material => GetMaterial()
   IF (.NOT.ASSOCIATED(Material)) THEN
      WRITE (Message,'(A)') 'No Material found '
      CALL FATAL('MeanBasalMelt',Message)
   END IF
   sm  = GetCReal(Material, 'sm melt Parameter', Found )
   IF (.NOT.Found) THEN
      WRITE(Message,'(A)')'Keyword >sm melt Parameter< not found in Material section'
      CALL FATAL('MeanBasalMelt', Message)
   END IF
   rm  = GetCReal(Material, 'rm melt Parameter', Found )
   IF (.NOT.Found) THEN
      WRITE(Message,'(A)')'Keyword >rm melt Parameter< not found in Material section'
      CALL FATAL('MeanBasalMelt', Message)
   END IF
   rs  = GetCReal(Material, 'rs melt Parameter', Found )
   IF (.NOT.Found) THEN
      WRITE(Message,'(A)')'Keyword >rs melt Parameter< not found in Material section'
      CALL FATAL('MeanBasalMelt', Message)
   END IF
   tspr  = GetCReal(Material, 'tspr melt Parameter', Found )
   IF (.NOT.Found) THEN
      WRITE(Message,'(A)')'Keyword >tspr melt Parameter< not found in Material section'
      CALL FATAL('MeanBasalMelt', Message)
   END IF
   taut  = GetCReal(Material, 'taut melt Parameter', Found )
   IF (.NOT.Found) THEN
      WRITE(Message,'(A)')'Keyword >taut melt Parameter< not found in Material section'
      CALL FATAL('MeanBasalMelt', Message)
   END IF
   dt  = GetCReal(Material, 'dt melt Parameter', Found )
   IF (.NOT.Found) THEN
      WRITE(Message,'(A)')'Keyword >dt melt Parameter< not found in Material section'
      CALL FATAL('MeanBasalMelt', Message)
   END IF

! Real time import
   Timevar => VariableGet( Model % Variables,'Time')
   t = TimeVar % Values(1)

! ablation is constant over a year. 
   IF ( t > tNext ) THEN
      tNext = t + 1.0_dp - Model % Solver % dt 
! Compute the ablation as the annual mean of melt as a function of Zs
! Computed each 100 m and will then be interpolated
! abla(1) for zs = 0.0 and abla(21) for zs = 3000 m
      DO ss = 1,nh
        zss = (ss-1.0)*100.0_dp
        ablazs(ss) = zss 
        abla(ss) = 0.0_dp
! r(zs,1) = r(zs,365) and dx = 1, then        
        DO dd = 2, 364
           rzst = max(0.0, 0.5*(rm + rs*sm)*(tanh((dd-tspr)/dt)-tanh((dd-taut)/dt)) - rs*zss)
           abla(ss) = abla(ss) + rzst 
        END DO
        !mm/a -> m/a
        abla(ss) = abla(ss)*1.0e-3
        write(*,*)zss,abla(ss)
      END DO
   CALL CubicSpline(nh,ablazs,abla,abla2zs)
   END IF
   
   abla2zs_p => abla2zs
   Source = max(InterpolateCurve(ablazs,abla,Zs,abla2zs_p),0.0)
   Source = Source + SSABasalMelt 
END FUNCTION MeanBasalMelt


