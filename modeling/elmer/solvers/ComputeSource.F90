! Compute the basal melt assumed to be the summ of the surface melt + (friction
! + geothermal) 
! r(zs,t) = max{0 ; (rm + rs*sm)/2*(tanh((t-tspr)/dt)-tanh((t-taut)/dt)) - rs*zs}
FUNCTION BasalSource (Model, nodenumber, Input) RESULT(Source)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   REAL(KIND=dp) :: Zs, Input, Source
   INTEGER :: nodenumber
   TYPE(Variable_t), POINTER :: TimeVar
   TYPE(ValueList_t), POINTER :: Material
   REAL(KIND=dp) :: t, day 
   REAL (KIND=dp) :: rm, rs, sm, tspr, taut, dt, ConstantBasalMelt
   INTEGER :: ss, dd 
   LOGICAL :: FirstTime=.TRUE., Found

   Material => GetMaterial()
   IF (.NOT.ASSOCIATED(Material)) THEN
      WRITE (Message,'(A)') 'No Material found '
      CALL FATAL('BasalSource',Message)
   END IF
   sm  = GetCReal(Material, 'sm melt Parameter', Found )
   IF (.NOT.Found) THEN
      WRITE(Message,'(A)')'Keyword >sm melt Parameter< not found in Material section'
      CALL FATAL('BasalSource', Message)
   END IF
   rm  = GetCReal(Material, 'rm melt Parameter', Found )
   IF (.NOT.Found) THEN
      WRITE(Message,'(A)')'Keyword >rm melt Parameter< not found in Material section'
      CALL FATAL('BasalSource', Message)
   END IF
   rs  = GetCReal(Material, 'rs melt Parameter', Found )
   IF (.NOT.Found) THEN
      WRITE(Message,'(A)')'Keyword >rs melt Parameter< not found in Material section'
      CALL FATAL('BasalSource', Message)
   END IF 
   tspr  = GetCReal(Material, 'tspr melt Parameter', Found )
   IF (.NOT.Found) THEN
      WRITE(Message,'(A)')'Keyword >tspr melt Parameter< not found in Material section'
      CALL FATAL('BasalSource', Message)
   END IF
   taut  = GetCReal(Material, 'taut melt Parameter', Found )
   IF (.NOT.Found) THEN
      WRITE(Message,'(A)')'Keyword >taut melt Parameter< not found in Material section'
      CALL FATAL('BasalSource', Message)
   END IF
   dt  = GetCReal(Material, 'dt melt Parameter', Found )
   IF (.NOT.Found) THEN
      WRITE(Message,'(A)')'Keyword >dt melt Parameter< not found in Material section'
      CALL FATAL('BasalSource', Message)
   END IF
   ConstantBasalMelt  = GetCReal(Material, 'Basal Melt Parameter', Found )
   IF (.NOT.Found) THEN
      WRITE(Message,'(A)')'Keyword >Basal melt Parameter< not found in Material section'
      CALL FATAL('BasalSource', Message)
   END IF


   Zs = Input

! Real time import
   Timevar => VariableGet( Model % Variables,'Time')
   t = TimeVar % Values(1)
! time in day of the year
   day = ( t - AINT(t))*365
  !print *, "Time in days", day
   

! should be about 200* ConstantBasalMelt for 3 cm/day 

  Source = ConstantBasalMelt + 125*ConstantBasalMelt * exp(-((day - 183.0)**2) / (20.0**2))




  ! IF (tspr < day .AND. day < taut) THEN
  !    Source = (ConstantBasalMelt * 40)*exp(-1000.0 / (3000.0 - (day - 182.0)**2))
  ! ELSE
  !    Source = ConstantBasalMelt
  ! END IF


  ! print *, "Computed basal source of", Source

END FUNCTION BasalSource
