MODULE VSTSDOF
USE BaseMethod, ONLY: Display, EqualLine, &
                      BlankLines, input, UpperCase, INV, &
                      Reallocate, CHAR_LF, stdout, &
                      TOSTRING, CHAR_SLASH
USE Lapack_Method, ONLY: Solve
USE GlobalData, ONLY: DFP, LGT, I4B, Equidistance, Monomial
USE BaseType, ONLY: ElemShapeData_, QuadraturePoint_, ReferenceLine_, &
                    TypeH1, TypeLIP => TypeLagrangeInterpolation
USE ElemshapeData_Method
USE QuadraturePoint_Method
USE ReferenceLine_Method
USE Utility, ONLY: OUTERPROD
USE TxtFile_Class
USE String_Class
USE GnuPlot_Class
USE UserFunction_Class
USE CSVFile_Class
USE ExceptionHandler_Class, ONLY: e
USE TomlUtility
USE tomlf, ONLY: toml_table, &
  & toml_serialize,  &
  & toml_get => get_value, &
  & toml_len => len, &
  & toml_array,  &
  & toml_stat
IMPLICIT NONE

PRIVATE
PUBLIC :: VSTSDOF_
CHARACTER(*), PARAMETER :: modName = "VSTSDOF_"
CHARACTER(*), PARAMETER :: prefix = "VSTSDOF"
REAL(DFP), PARAMETER :: default_alpha = 1.0_DFP
INTEGER(I4B), PARAMETER :: default_order = 1
CHARACTER(*), PARAMETER :: default_result_dir = "./results"
CHARACTER(*), PARAMETER :: default_filename = "VSTSDOF"
INTEGER(I4B), PARAMETER :: default_verbosity = 0

!----------------------------------------------------------------------------
!
!----------------------------------------------------------------------------

!> author: Shion Shimizu
! date:   2024-05-10
! summary:  Space-Time Algorithm type is defined here

TYPE :: VSTSDOF_
  REAL(dfp), ALLOCATABLE :: ct(:, :), mt(:, :), kt(:, :)
  REAL(DFP), ALLOCATABLE :: tanmat(:, :)
  REAL(DFP), ALLOCATABLE :: rhs(:), sol(:)
  REAL(DFP) :: u0, v0, a0
  REAL(DFP), ALLOCATABLE :: uvar(:), vvar(:), avar(:)

  REAL(DFP), ALLOCATABLE :: rhs_u1(:, :)
  !! size will be (3, obj%nnt)
  !! NOTE: be careful for the order
  !! col1 if for rhs 1 (equation at time tn)
  !! col2 is for rhs 2 (equation at time tn+dt)
  !! col 3 to obj%nnt for rhs 3 to obj%nnt
  !! (equation at time tn < t < tn+dt)
  !!
  !! coefficient for rhs for (displacement/dt) at time tn
  !!
  !! rhs1 += (rhs_u1(1,1)*M + rhs_u1(2,1)*C*dt + rhs_u1(3,1)*K*dt^2)*u/dt
  !! rhs2 += (rhs_u1(1,2)*M + rhs_u1(2,2)*C*dt + rhs_u1(3,2)*K*dt^2)*u/dt
  !!
  !! rhs_u1(1,a) = coefficient for M
  !! rhs_u1(2,a) = coefficient for C*dt
  !! rhs_u1(3,a) = coefficient for K*dt^2

  REAL(DFP), ALLOCATABLE :: rhs_v1(:, :)
  !! size will be (3, obj%nnt)
  !! the order rule is the same as one of rhs_u1
  !!
  !! coefficient for rhs for (velocity) at time tn
  !!
  !! rhs1 += (rhs_v1(1,1)*M + rhs_v1(2,1)*C*dt + rhs_v1(3,1)*K*dt^2)*v
  !! rhs2 += (rhs_v1(1,2)*M + rhs_v1(2,2)*C*dt + rhs_v1(3,2)*K*dt^2)*v
  !!
  !! rhs_v1(1,a) = coefficient for M
  !! rhs_v1(2,a) = coefficient for C*dt
  !! rhs_v1(3,a) = coefficient for K*dt^2

  REAL(DFP), ALLOCATABLE :: disp(:)
  !! size will be (2+obj%nnt)
  !! disp coefficient for displacement update
  !! disp = disp(1)*u0 + disp(2) * v0 * dt + disp(3:) * sol(:) * dt
  !!
  !! disp(1) coefficient of u0
  !! disp(2) coefficient of v0 * dt
  !! disp(3:) coefficients of sol(:) * dt

  REAL(DFP), ALLOCATABLE :: vel(:)
  !! size will be (2+obj%nnt)
  !! vel coefficient for velocity update
  !! vel = vel(1)*u0/dt + vel(2) * v0 + vel(3:) * sol(:)
  !!
  !! vel(1) coefficient of u0 / dt
  !! vel(2) coefficient of v0
  !! vel(3:) coefficient of sol(:)

  INTEGER(I4B) :: totalTimeNodes = 0
  !! total nodes in time domain

  INTEGER(I4B) :: totalTimeElements = 0
  !! total elements in time domain

  REAL(DFP), ALLOCATABLE :: timeElemLength(:)
  !! length of each time element
  !! the size should be totalTimeElements

  INTEGER(I4B), ALLOCATABLE :: timeOrder(:)
  !! time order of each element

  REAL(DFP), ALLOCATABLE :: alpha(:)
  !! parameter for lcvst fem

  INTEGER(I4B) :: maxTimeOrder = 0
  !! maximum time order

  INTEGER(I4B) :: currentTimeStep = 1
  !! current time step

  REAL(DFP) :: currentTime = 0.0_DFP
  !! current time

  REAL(DFP) :: timeDomain(2) = 0.0_DFP
  !! Total length of time domain

  INTEGER(I4B) :: verbosity = 0
  !! verbosity level
  !! 0 means minimum

  REAL(DFP) :: naturalFrequency
  !! natural frequency (omega)

  REAL(DFP) :: dampingRatio
  !! damping ratio (xi)

  TYPE(String) :: result_dir
  !! Result directory name

  TYPE(String) :: filename
  !! Filename

  TYPE(GnuPlot_) :: plot
  !! for plotting

  TYPE(UserFunction_), POINTER :: initialVel => NULL()
  !! initial condition for velocity

  TYPE(UserFunction_), POINTER :: initialDisp => NULL()
  !! initial condition for displacement

  TYPE(CSVFile_) :: dispfile, velfile, accfile, datafile
  !! file to write displacement, velocity, acceleration

  LOGICAL(LGT) :: saveData(4)
  !! boolean to decide write the data of
  !! diaplacement, velocity, acceleration, and all

  LOGICAL(LGT) :: plotData(3)
  !! boolean to decide plot the data of
  !! diaplacement, velocity, acceleration
CONTAINS
  PROCEDURE, PUBLIC, PASS(obj) :: DEALLOCATE => obj_Deallocate
  PROCEDURE, PUBLIC, PASS(obj) :: Display => obj_Display
  PROCEDURE, PUBLIC, PASS(obj) :: VSTFEM => obj_VSTFEM
  PROCEDURE, PUBLIC, PASS(obj) :: SetInitialDisplacement => obj_SetInitialDisplacement
  PROCEDURE, PUBLIC, PASS(obj) :: SetInitialVelocity => obj_SetInitialVelocity
  PROCEDURE, PUBLIC, PASS(obj) :: GetTanmat => obj_GetTanmat
  PROCEDURE, PUBLIC, PASS(obj) :: getrhs => obj_getrhs
  PROCEDURE, PUBLIC, PASS(obj) :: solve => obj_solve
  PROCEDURE, PUBLIC, PASS(obj) :: update => obj_update
  PROCEDURE, PUBLIC, PASS(obj) :: Set => obj_Set
  PROCEDURE, PUBLIC, PASS(obj) :: Run => obj_Run
  PROCEDURE, PUBLIC, PASS(obj) :: writeData => obj_writeData
  PROCEDURE, PASS(obj) :: ImportFromToml1 => obj_ImportFromToml1
  PROCEDURE, PASS(obj) :: ImportFromToml2 => obj_ImportFromToml2
  GENERIC, PUBLIC :: ImportFromToml => ImportFromToml1, ImportFromToml2
END TYPE VSTSDOF_

CONTAINS

!----------------------------------------------------------------------------
!                                                               Deallocate
!----------------------------------------------------------------------------

SUBROUTINE obj_Deallocate(obj)
  CLASS(VSTSDOF_), INTENT(INOUT) :: obj

  IF (ALLOCATED(obj%tanmat)) DEALLOCATE (obj%tanmat)

  IF (ALLOCATED(obj%rhs_u1)) DEALLOCATE (obj%rhs_u1)

  IF (ALLOCATED(obj%rhs_v1)) DEALLOCATE (obj%rhs_v1)

  IF (ALLOCATED(obj%disp)) DEALLOCATE (obj%disp)

  IF (ALLOCATED(obj%vel)) DEALLOCATE (obj%vel)

END SUBROUTINE obj_Deallocate

!----------------------------------------------------------------------------
!                                                                  Display
!----------------------------------------------------------------------------

SUBROUTINE obj_Display(obj, msg, unitno)
  CLASS(VSTSDOF_), INTENT(in) :: obj
  CHARACTER(*), INTENT(in) :: msg
  INTEGER(I4B), OPTIONAL, INTENT(IN) :: unitno
#ifdef DEBUG_VER
  CHARACTER(*), PARAMETER :: myName = "obj_Display()"
#endif

  LOGICAL(LGT) :: isok

#ifdef DEBUG_VER
  CALL e%RaiseInformation(modName//'::'//myName//' - '// &
                          '[START] ')
#endif

  CALL Display(msg)
  CALL Display(obj%totalTimeNodes, "totalTimeNodes: ", unitno=unitno)

  CALL Display(obj%maxTimeOrder, "maxTimeOrder: ", unitno=unitno)
  CALL Display(obj%timeDomain, "timeDomain: ", unitno=unitno)

  isok = ALLOCATED(obj%timeOrder)
  CALL Display(isok, "timeOrder: ALLOCATED: ", unitno=unitno)
  IF (isok) THEN
    CALL Display(obj%timeOrder, "timeOrder: ", unitno=unitno)
  END IF

  isok = ALLOCATED(obj%timeElemLength)
  CALL Display(isok, "timeElemLength: ALLOCATED: ", unitno=unitno)
  IF (isok) THEN
    CALL Display(obj%timeElemLength, "timeElemLength: ", unitno=unitno)
  END IF

  CALL Display(obj%result_dir%chars(), "result_dir: ", unitno=unitno)
  CALL Display(obj%filename%chars(), "filename: ", unitno=unitno)

  isok = ASSOCIATED(obj%initialDisp)
  CALL Display(isok, "initialDisp ASSOCIATED: ", unitno=unitno)

  isok = ASSOCIATED(obj%initialVel)
  CALL Display(isok, "initialVel ASSOCIATED: ", unitno=unitno)

  CALL Display("Output Controll")
  CALL Display(obj%saveData, "saveData: ", unitno=unitno)
  CALL Display(obj%plotData, "plotData: ", unitno=unitno)

#ifdef DEBUG_VER
  CALL e%RaiseInformation(modName//'::'//myName//' - '// &
                          '[END] ')
#endif

END SUBROUTINE obj_Display

!----------------------------------------------------------------------------
!                                                             ImportFromToml
!----------------------------------------------------------------------------

!> author: Shion Shimizu
! date:  2024-05-11
! summary:  Initiate all paramters of VSTAlgo from given toml table

SUBROUTINE obj_ImportFromToml1(obj, table)
  CLASS(VSTSDOF_), INTENT(inout) :: obj
  TYPE(toml_table), INTENT(inout) :: table

  CHARACTER(*), PARAMETER :: myName = "obj_ImportFromToml1()"
  INTEGER(I4B) :: origin, stat, tsize
  LOGICAL(LGT) :: isok, abool
  TYPE(String) :: astr
  INTEGER(I4B), ALLOCATABLE :: tempintvec(:)
  REAL(DFP), ALLOCATABLE :: temprealvec(:)
  TYPE(toml_table), POINTER :: node
  TYPE(toml_array), POINTER :: array
  LOGICAL(LGT), ALLOCATABLE :: tempboolvec(:)

#ifdef DEBUG_VER
  CALL e%RaiseInformation(modName//'::'//myName//' - '// &
                          '[START]')
#endif

  CALL obj%DEALLOCATE()

#ifdef DEBUG_VER
  CALL Display(myName//" result_dir")
#endif

  CALL GetValue(table=table, key="result_dir", VALUE=obj%result_dir, &
     default_value=default_result_dir, origin=origin, stat=stat, isfound=isok)

#ifdef DEBUG_VER
  CALL Display(myName//" filename")
#endif

  CALL GetValue(table=table, key="filename", VALUE=obj%filename, &
       origin=origin, stat=stat, isfound=isok, default_value=default_filename)

  CALL GetValue(table=table, key="verbosity", VALUE=obj%verbosity, &
      origin=origin, stat=stat, isfound=isok, default_value=default_verbosity)

#ifdef DEBUG_VER
  CALL Display(myName//" totalTimeNodes")
#endif

  CALL GetValue(table=table, key="totalTimeNodes", VALUE=obj%totalTimeNodes, &
                default_value=0_I4B, origin=origin, stat=stat, isfound=isok)

#ifdef DEBUG_VER
  CALL Display(myName//" totalTimeElements")
#endif

  CALL GetValue(table=table, key="totalTimeElements", &
                VALUE=obj%totalTimeElements, &
                default_value=0_I4B, origin=origin, stat=stat, isfound=isok)

#ifdef DEBUG_VER
  CALL Display(myName//" timeDomain")
#endif

  CALL GetValue_(table=table, key="timeDomain", tsize=tsize, &
                 VALUE=obj%timeDomain, origin=origin, stat=stat, isfound=isok)
  isok = tsize .EQ. 2
  CALL AssertError1(isok, myname, "timeDomain should have 2 values")

!INFO: timeOrder
#ifdef DEBUG_VER
  CALL Display(myName//" timeOrder")
#endif

  CALL GetValue(table=table, key="timeOrder", VALUE=tempintvec, &
                origin=origin, stat=stat, isfound=isok)
  CALL AssertError1(isok, myname, "timeOrder not found")

#ifdef DEBUG_VER
  CALL Display(obj%totalTimeElements, myname//" totalTimeElements: ")
#endif

  CALL Reallocate(obj%timeOrder, obj%totalTimeElements)

  abool = SIZE(tempintvec) .EQ. 1
  IF (abool) THEN
    obj%timeOrder(:) = tempintvec(1)
  ELSE
    isok = SIZE(tempintvec) .EQ. obj%totalTimeElements
    CALL AssertError1(isok, myname, "timeOrder should have "// &
                      "totalTimeElements values")
    obj%timeOrder(:) = tempintvec(1:obj%totalTimeElements)
  END IF

  obj%maxTimeOrder = MAXVAL(obj%timeOrder)

! !INFO: timeElemLength
#ifdef DEBUG_VER
  CALL Display(myName//" timeElemLength")
#endif

  CALL GetValue(table=table, key="timeElemLength", VALUE=temprealvec, &
                origin=origin, stat=stat, isfound=isok)
  CALL AssertError1(isok, myname, "timeElemLength not found")

  CALL Reallocate(obj%timeElemLength, obj%totalTimeElements)

  abool = SIZE(temprealvec) .EQ. 1
  IF (abool) THEN
    obj%timeElemLength = temprealvec(1)
  ELSE
    isok = SIZE(temprealvec) .EQ. obj%totalTimeElements
    CALL AssertError1(isok, myname, "timeElemLength should have "// &
                      "totalTimeElements values")
    obj%timeElemLength(:) = temprealvec(1:obj%totalTimeElements)
  END IF

!INFO: alpha
#ifdef DEBUG_VER
  CALL Display(myName//" alpha")
#endif

  CALL GetValue(table=table, key="alpha", VALUE=temprealvec, &
                origin=origin, stat=stat, isfound=isok)
  CALL AssertError1(isok, myname, "alpha not found")

#ifdef DEBUG_VER
  CALL Display(obj%totalTimeElements, myname//" totalTimeElements: ")
#endif

  CALL Reallocate(obj%alpha, obj%totalTimeElements)

  abool = SIZE(temprealvec) .EQ. 1
  IF (abool) THEN
    obj%alpha(:) = temprealvec(1)
  ELSE
    isok = SIZE(temprealvec) .EQ. obj%totalTimeElements
    CALL AssertError1(isok, myname, "alpha should have "// &
                      "totalTimeElements values")
    obj%alpha(:) = temprealvec(1:obj%totalTimeElements)
  END IF

#ifdef DEBUG_VER
  CALL Display(myName//" natural Frequency")
#endif

  CALL GetValue(table=table, key="naturalFrequency", &
                VALUE=obj%naturalFrequency, default_value=1.0_DFP, &
                origin=origin, stat=stat)

#ifdef DEBUG_VER
  CALL Display(myName//" damping Ratio")
#endif

  CALL GetValue(table=table, key="dampingRatio", &
                VALUE=obj%dampingRatio, default_value=0.0_DFP, &
                origin=origin, stat=stat)

! INFO: saveData
#ifdef DEBUG_VER
  CALL Display(myName//" saveData")
#endif

  obj%saveData = .TRUE.
  CALL toml_get(table, "saveData", array, origin=origin, stat=stat)
  CALL toml_get(array, tempboolvec, origin=origin, stat=stat)

  IF (SIZE(tempboolvec) .EQ. 4) THEN
    obj%saveData = tempboolvec
  ELSE
    CALL AssertError1(.FALSE., myname, "saveData should have 4 values")
  END IF
  array => NULL()
  DEALLOCATE (tempboolvec)

! INFO: plotData
#ifdef DEBUG_VER
  CALL Display(myName//" plotData")
#endif

  obj%plotData = .TRUE.
  CALL toml_get(table, "plotData", array, origin=origin, stat=stat)
  CALL toml_get(array, tempboolvec, origin=origin, stat=stat)

  IF (SIZE(tempboolvec) .EQ. 3) THEN
    obj%plotData = tempboolvec
  ELSE
    CALL AssertError1(.FALSE., myname, "plotData should have 3 values")
  END IF
  array => NULL()

!INFO: initialVel
  astr = "initialVel"
#ifdef DEBUG_VER
  CALL Display(myName//astr%chars())
#endif
  node => NULL()
  CALL toml_get(table, astr%chars(), node, origin=origin, requested=.FALSE., &
                stat=stat)
  isok = ASSOCIATED(node)
  IF (isok) THEN
    ALLOCATE (obj%initialVel)
    CALL obj%initialVel%ImportFromToml(table=node)
  END IF
  node => NULL()

!INFO: initialDisp
  astr = "initialDisp"
#ifdef DEBUG_VER
  CALL Display(myName//astr%chars())
#endif
  node => NULL()
  CALL toml_get(table, astr%chars(), node, origin=origin, requested=.FALSE., &
                stat=stat)
  isok = ASSOCIATED(node)
  IF (isok) THEN
    ALLOCATE (obj%initialDisp)
    CALL obj%initialDisp%ImportFromToml(table=node)
  END IF
  node => NULL()

#ifdef DEBUG_VER
  CALL e%RaiseInformation(modName//'::'//myName//' - '// &
                          '[END] ')
#endif
END SUBROUTINE obj_ImportFromToml1

!----------------------------------------------------------------------------
!                                                            ImportFromToml
!----------------------------------------------------------------------------

!> author: Shion Shimizu
! date:   2024-05-11
! summary:  Initiate kernel from the toml file

SUBROUTINE obj_ImportFromToml2(obj, tomlName, afile,  &
  & filename, printToml)
  CLASS(VSTSDOF_), INTENT(INOUT) :: obj
  CHARACTER(*), INTENT(IN) :: tomlName
  TYPE(TxtFile_), OPTIONAL, INTENT(INOUT) :: afile
  CHARACTER(*), OPTIONAL, INTENT(IN) :: filename
  LOGICAL(LGT), OPTIONAL, INTENT(IN) :: printToml
  ! internal variables
  CHARACTER(*), PARAMETER :: myName = "obj_ImportFromToml2()"
  TYPE(toml_table), ALLOCATABLE :: table
  TYPE(toml_table), POINTER :: node
  INTEGER(I4B) :: origin, stat

#ifdef DEBUG_VER
  CALL e%RaiseInformation(modName//'::'//myName//' - '// &
    & '[START]')
#endif

  CALL GetValue(table=table, afile=afile, filename=filename)

  node => NULL()
  CALL toml_get(table, tomlName, node, origin=origin, requested=.FALSE.,  &
      & stat=stat)

  IF (.NOT. ASSOCIATED(node)) THEN
    CALL e%RaiseError(modName//'::'//myName//' - '// &
      & '[CONFIG ERROR] :: following error occured while reading '//  &
      & 'the toml file :: cannot find ['//tomlName//"] table in config.")
  END IF

  CALL obj%ImportFromToml(table=node)

#ifdef DEBUG_VER
  IF (PRESENT(printToml)) THEN
    CALL Display(toml_serialize(node), "toml config = "//CHAR_LF,  &
      & unitNo=stdout)
  END IF

  CALL e%RaiseInformation(modName//'::'//myName//' - '// &
    & '[END]')
#endif

END SUBROUTINE obj_ImportFromToml2

!----------------------------------------------------------------------------
!                                                                   VSTFEM
!----------------------------------------------------------------------------

SUBROUTINE obj_VSTFEM(obj, order, alpha)
  CLASS(VSTSDOF_), INTENT(INOUT) :: obj
  REAL(DFP), OPTIONAL, INTENT(IN) :: alpha
  INTEGER(I4B), OPTIONAL, INTENT(IN) :: order

  TYPE(ElemShapeData_) :: elemsd1, elemsd2
  TYPE(ReferenceLine_) :: refline
  TYPE(QuadraturePoint_) :: quad
  REAL(DFP), ALLOCATABLE :: pt(:, :), wt(:), temp(:)
  REAL(DFP), ALLOCATABLE :: c_t(:, :), m_t(:, :), inv_c_t(:, :)
  REAL(DFP), ALLOCATABLE :: tildeT_uv(:, :), tildeT_v(:, :), temp_N(:, :)
  REAL(DFP), ALLOCATABLE :: kmat(:, :), coef(:, :)
  INTEGER(I4B) :: ii, jj, kk, aint, nnt, order0
  REAL(DFP), ALLOCATABLE :: xij(:, :), int_points(:, :)
  REAL(DFP) :: alpha0, ja
  LOGICAL(LGT) :: isOK

  INTEGER(I4B), PARAMETER :: ipT = Equidistance, basisT = Monomial
  CHARACTER(*), PARAMETER :: quadT = "GaussLegendre"
  CHARACTER(*), PARAMETER :: myName = "obj_VSTFEM"

  alpha0 = Input(default=default_alpha, option=alpha)
  order0 = Input(default=default_order, option=order)

  isOK = order0 .GT. 0_I4B
  IF (.NOT. isOK) THEN
    CALL e%RaiseError(modName//'::'//myName//' - '// &
      & '[INTERNAL ERROR] :: order should be greater 0')
    RETURN
  END IF

  !----------------------------------------------------------------------------
  !                                                                uvst
  !----------------------------------------------------------------------------

  ! order + 1 should be equivalent to nnt
  nnt = order0 + 1
  CALL Reallocate(c_t, nnt, nnt)
  CALL Reallocate(inv_c_t, nnt, nnt)
  CALL reallocate(kmat, nnt, nnt)
  CALL reallocate(coef, nnt, nnt)
  CALL Reallocate(m_t, nnt, nnt)

  CALL Reallocate(obj%ct, nnt, nnt)
  CALL Reallocate(obj%mt, nnt, nnt)
  CALL Reallocate(obj%kt, nnt, nnt)

  CALL Initiate(obj=refline, nsd=1_I4B)
  CALL Initiate(obj=quad, refelem=refline, &
                order=2 * order0, QuadratureType=quadT)
  CALL GetQuadraturepoints(obj=quad, points=pt, weights=wt)
  CALL Initiate(obj=elemsd1, quad=quad, refelem=refline, &
                BaseContinuity=TypeH1, baseInterpolation=TypeLIP, &
                ipType=ipT, basisType=basisT, order=order0)

  DO ii = 1, nnt
    c_t = c_t + wt(ii) &
         & * OUTERPROD(elemsd1%N(:, ii), elemsd1%dNdXi(:, 1, ii))
    m_t = m_t + wt(ii) &
         & * OUTERPROD(elemsd1%N(:, ii), elemsd1%N(:, ii))
  END DO

  ! contribution from jump term
  c_t(1, 1) = c_t(1, 1) + 1.0_DFP ! coef for mass
  CALL Inv(invA=inv_c_t, A=c_t)
  coef = MATMUL(inv_c_t, m_t * 0.50_DFP)

  tildeT_uv = elemsd1%N
  tildeT_uv = 0.0_DFP

  ! coef for stiff
  DO ii = 1, order + 1
    tildeT_uv(:, ii) = tildeT_uv(:, ii) + MATMUL(elemsd1%N(:, ii), coef)
  END DO

  temp_N = elemsd1%N
  tildeT_v = elemsd1%N
  tildeT_v = 0.0_DFP
  CALL Reallocate(xij, 1, 2)
  aint = CEILING((REAL(order) + 1.0) / 2.0)

  IF (aint .EQ. 1) aint = 2

  CALL Reallocate(int_points, 2, aint)
  CALL Initiate(obj=quad, refelem=refline, &
                order=aint, QuadratureType=quadT)
  CALL Initiate(obj=elemsd1, quad=quad, refelem=refline, &
                BaseContinuity=TypeH1, baseInterpolation=TypeLIP, &
                ipType=ipT, basisType=basisT, order=1_I4B)

  xij(1, 1) = -1.0_DFP
  int_points(2, :) = quad%points(2, :)

  DO ii = 1, SIZE(pt, 2)
    xij(1, 2) = pt(1, ii)
    CALL set(obj=elemsd1, val=xij, N=elemsd1%N, dNdXi=elemsd1%dNdXi)
    int_points(1, :) = elemsd1%coord(1, :)
    ja = elemsd1%jacobian(1, 1, 1)
    CALL Initiate(obj=quad, points=int_points)
    CALL Initiate(obj=elemsd2, quad=quad, refelem=refline, &
                  BaseContinuity=TypeH1, baseInterpolation=TypeLIP, &
                  ipType=ipT, basisType=basisT, order=order0)

    DO jj = 1, nnt
      DO kk = 1, SIZE(int_points, 2)
        tildeT_v(jj, ii) = tildeT_v(jj, ii) + &
                         0.50_DFP * ja * int_points(2, kk) * elemsd2%N(jj, kk)
      END DO
    END DO
  END DO

  DO ii = 1, nnt
    temp = (1.0_DFP - alpha0) * tildeT_v(:, ii) + alpha0 * tildeT_uv(:, ii)
    kmat = kmat + wt(ii) * OUTERPROD(temp_N(:, ii), temp)
  END DO

  obj%ct = c_t
  obj%mt = m_t * 0.50_DFP
  obj%kt = kmat * 0.50_DFP

  CALL Reallocate(obj%rhs_u1, 3, nnt)
  CALL Reallocate(obj%rhs_v1, 3, nnt)
  CALL Reallocate(obj%disp, 2 + nnt)
  CALL Reallocate(obj%vel, 2 + nnt)

  obj%rhs_u1 = 0.0_DFP
  obj%rhs_v1 = 0.0_DFP
  DO ii = 1, nnt
    obj%rhs_u1(3, ii) = -0.5_DFP * DOT_PRODUCT(wt(:), temp_N(ii, :))
  END DO
  obj%rhs_v1(1, 1) = 1.0_DFP

  obj%disp = 0.0_DFP
  obj%vel = 0.0_DFP
  obj%disp(1) = 1.0_DFP
  obj%vel(4) = 1.0_DFP

  DO ii = 1, nnt
    obj%disp(2 + ii) = coef(2, ii)
  END DO

END SUBROUTINE obj_VSTFEM

!----------------------------------------------------------------------------
!
!----------------------------------------------------------------------------

SUBROUTINE obj_Set(obj)
  CLASS(VSTSDOF_), INTENT(inout) :: obj
#ifdef DEBUG_VER
  CHARACTER(*), PARAMETER :: myName = "obj_Set()"
#endif

#ifdef DEBUG_VER
  CALL e%RaiseInformation(modName//'::'//myName//' - '// &
                          '[START] ')
#endif

  obj%currentTimeStep = 1
  obj%currentTime = obj%timeDomain(1)

! CALL SetTotalDOFTime(obj=obj)

  ALLOCATE (obj%uvar(obj%totalTimeElements))
  ALLOCATE (obj%vvar(obj%totalTimeElements))
  ALLOCATE (obj%avar(obj%totalTimeElements))

  obj%uvar = 0.0_DFP
  obj%vvar = 0.0_DFP
  obj%avar = 0.0_DFP

#ifdef DEBUG_VER
  CALL e%RaiseInformation(modName//'::'//myName//' - '// &
                          '[END] ')
#endif

END SUBROUTINE obj_Set
!----------------------------------------------------------------------------
!
!----------------------------------------------------------------------------

SUBROUTINE obj_getTanmat(obj, timeElemNum)
  CLASS(VSTSDOF_), INTENT(INOUT) :: obj
  INTEGER(I4B), INTENT(IN) :: timeElemNum
  INTEGER(I4B) :: nnt, order
  REAL(DFP) :: alpha, scale, dt, dts, dt_by_2, &
               dts_by_2

  order = obj%timeOrder(obj%currentTimeStep)
  nnt = order + 1_I4B
  alpha = obj%alpha(obj%currentTimeStep)
  dt = obj%timeElemLength(timeElemNum)
  dts = dt * dt

  CALL Reallocate(obj%tanmat, nnt, nnt)
  CALL Reallocate(obj%rhs, nnt)
  CALL Reallocate(obj%sol, nnt)

  CALL obj%VSTFEM(order, alpha)

  obj%tanmat = obj%ct

  scale = dt * 2.0_DFP * obj%naturalFrequency * obj%dampingRatio
  obj%tanmat = obj%tanmat + scale * obj%mt

  scale = dts * obj%naturalFrequency**2
  obj%tanmat = obj%tanmat + scale * obj%kt

END SUBROUTINE obj_getTanmat

!----------------------------------------------------------------------------
!
!----------------------------------------------------------------------------

SUBROUTINE obj_getRHS(obj, timeElemNum)
  CLASS(VSTSDOF_), INTENT(INOUT) :: obj
  INTEGER(I4B), INTENT(IN) :: timeElemNum
  INTEGER(I4B) :: nnt
  REAL(DFP) :: dt, scale

  nnt = obj%timeOrder(obj%currentTimeStep) + 1_I4B
  dt = obj%timeElemLength(timeElemNum)

  obj%rhs = 0.0_DFP
  obj%rhs(1) = obj%rhs(1) + obj%rhs_v1(1, 1) * obj%v0

  scale = dt * obj%naturalFrequency**2
  obj%rhs = obj%rhs + scale * obj%rhs_u1(3, :) * obj%u0

END SUBROUTINE obj_getRHS

!----------------------------------------------------------------------------
!
!----------------------------------------------------------------------------

SUBROUTINE obj_solve(obj)
  CLASS(VSTSDOF_), INTENT(INOUT) :: obj

  CALL Solve(X=obj%sol, A=obj%tanmat, B=obj%rhs)

END SUBROUTINE obj_solve

!----------------------------------------------------------------------------
!
!----------------------------------------------------------------------------

SUBROUTINE obj_update(obj)
  CLASS(VSTSDOF_), INTENT(inout) :: obj
  INTEGER(I4B) :: ii, nnt
  REAL(DFP) :: dt

  nnt = obj%timeOrder(obj%currentTimeStep) + 1_I4B

  dt = obj%timeElemLength(obj%currentTimeStep)
  obj%currentTime = obj%currentTime + dt
  obj%currentTimeStep = obj%currentTimeStep + 1_I4B

  obj%v0 = 0.0_DFP
  obj%u0 = obj%u0 * obj%disp(1)

  DO ii = 1, nnt
    obj%v0 = obj%v0 + obj%vel(2_I4B + ii) * obj%sol(ii)
    obj%u0 = obj%u0 + dt * obj%disp(2_I4B + ii) * obj%sol(ii)
  END DO

  obj%uvar(obj%currentTimeStep - 1) = obj%u0
  obj%vvar(obj%currentTimeStep - 1) = obj%v0

END SUBROUTINE obj_update

!----------------------------------------------------------------------------
!                                                         SetInitialVelocity
!----------------------------------------------------------------------------

SUBROUTINE obj_SetInitialVelocity(obj)
  CLASS(VSTSDOF_), INTENT(inout) :: obj
#ifdef DEBUG_VER
  CHARACTER(*), PARAMETER :: myName = "obj_SetInitialVelocity()"
#endif

  LOGICAL(LGT) :: isok

#ifdef DEBUG_VER
  CALL e%RaiseInformation(modName//'::'//myName//' - '// &
                          '[START] ')
#endif

  isok = ASSOCIATED(obj%initialVel)
  IF (.NOT. isok) RETURN

  CALL obj%initialVel%Get(val=obj%v0)

#ifdef DEBUG_VER
  CALL e%RaiseInformation(modName//'::'//myName//' - '// &
                          '[END] ')
#endif

END SUBROUTINE obj_SetInitialVelocity

!----------------------------------------------------------------------------
!                                                      SetInitialDisplacement
!----------------------------------------------------------------------------

SUBROUTINE obj_SetInitialDisplacement(obj)
  CLASS(VSTSDOF_), INTENT(inout) :: obj
#ifdef DEBUG_VER
  CHARACTER(*), PARAMETER :: myName = "obj_SetInitialDisplacement()"
#endif

  LOGICAL(LGT) :: isok

#ifdef DEBUG_VER
  CALL e%RaiseInformation(modName//'::'//myName//' - '// &
                          '[START] ')
#endif

  isok = ASSOCIATED(obj%initialDisp)
  IF (.NOT. isok) RETURN

  CALL obj%initialDisp%Get(val=obj%u0)

#ifdef DEBUG_VER
  CALL e%RaiseInformation(modName//'::'//myName//' - '// &
                          '[END] ')
#endif

END SUBROUTINE obj_SetInitialDisplacement

!----------------------------------------------------------------------------
!
!----------------------------------------------------------------------------

SUBROUTINE obj_Run(obj)
  CLASS(VSTSDOF_), INTENT(inout) :: obj
#ifdef DEBUG_VER
  CHARACTER(*), PARAMETER :: myName = "obj_Run()"
#endif

  INTEGER(I4B) :: ielTime

#ifdef DEBUG_VER
  CALL e%RaiseInformation(modName//'::'//myName//' - '// &
                          '[START] ')
#endif

  CALL obj%SetInitialVelocity()
  CALL obj%SetInitialDisplacement()

  DO ielTime = 1, obj%totalTimeElements
    CALL Display(obj%currentTime, myname//" t1: ")
    CALL obj%GetTanmat(timeElemNum=ielTime)
    CALL obj%GetRHS(timeElemNum=ielTime)
    CALL obj%Solve()
    CALL obj%Update()
  END DO

  CALL obj%WriteData()

#ifdef DEBUG_VER
  CALL e%RaiseInformation(modName//'::'//myName//' - '// &
                          '[END] ')
#endif

END SUBROUTINE obj_Run

!----------------------------------------------------------------------------
!
!----------------------------------------------------------------------------

SUBROUTINE obj_WriteData(obj)
  CLASS(VSTSDOF_), INTENT(inout) :: obj
#ifdef DEBUG_VER
  CHARACTER(*), PARAMETER :: myName = "obj_WriteData()"
#endif

  REAL(DFP) :: xlim(2), ylim(2)

  REAL(DFP), ALLOCATABLE :: timeData(:)

  INTEGER(I4B) :: ii

  CHARACTER(:), ALLOCATABLE :: filename_disp, filename_vel, filename_acc, &
                               aline, filename_data
  LOGICAL(LGT) :: abool1, abool2, abool3

  abool1 = ANY(obj%saveData)
  abool2 = ANY(obj%plotData)

  IF (.NOT. abool1 .AND. .NOT. abool2) RETURN

#ifdef DEBUG_VER
  CALL e%RaiseInformation(modName//'::'//myName//' - '// &
                          '[START] ')
#endif

  ALLOCATE (timeData(obj%totalTimeElements))
  timeData(1) = obj%timeElemLength(1)
  DO ii = 2, obj%totalTimeElements
    timeData(ii) = timedata(ii - 1) + obj%timeElemLength(ii)
  END DO

  filename_disp = obj%result_dir//CHAR_SLASH//obj%filename//'_disp'

  filename_vel = obj%result_dir//CHAR_SLASH//obj%filename//'_vel'

  filename_acc = obj%result_dir//CHAR_SLASH//obj%filename//'_acc'

  filename_data = obj%result_dir//CHAR_SLASH//obj%filename//'_data'

  abool1 = obj%saveData(1) .OR. obj%plotData(1) .OR. obj%saveData(4)
  abool2 = obj%saveData(2) .OR. obj%plotData(2) .OR. obj%saveData(4)
  abool3 = obj%saveData(3) .OR. obj%plotData(3) .OR. obj%saveData(4)

! disp
  IF (obj%saveData(1)) THEN
#ifdef DEBUG_VER
    CALL Display("Writing data to file: "//filename_disp//".csv")
#endif
    CALL obj%dispfile%Initiate(filename=filename_disp//".csv", unit=100, &
                               status="REPLACE", action="WRITE", &
                               comment="#", separator=",")
    CALL obj%dispfile%OPEN()

    aline = "time, disp"
    CALL obj%dispfile%WRITE(aline)

    DO ii = 1, obj%totalTimeElements
      aline = tostring(timedata(ii))//", "//tostring(obj%uvar(ii))
      CALL obj%dispfile%WRITE(aline)
    END DO
    CALL obj%dispfile%DEALLOCATE()
  END IF

! vel
  IF (obj%saveData(2)) THEN
#ifdef DEBUG_VER
    CALL Display("Writing data to file: "//filename_vel//".csv")
#endif
    CALL obj%velfile%Initiate(filename=filename_vel//".csv", unit=100, &
                              status="REPLACE", action="WRITE", &
                              comment="#", separator=",")
    CALL obj%velfile%OPEN()

    aline = "time, vel"
    CALL obj%velfile%WRITE(aline)

    DO ii = 1, obj%totalTimeElements
      aline = tostring(timedata(ii))//", "//tostring(obj%vvar(ii))
      CALL obj%velfile%WRITE(aline)
    END DO
    CALL obj%velfile%DEALLOCATE()
  END IF

! acc
  IF (obj%saveData(3)) THEN
#ifdef DEBUG_VER
    CALL Display("Writing data to file: "//filename_acc//".csv")
#endif
    CALL obj%accfile%Initiate(filename=filename_acc//".csv", unit=100, &
                              status="REPLACE", action="WRITE", &
                              comment="#", separator=",")
    CALL obj%accfile%OPEN()

    aline = "time, acc"
    CALL obj%accfile%WRITE(aline)

    DO ii = 1, obj%totalTimeElements
      aline = tostring(timedata(ii))//", "//tostring(obj%avar(ii))
      CALL obj%accfile%WRITE(aline)
    END DO
    CALL obj%accfile%DEALLOCATE()
  END IF

! write all data
  IF (obj%saveData(4)) THEN
#ifdef DEBUG_VER
    CALL Display("Writing data to file: "//filename_data//".csv")
#endif
    CALL obj%datafile%Initiate(filename=filename_data//".csv", &
                 status="REPLACE", action="WRITE", comment="#", separator=",")
    CALL obj%datafile%OPEN()

    aline = "time, disp, vel, acc"
    CALL obj%datafile%WRITE(aline)
    DO ii = 1, obj%totalTimeElements
      aline = tostring(timedata(ii))//", "//tostring(obj%uvar(ii))// &
              ", "//tostring(obj%vvar(ii))//", "//tostring(obj%avar(ii))
      CALL obj%accfile%WRITE(aline)
    END DO
    CALL obj%datafile%DEALLOCATE()
  END IF

#ifdef DEBUG_VER
  CALL Display("Done writing files csvfiles")
#endif

! plotting
  IF (obj%plotData(1)) THEN
    CALL obj%plot%filename(filename_disp//'.plt')
    CALL obj%plot%options('set terminal pngcairo; set output "' &
                          //filename_disp//'.png"')
    xlim = obj%timeDomain
    ylim = [MINVAL(obj%uvar), MAXVAL(obj%uvar)]
    xlim(1) = xlim(1) - 0.1 * (xlim(2) - xlim(1))
    xlim(2) = xlim(2) + 0.1 * (xlim(2) - xlim(1))
    ylim(1) = ylim(1) - 0.1 * (ylim(2) - ylim(1))
    ylim(2) = ylim(2) + 0.1 * (ylim(2) - ylim(1))

    CALL obj%plot%xlim(xlim)
    CALL obj%plot%ylim(ylim)
    CALL obj%plot%xlabel('t')
    CALL obj%plot%ylabel('u')
    CALL obj%plot%plot(x1=timedata, y1=obj%uvar)
    CALL obj%plot%reset()

  END IF

  IF (obj%plotData(2)) THEN
    CALL obj%plot%filename(filename_vel//'.plt')
    CALL obj%plot%options('set terminal pngcairo; set output "'// &
                          filename_vel//'.png"')
    xlim = obj%timeDomain
    ylim(1) = MINVAL(obj%vvar)
    ylim(2) = MAXVAL(obj%vvar)
    xlim(1) = xlim(1) - 0.1 * (xlim(2) - xlim(1))
    xlim(2) = xlim(2) + 0.1 * (xlim(2) - xlim(1))
    ylim(1) = ylim(1) - 0.1 * (ylim(2) - ylim(1))
    ylim(2) = ylim(2) + 0.1 * (ylim(2) - ylim(1))

    CALL obj%plot%xlim(xlim)
    CALL obj%plot%ylim(ylim)
    CALL obj%plot%xlabel('t')
    CALL obj%plot%ylabel('v')
    CALL obj%plot%plot(x1=timedata, y1=obj%vvar)
    CALL obj%plot%reset()
  END IF

! IF (obj%plotData(3)) THEN
!   CALL obj%plot%filename(filename_acc//'.plt')
! CALL obj%plot%options('set terminal pngcairo; set output "'//filename_acc//'.png"')
!   xlim = obj%timeDomain
!   ylim(1) = MINVAL(obj%avar)
!   ylim(2) = MAXVAL(obj%avar)
!   xlim(1) = xlim(1) - 0.1 * (xlim(2) - xlim(1))
!   xlim(2) = xlim(2) + 0.1 * (xlim(2) - xlim(1))
!   ylim(1) = ylim(1) - 0.1 * (ylim(2) - ylim(1))
!   ylim(2) = ylim(2) + 0.1 * (ylim(2) - ylim(1))
!
!   CALL obj%plot%xlim(xlim)
!   CALL obj%plot%ylim(ylim)
!   CALL obj%plot%xlabel('t')
!   CALL obj%plot%ylabel('a')
!   CALL obj%plot%plot(x1=timeData, y1=obj%avar)
!   CALL obj%plot%reset()
! END IF

  DEALLOCATE (timedata)

#ifdef DEBUG_VER
  CALL e%RaiseInformation(modName//'::'//myName//' - '// &
                          '[END] ')
#endif

END SUBROUTINE obj_WriteData

!----------------------------------------------------------------------------
!
!----------------------------------------------------------------------------

#include "./include/errors.F90"

END MODULE VSTSDOF
