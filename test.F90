PROGRAM main
USE GlobalData
USE Helmholtz1DFEM_Class
USE FPL
USE ExceptionHandler_Class, ONLY: e, EXCEPTION_INFORMATION
USE Display_Method
USE BaseInterpolation_Method, ONLY: BaseInterpolation_ToString
USE String_Class
USE LagrangePolynomialUtility, ONLY: InterpolationPoint, &
                                     InterpolationPoint_
USE BaseType, ONLY: elem => TypeElemNameOpt

IMPLICIT NONE

TYPE(Helmholtz1DFEM_) :: obj
CHARACTER(*), PARAMETER :: testname = "test", &
                           tomlName = "kernel", &
                           filename = "helmholtz_config.toml"
INTEGER(I4B) :: aint, elemNum, inds(2), ii
INTEGER(I4B), ALLOCATABLE :: con(:)
REAL(DFP) :: xij(1, 2)
REAL(DFP), ALLOCATABLE :: ans(:, :)
TYPE(String) :: astr

CALL e%setQuietMode(.TRUE.)

CALL FPL_Init

CALL obj%ImportFromToml(tomlName=tomlName, filename=filename)

CALL obj%Set()

elemNum = 2
ALLOCATE (con(obj%spaceOrder(elemNum) + 1))
CALL obj%InitiateConnectivity()
CALL obj%InitiateFields()
CALL obj%GetConnectivity(spaceElemNum=elemNum, &
                         ans=con, tsize=aint)

CALL Display(elemNum, " element number ::")
CALL Display(obj%spaceOrder(1), " space order ::")
CALL Display(con, "connectivity ::")

xij(1, 1) = 0.0
xij(1, 2) = 1.0
CALL obj%SetQuadForSpace(elemNum)
CALL obj%SetElemsdForSpace(elemNum, xij)

astr = BaseInterpolation_ToString(obj%ipTypeForSpace)
CALL Display(astr, " ip type ::")
ALLOCATE (ans(1, MAXVAL(obj%spaceOrder) + 1))
DO ii = 1, 5
  xij(1, 2) = xij(1, 1) + 1.0
  CALL InterpolationPoint_(order=obj%spaceOrder(ii), elemType=elem%line, &
                          ipType=obj%ipTypeForSpace, xij=xij, layout="VEFC", &
                           ans=ans, nrow=inds(1), ncol=inds(2))
  CALL Display(ans, " interpolation points ::")
  xij(1, 1) = xij(1, 2)
END DO

DEALLOCATE (con)
CALL FPL_Finalize

END PROGRAM main
