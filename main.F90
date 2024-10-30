PROGRAM main
USE GlobalData
USE ElastoDynamics1DSTFEM_Class
USE FPL
USE ExceptionHandler_Class, ONLY: e, EXCEPTION_INFORMATION
USE Display_Method
USE ElemshapeData_Method

IMPLICIT NONE

TYPE(ElastoDynamics1DSTFEM_) :: obj
INTEGER(I4B) :: ii
CHARACTER(*), PARAMETER :: testname = "test", &
                           tomlName = "kernel", &
                           filename = "config.toml"

CALL FPL_Init

CALL obj%ImportFromToml(tomlName=tomlName, filename=filename)

IF (obj%verbosity .EQ. 0) THEN
  CALL e%setQuietMode(EXCEPTION_INFORMATION, .TRUE.)
END IF

CALL obj%Set()
! CALL obj%display("st kernel")

! DO ii = 1, obj%totalSpaceElements
!   CALL Display(obj%elasticModulus(ii), "E of "//tostring(ii)//"::", &
!                advance="NO")
!   CALL Display(obj%density(ii), "density of "//tostring(ii)//"::", &
!                advance="NO")
!   CALL Display(obj%spaceElemLength(ii), "length of "//tostring(ii)//"::", &
!                advance="NO")
!   CALL Display(obj%spaceOrder(ii), "order of "//tostring(ii)//"::")
! END DO
! CALL Display(obj%rayleighBeta, " rayleighBeta ::")
! CALL Display(obj%rayleighAlpha, " rayleighAlpha ::")
CALL obj%Run()

CALL FPL_Finalize

END PROGRAM main
