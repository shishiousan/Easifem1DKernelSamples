PROGRAM main
USE GlobalData
USE Helmholtz1DFEM_Class
USE FPL
USE ExceptionHandler_Class, ONLY: e, EXCEPTION_INFORMATION
USE Display_Method

IMPLICIT NONE

TYPE(Helmholtz1DFEM_) :: obj
CHARACTER(*), PARAMETER :: testname = "test", &
                           tomlName = "kernel", &
                           filename = "helmholtz_config.toml"

CALL e%setQuietMode(.TRUE.)

CALL FPL_Init

CALL obj%ImportFromToml(tomlName=tomlName, filename=filename)

CALL obj%Set()

CALL obj%Run()

CALL FPL_Finalize

END PROGRAM main
