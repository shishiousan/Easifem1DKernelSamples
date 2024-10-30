PROGRAM main
USE GlobalData
USE ElastoDynamics1DSDFEM_Class
USE FPL
USE ExceptionHandler_Class, ONLY: e, EXCEPTION_INFORMATION
USE Display_Method

IMPLICIT NONE

TYPE(ElastoDynamics1DSDFEM_) :: obj
CHARACTER(*), PARAMETER :: testname = "test", &
                           tomlName = "kernel", &
                           filename = "config_ricker.toml"

CALL FPL_Init

CALL obj%ImportFromToml(tomlName=tomlName, filename=filename)

IF (obj%verbosity .EQ. 0) THEN
  CALL e%setQuietMode(EXCEPTION_INFORMATION, .TRUE.)
END IF

CALL obj%Set()
CALL obj%Display(testname)

CALL obj%Run()

CALL FPL_Finalize

END PROGRAM main
