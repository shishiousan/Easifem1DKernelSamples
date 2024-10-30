PROGRAM main
USE easifemBase
IMPLICIT NONE
INTEGER(I4B) :: aint
LOGICAL(LGT) :: abools(4)
CALL execute_command_line(command="gnuplot ./results/helm-1_disp.plt 2>/dev/null", &
                          exitstat=aint)

IF (aint .NE. 0) THEN
  CALL Display("gnuplot error")
END IF

abools = [.TRUE., .TRUE., .FALSE., .FALSE.]

aint = COUNT(abools)

CALL Display(aint, " the number of true in abools ::")

END PROGRAM main

