PROGRAM main
USE easifemBase
USE tomlf, ONLY: toml_table, &
                 toml_get => get_value, &
                 toml_stat, &
                 toml_array, &
                 toml_len => len
USE TomlUtility, ONLY: GetValue
USE TxtFile_Class
USE CSVFile_Class
IMPLICIT NONE
TYPE(toml_table), ALLOCATABLE :: table
TYPE(toml_table), POINTER :: node
TYPE(toml_array), POINTER :: array
REAL(DFP), ALLOCATABLE :: vec(:)
INTEGER(I4B), ALLOCATABLE :: intvec(:)
INTEGER(I4B) :: stat, origin, iostat, ii, tsize
LOGICAL(LGT) :: isok, abool
TYPE(String) :: filename, ext
REAL(DFP) :: temp
TYPE(TxtFile_) :: atxtfile
TYPE(CSVFile_) :: acsvfile
CHARACTER(512) :: iomsg
CHARACTER(*), PARAMETER :: Key = "elasticModulus"

CALL GetValue(table=table, filename="config.toml")

node => NULL()
CALL toml_get(table, "kernel", node)

CALL toml_get(node, key, array, origin=origin, stat=stat, &
              requested=.FALSE.)

IF (stat .EQ. toml_stat%success) THEN
  tsize = toml_len(array)
  CALL Reallocate(vec, tsize)
  DO ii = 1, tsize
    CALL toml_get(array, ii, vec(ii))
  END DO
  NULLIFY (array)
  ! return
ELSE
  CALL toml_get(node, key, temp, origin=origin, stat=stat)

  IF (stat .EQ. toml_stat%success) THEN
    CALL Reallocate(vec, 1)
    vec(1) = temp
    ! return
  END IF

END IF

IF (stat .NE. toml_stat%success) THEN

  CALL toml_get(node, key, filename%raw, origin=origin, stat=stat)

  ext = filename%extension()
  SELECT CASE (ext%chars())
  CASE (".txt")
    IF (stat .EQ. toml_stat%success) THEN
      CALL atxtfile%Initiate(filename=filename%Chars(), &
                             action="READ", status="OLD")
      CALL atxtfile%OPEN()
      CALL atxtfile%READ(val=vec, iostat=iostat, iomsg=iomsg)

      abool = iostat .NE. 0 .AND. (.NOT. atxtfile%isEOF())
      IF (abool) THEN
        STOP "error occur "
        filename = ""
        RETURN
      END IF

      CALL atxtfile%DEALLOCATE()
      filename = ""
    END IF
  CASE (".csv")
    IF (stat .EQ. toml_stat%success) THEN
      CALL acsvfile%Initiate(filename=filename%Chars(), &
                             action="READ", status="OLD", &
                             delimiter=",")
      CALL acsvfile%OPEN()
      CALL acsvfile%SetHeaderIndx(1)
      CALL acsvfile%READ()
      CALL Display(acsvfile%header(1), "header(1) ::")
      CALL Display(acsvfile%header(2), "header(2) ::")
      CALL Display(acsvfile%header(3), "header(3) ::")
      CALL Display(acsvfile%Getncols(), "number of columns")
      CALL acsvfile%get(1, val=intvec)
      CALL acsvfile%get(3, val=vec)
      CALL Display(intvec, "start::")
      CALL Display(vec, "values::")
      CALL acsvfile%DEALLOCATE()
      filename = ""
    END IF
  END SELECT

END IF

IF (ALLOCATED(vec)) THEN
  CALL Display(vec, " key ::")
  CALL Display(stat, " stat ::")
  DEALLOCATE (vec)
END IF

END PROGRAM main

