PROGRAM PLASMA 
USE NUMERICS
USE UTILS
USE GODUNOV
USE SEMILAG

IMPLICIT NONE

!CALL TEST_INIT(0.5_RP,1000,1000)
CALL INITILISATION(25.5_rp,750,750)
!CALL GS
CALL SL_SOLVER


END PROGRAM PLASMA