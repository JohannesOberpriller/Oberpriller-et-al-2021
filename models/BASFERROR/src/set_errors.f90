Subroutine set_process_errors(err)

use errors
implicit none


real :: err(5)
RUNOFFERR = err(1)
RESPERR = err(2)
RLEACHERR = err(3)
fLUETERR = err(4)
FRERR = err(5)

end Subroutine set_process_errors

Subroutine set_state_errors(err)

use errors
implicit none

real :: err(3)
WAERR = err(1)
NPPERR = err(2)
NSOMSERR = err(3)

end Subroutine set_state_errors
