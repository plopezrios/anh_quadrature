MODULE lsf
  IMPLICIT NONE
  PRIVATE
  PUBLIC lsf_res
  PUBLIC lsf_xpower_expval
  DOUBLE PRECISION, ALLOCATABLE :: lsf_xpower_expval(:)

CONTAINS


  SUBROUTINE lsf_res (nxy, nparam, param, nf, res, uiparm, urparm, ufparm)
    !--------------------------------!
    ! Evaluate residuals for NL2SOL. !
    !--------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nxy,nparam
    DOUBLE PRECISION,INTENT(in) :: param(:)
    INTEGER,INTENT(inout) :: nf
    DOUBLE PRECISION,INTENT(out) :: res(:)
    INTEGER,INTENT(inout),OPTIONAL :: uiparm(:)
    DOUBLE PRECISION,INTENT(inout),OPTIONAL :: urparm(:),ufparm
    DOUBLE PRECISION grid_x(nparam/2), grid_P(nparam/2)
    INTEGER i

    ! Silence compiler warnings.
    if(present(uiparm).or.present(urparm).or.present(ufparm))nf=nf

    ! Evaluate residual (fit function minus y from data) for each data point.
    grid_x = param(1:nparam/2)
    grid_P = param(nparam/2+1:nparam)
    if (any(grid_P<=0.d0)) then
      nf = -1
      return
    endif
    res(1:nxy) = (/ ( sum( grid_P*grid_x**dble(i-1) ) - &
       &              lsf_xpower_expval(i-1), i=1,nxy ) /)

  END SUBROUTINE lsf_res


END MODULE lsf
