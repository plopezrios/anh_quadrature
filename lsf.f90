MODULE lsf
  IMPLICIT NONE
  PRIVATE
  PUBLIC lsf_res
  PUBLIC lsf_xpower_expval, lsf_xscale, lsf_norder, lsf_orbcoeff
  INTEGER :: lsf_norder = 0
  DOUBLE PRECISION :: lsf_xscale = 1.d0
  DOUBLE PRECISION, ALLOCATABLE :: lsf_xpower_expval(:), lsf_orbcoeff(:)

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
    DOUBLE PRECISION grid_x(nparam/2), grid_P(nparam/2), &
       &hbasis(0:lsf_norder), x
    INTEGER i, ngrid

    ! Silence compiler warnings.
    if(present(uiparm).or.present(urparm).or.present(ufparm))nf=nf

    ! Evaluate residual (fit function minus y from data) for each data point.
    ngrid = nparam/2
    grid_x = param(1:ngrid)
    grid_P = param(ngrid+1:2*ngrid)
    if (any(grid_P<=0.d0) .or. &
       &any(grid_x(1:ngrid-1)>=grid_x(2:ngrid))) then
      nf = -1
      return
    endif
    do i = 1, ngrid
      x = grid_x(i)*lsf_xscale
      call eval_hermite_poly_norm (lsf_norder, x, hbasis)
      grid_P(i) = grid_P(i)*exp(-x*x) * &
         &sum(lsf_orbcoeff(0:lsf_norder)*hbasis(0:lsf_norder))**2
    enddo ! i
    grid_P = grid_P/sum(grid_P)
    res(1:nxy) = (/ ( sum( grid_P*grid_x**dble(i-1) ) - &
       &              lsf_xpower_expval(i-1)/lsf_xscale**dble(i-1), &
       &              i=1,nxy ) /)

  END SUBROUTINE lsf_res


  SUBROUTINE eval_hermite_poly_norm (n, x, h)
    !---------------------------------------------------!
    ! Evaluate N_n(x) for n = 0:N at x=X, returning the !
    ! values in H(0:N).                                 !
    !---------------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: n
    DOUBLE PRECISION, INTENT(in) :: x
    DOUBLE PRECISION, INTENT(inout) :: h(0:n)
    ! Numerical constants.
    DOUBLE PRECISION, PARAMETER :: pi = 4.d0*atan(1.d0)
    DOUBLE PRECISION, PARAMETER :: inv_pi_one_fourth = pi**(-0.25d0)
    DOUBLE PRECISION, PARAMETER :: sqrt2 = sqrt(2.d0)
    ! Local variables.
    DOUBLE PRECISION t0, t1, t2, sqrti, sqrti_1, sqrt2_x
    INTEGER i
    if (n<0) return
    t1 = inv_pi_one_fourth
    h(0) = t1
    if (n<1) return
    sqrt2_x = sqrt2*x
    t0 = inv_pi_one_fourth*sqrt2_x
    sqrti = 1.d0
    h(1) = t0
    do i = 2, n
      t2 = t1
      t1 = t0
      sqrti_1 = sqrti
      sqrti = sqrt(dble(i))
      t0 = (sqrt2_x*t1 - sqrti_1*t2)/sqrti
      h(i) = t0
    enddo ! i
  END SUBROUTINE eval_hermite_poly_norm


END MODULE lsf
