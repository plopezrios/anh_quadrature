PROGRAM anh_quadrature
  !-----------------------------------------------------------------!
  ! ANH_QUADRATURE                                                  !
  ! ==============                                                  !
  ! Given a one-dimensional potential V(x), obtain the ground state !
  ! wave function as a linear combination of eigenfunctions of the  !
  ! harmonic oscillator and evaluate the optimal quadrature grid    !
  ! for integration of expectations values.                         !
  !                                                                 !
  ! PLR 08.2017                                                     !
  !-----------------------------------------------------------------!
  IMPLICIT NONE

  call main()


CONTAINS


  SUBROUTINE main()
    !--------------!
    ! Main driver. !
    !--------------!
    IMPLICIT NONE
    ! System definition.
    INTEGER norder_v
    DOUBLE PRECISION, ALLOCATABLE :: vcoeff(:)
    DOUBLE PRECISION omega, inv_sqrt_omega
    ! Variables defining ground state of Hamiltonian.
    INTEGER norder
    DOUBLE PRECISION e0
    DOUBLE PRECISION, ALLOCATABLE :: orbcoeff(:)
    ! Maximum expansion order.
    INTEGER, PARAMETER :: MAX_NORDER = 200
    ! The wave function converges at a given expansion order if its
    ! variational energy is within relative REL_TOL or absolute tolerance
    ! ABS_TOL of all of the variational energies obtained with the previous
    ! NITER_CONVERGE expansion orders.
    INTEGER, PARAMETER :: NITER_CONVERGE = 6
    DOUBLE PRECISION, PARAMETER :: REL_TOL = sqrt(epsilon(1.d0))
    DOUBLE PRECISION, PARAMETER :: ABS_TOL = 1.d-6
    ! Variables for plotting.
    DOUBLE PRECISION t1, orbnorm
    DOUBLE PRECISION, ALLOCATABLE :: vcoeff_unnorm(:), hbasis(:)
    DOUBLE PRECISION u, x, vu, psiu
    ! Numerical constants.
    DOUBLE PRECISION, PARAMETER :: log2 = log(2.d0)
    DOUBLE PRECISION, PARAMETER :: pi = 4.d0*atan(1.d0)
    ! Mist local variables.
    CHARACTER(2048) line
    INTEGER i, ierr, ivec
    DOUBLE PRECISION e0vec(0:NITER_CONVERGE-1)

    ! FIXME - ask for V in natural polynomial form, find optimal OMEGA (by some
    ! criterion) and re-represent in Hermite polynomials internally.

    ! Get scaling factor.
    write(6,*) 'Enter omega (a.u.):'
    read(5,*,iostat=ierr) omega
    if (ierr/=0) call quit()
    if (omega<=0.d0) call quit('Omega must be positive.')
    inv_sqrt_omega = 1.d0/sqrt(omega)

    ! Get V coefficients and norder.
    write(6,*) 'Enter coefficients of expansion of V(x/sqrt(omega)) in &
       &Hermite polynomials,'
    write(6,*)'V(x/sqrt(omega)) = sum_k c_k H_k(x)    (one line):'
    read(5,'(a)') line
    norder_v = -1
    do
      read(line,*,iostat=ierr) (t1, i=0,norder_v+1)
      if(ierr/=0)exit
      norder_v = norder_v+1
    enddo
    if (norder_v<0) call quit('Could not parse coefficients.')
    allocate(vcoeff(0:norder_v))
    read(line,*) vcoeff(0:norder_v)

    ! Get good approximation to ground state wave function.
    ! FIXME/wtf - all LAPACK routines hang at norder=202 regardless of vcoeff
    e0vec(0:NITER_CONVERGE-1) = (/ (dble(i)*1.d7, i=1,NITER_CONVERGE) /)
    ivec = 0
    do norder = max(norder_v,2), MAX_NORDER
      if (allocated(orbcoeff)) deallocate (orbcoeff)
      allocate (orbcoeff(0:norder))
      orbcoeff = 0.d0
      call get_ground_state (norder, norder_v, vcoeff, omega, e0, orbcoeff)
      if (all(abs(e0vec-e0)<REL_TOL*abs(e0))) exit
      if (all(abs(e0vec-e0)<ABS_TOL)) exit
      ivec = modulo(ivec-1,NITER_CONVERGE)
      e0vec(ivec) = e0
    enddo ! norder
    if (norder>MAX_NORDER) call quit('Failed to converge energy to target &
       &accuracy.')

    ! Report ground-state energy.
    write(6,*) 'Ground state energy converges at expansion order '//&
       &trim(i2s(norder))//':'
    write(6,*) '  E0 = ', e0

    ! Plot V(x) and Psi(x).
    allocate (hbasis(0:norder), vcoeff_unnorm(0:norder_v))
    orbnorm = (omega/pi)**0.25d0
    vcoeff_unnorm(0:norder_v) = &
       &(/ ( vcoeff(i) * exp( 0.5d0 * (eval_log_fact(i)+dble(i)*log2) ), &
       &     i=0,norder_v ) /)
    do i = -100, 100
      x = dble(i)/dble(100) * 5.d0
      u = x*inv_sqrt_omega
      call eval_hermite_poly_norm (norder, x, hbasis)
      vu = sum(vcoeff_unnorm(0:norder_v)*hbasis(0:norder_v))
      t1 = orbnorm*exp(-0.5d0*x*x)
      psiu = t1*sum(orbcoeff(0:norder)*hbasis(0:norder))
      write(11,*) u, vu, psiu
    enddo ! i
    deallocate (hbasis, vcoeff_unnorm)

    ! Solve for the quadrature grid.
    ! FIXME - write

  END SUBROUTINE main


  SUBROUTINE get_ground_state (norder, norder_v, vcoeff, omega, e0, orbcoeff)
    !---------------------------------------------------------!
    ! Given a one-dimensional (anharmonic) potential,         !
    !                                                         !
    !   V(x) = Sum_i VCOEFF(i)*H_i(sqrt(OMEGA)*x) ,           !
    !                                                         !
    ! where H_i is the i-th Hermite polynomial, construct the !
    ! matrix elements of the Hamiltonian                      !
    !                                                         !
    !   H(x) = -1/2 d/dx^2 + V(x) ,                           !
    !                                                         !
    ! in the basis of the eigenfunctions of the harmonic      !
    ! oscillator,                                             !
    !                                                         !
    !   phi_i(x) = (2^i i!)^-1/2 (OMEGA/pi)^1/4 *             !
    !              exp( -OMEGA * x^2 / 2) *                   !
    !              H_i( sqrt(OMEGA)*x ) ,                     !
    !                                                         !
    ! and solve for the coefficients of the normalized trial  !
    ! wave function,                                          !
    !                                                         !
    !   Psi(x) = Sum_i ORBCOEFF(i)*phi_i(x) ,                 !
    !                                                         !
    ! and for the associated variational energy E0.           !
    !---------------------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: norder, norder_v
    DOUBLE PRECISION, INTENT(in) :: vcoeff(0:norder_v), omega
    DOUBLE PRECISION, INTENT(inout) :: e0, orbcoeff(0:norder)
    ! Eigenproblem arrays.
    DOUBLE PRECISION alpha(0:norder), hmatrix(0:norder,0:norder)
    ! Buffer for logarithms of factorials.
    DOUBLE PRECISION log_fact(0:norder)
    ! LAPACK work arrays.
    DOUBLE PRECISION, ALLOCATABLE :: lapack_work(:)
    INTEGER, ALLOCATABLE :: lapack_iwork(:)
    INTEGER lapack_lwork, lapack_liwork
    ! Misc local variables.
    DOUBLE PRECISION t1
    INTEGER i, j, ierr

    ! Subtract harmonic potential from vcoeff and store coeffs in alpha.
    alpha = 0.d0
    alpha(0:norder_v) = vcoeff(0:norder_v)
    alpha(0) = alpha(0) - 0.25d0*omega
    alpha(2) = alpha(2) - 0.125d0*omega

    ! Get numerical constants to speed up operations.
    log_fact(0:norder) = eval_log_fact( (/ (i, i=0,norder) /) )

    ! Populate Hamiltonian matrix in the basis of harmonic-oscillator
    ! eigenfunctions.
    do i = 0, norder
      do j = i, norder
        hmatrix(i,j) = eval_hmatrix (norder, i, j, omega, alpha, log_fact)
        if (j>i) hmatrix(j,i) = hmatrix(i,j)
      enddo ! j
    enddo ! i

    ! Version with DSYEVD
    ! Diagonalize Hamiltonian.  NB, alpha now used for eigenvalues.
    lapack_lwork = 1
    lapack_liwork = 1
    allocate(lapack_work(lapack_lwork), lapack_iwork(lapack_liwork))
    lapack_lwork = -1
    lapack_liwork = -1
    call dsyevd ('V', 'U', norder+1, hmatrix, norder+1, alpha, &
       &lapack_work, lapack_lwork, lapack_iwork, lapack_liwork, ierr)
    if (ierr/=0) call quit ('LAPACK error '//trim(i2s(ierr))//'.')
    lapack_lwork = nint(lapack_work(1))
    lapack_liwork = lapack_iwork(1)
    deallocate(lapack_work, lapack_iwork)
    allocate(lapack_work(lapack_lwork), lapack_iwork(lapack_liwork))
    call dsyevd ('V', 'U', norder+1, hmatrix, norder+1, alpha, &
       &lapack_work, lapack_lwork, lapack_iwork, lapack_liwork, ierr)
    if (ierr/=0) call quit ('LAPACK error '//trim(i2s(ierr))//'.')
    deallocate(lapack_work, lapack_iwork)
    orbcoeff(0:norder) = hmatrix(0:norder,0)

    ! Return ground-state components.
    e0 = alpha(0)

    ! Normalize coefficients so L2 norm is 1.
    t1 = 1.d0/sqrt(sum(orbcoeff(0:norder)**2))
    orbcoeff(0:norder) = t1*orbcoeff(0:norder)

  END SUBROUTINE get_ground_state


  DOUBLE PRECISION FUNCTION eval_hmatrix (norder, i, j, omega, alpha, &
     &log_fact)
    !------------------------------------------------------------!
    ! Evaluate the matrix element of the Hamiltonian between the !
    ! i-th and j-th eigenfunctions of the harmonic oscillator.   !
    !------------------------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: norder, i, j
    DOUBLE PRECISION, INTENT(in) :: omega, alpha(0:norder), &
       &log_fact(0:norder)
    ! Numerical constants.
    DOUBLE PRECISION, PARAMETER :: log2 = log(2.d0)
    ! Local variables.
    DOUBLE PRECISION t1, t2, log_sqrt_fact_i_fact_j
    INTEGER k
    eval_hmatrix = 0.d0
    log_sqrt_fact_i_fact_j = 0.5d0*(log_fact(i)+log_fact(j))
    do k = i+j, abs(i-j), -2
      if (k>norder) cycle
      t1 = alpha(k)
      if (k==0) t1 = t1 + omega*(dble(i)+0.5d0)
      t2 = log_sqrt_fact_i_fact_j + log_fact(k) + 0.5d0*dble(k)*log2 -&
        &log_fact((j+k-i)/2) - log_fact((k+i-j)/2) - log_fact((i+j-k)/2)
      eval_hmatrix = eval_hmatrix + t1*exp(t2)
    enddo ! k
  END FUNCTION eval_hmatrix


  ELEMENTAL DOUBLE PRECISION FUNCTION eval_log_fact(k)
    !----------------------------------------------!
    ! Returns the logarithm of the factorial of k. !
    !----------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: k
    ! Local variables.
    INTEGER i
    eval_log_fact = 0.d0
    do i = 2, k
      eval_log_fact = eval_log_fact + log(dble(i))
    enddo ! i
  END FUNCTION eval_log_fact


  SUBROUTINE eval_hermite_poly_norm (n, x, h)
    !---------------------------------------------------!
    ! Evaluate "normalized" Hermite polyonials,         !
    !   h_n(x) = H_n(x) / sqrt(2^n n!) ,                !
    ! for n=0:N at x=X, returning the values in h(0:N). !
    !---------------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: n
    DOUBLE PRECISION, INTENT(in) :: x
    DOUBLE PRECISION, INTENT(inout) :: h(0:n)
    ! Numerical constants.
    DOUBLE PRECISION, PARAMETER :: sqrt2 = sqrt(2.d0)
    ! Local variables.
    DOUBLE PRECISION t0, t1, t2, sqrti, sqrti_1, sqrt2_x
    INTEGER i
    if (n<0) return
    t1 = 1.d0
    h(0) = t1
    if (n<1) return
    sqrt2_x = sqrt2*x
    t0 = sqrt2_x
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


  ! String utilities.


  CHARACTER(12) FUNCTION i2s(n)
    !-----------------------------------------------------------------------!
    ! Convert integers to left justified strings that can be printed in the !
    ! middle of a sentence without introducing large amounts of white space.!
    !-----------------------------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: n
    INTEGER i,j
    INTEGER,PARAMETER :: ichar0=ichar('0')
    i2s=''
    i=abs(n)
    do j=len(i2s),1,-1
      i2s(j:j)=achar(ichar0+mod(i,10))
      i=i/10
      if(i==0)exit
    enddo ! j
    if(n<0)then
      i2s='-'//adjustl(i2s(2:12))
    else
      i2s=adjustl(i2s)
    endif ! n<0
  END FUNCTION i2s


  ! Generic utilities.


  SUBROUTINE quit (msg)
    !---------------------!
    ! Quit with an error. !
    !---------------------!
    IMPLICIT NONE
    CHARACTER(*), INTENT(in), OPTIONAL :: msg
    if (present(msg)) then
      write(6,*)'ERROR : '//msg
    else
      write(6,*)'Quitting.'
    endif
    stop
  END SUBROUTINE quit


END PROGRAM anh_quadrature
