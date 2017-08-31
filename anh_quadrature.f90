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
    DOUBLE PRECISION, ALLOCATABLE :: vcoeff_nat(:), vcoeff(:)
    DOUBLE PRECISION ucentre, omega
    ! Variables defining ground state of Hamiltonian.
    INTEGER norder
    DOUBLE PRECISION e0, vratio
    DOUBLE PRECISION, ALLOCATABLE :: orbcoeff(:)
    ! The wave function converges at a given expansion order if its
    ! associated virial ratio is within VIRIAL_TOL of unity.
    INTEGER, PARAMETER :: MAX_NORDER = 200 ! FIXME - dsyev* hang at norder=202
    DOUBLE PRECISION, PARAMETER :: VIRIAL_TOL = 1.d-9
    ! Variables for quadrature grid.
    INTEGER ngrid
    DOUBLE PRECISION, ALLOCATABLE :: xpower_expval(:)
    DOUBLE PRECISION, ALLOCATABLE :: grid_x(:), grid_p(:)
    ! Evaluation of expectation value of v(x).
    DOUBLE PRECISION x, fx, fexpval, dx, xl, xr, psix2
    DOUBLE PRECISION, ALLOCATABLE :: hbasis(:)
    INTEGER, PARAMETER :: NPOINT_BRUTE_FORCE = 10000
    INTEGER, PARAMETER :: MAX_NGRID = 20
    DOUBLE PRECISION, PARAMETER :: PSI2_TOL = 1.d-12
    ! Mist local variables.
    CHARACTER(2048) line
    INTEGER i, iexp, igrid, ierr
    DOUBLE PRECISION t1

    ! Get V coefficients and norder.
    write(6,*) 'Enter coefficients c_k of expansion of V(u) in natural powers,'
    write(6,*)
    write(6,*) '  V(u) = sum_k=0^n c_k u^k       (a.u., in one line):'
    write(6,*)
    read(5,'(a)',iostat=ierr) line
    if (ierr/=0) call quit()
    write(6,*)
    norder_v = -1
    do
      read(line,*,iostat=ierr) (t1, i=0,norder_v+1)
      if(ierr/=0)exit
      norder_v = norder_v+1
    enddo
    if (norder_v<0) call quit ('Could not parse coefficients.')
    if (mod(norder_v,2)/=0) call quit ('Leading term of V(u) must be an even &
      &power.')
    allocate(vcoeff_nat(0:norder_v))
    read(line,*) vcoeff_nat(0:norder_v)
    if (vcoeff_nat(norder_v)<epsilon(1.d0)) &
       &call quit ('Leading term of V(u) must have a positive coefficient.')

    ! Internally, we work in a dimensionless scale where x = sqrt(omega)*u
    ! and energies are e = E/omega, where omega is optimized so as to yield
    ! the best solution at a fixed expansion order.  Also, we shift the
    ! potential self-consistently so that <x> = 0 at this expansion order.
    ucentre = 0.d0
    call obtain_ucentre_omega (norder_v, vcoeff_nat, 20, ucentre, omega)
    write(6,*) 'Centre          = ', ucentre
    write(6,*) 'Omega (a.u.)    = ', omega

    ! Get internal representation of V in normalized Hermite polynomials.
    allocate(vcoeff(0:norder_v))
    call transform_potential (norder_v, vcoeff_nat, ucentre, omega, vcoeff)

    ! Converge trial ground-state wave function with expansion order.
    do norder = max(norder_v,2), MAX_NORDER
      if (allocated(orbcoeff)) deallocate (orbcoeff)
      allocate (orbcoeff(0:norder))
      orbcoeff = 0.d0
      call get_ground_state (norder, norder_v, vcoeff, e0, orbcoeff, vratio)
      if (abs(vratio-1.d0)<VIRIAL_TOL) exit
    enddo ! norder
    if (norder>MAX_NORDER) then
      write(6,*) 'Virial ratio    = ',vratio
      call quit('Failed to converge virial ratio to target accuracy.')
    endif

    ! Solve again to make plot, and report.
    call get_ground_state (norder, norder_v, vcoeff, e0, orbcoeff, vratio, &
       &nplot=10)
    write(6,*) 'Expansion order = '//trim(i2s(norder))
    write(6,*) 'E0 (a.u.)       = ', e0*omega
    write(6,*) 'Virial ratio    = ', vratio

    ! Brute-force integration of target function.
    ! Locate integration range.
    allocate(hbasis(0:norder))
    x = -0.5d0
    do
      call eval_hermite_poly_norm (norder, x, hbasis)
      t1 = sum(orbcoeff(0:norder)*hbasis(0:norder))
      if (exp(-x*x)<PSI2_TOL*t1*t1) exit
      x = 1.5d0*x
    enddo
    xl = x
    x = 0.5d0
    do
      call eval_hermite_poly_norm (norder, x, hbasis)
      t1 = sum(orbcoeff(0:norder)*hbasis(0:norder))
      if (exp(-x*x)<PSI2_TOL*t1*t1) exit
      x = 1.5d0*x
    enddo
    xr = x
    dx = (xr-xl)/dble(NPOINT_BRUTE_FORCE)
    fexpval = 0.d0
    do i = 0, NPOINT_BRUTE_FORCE
      x = xl + (dble(i)/dble(NPOINT_BRUTE_FORCE))*(xr-xl)
      call eval_hermite_poly_norm (norder, x, hbasis)
      t1 = sum(orbcoeff(0:norder)*hbasis(0:norder))
      psix2 = exp(-x*x)*t1*t1
      !fx = cos(x/sqrt(omega))
      !fx = 1.d0 / (1.d0 + x**2/omega)
      fx = exp(-(x/sqrt(omega)-0.5d0)**2)
      fexpval = fexpval + fx*psix2*dx
    enddo ! i
    write(6,*) '<f>             = ', fexpval
    write(6,*) 'Int. range      = [', xl, ':', xr,']'
    write(6,*)
    write(13,*) 0, fexpval
    write(13,*) MAX_NGRID, fexpval
    write(13,'(a)') '&'
    deallocate(hbasis)

    ! Loop over quadrature grid sizes.
    do ngrid = 2, MAX_NGRID
      allocate(xpower_expval(0:2*ngrid-1), grid_x(ngrid), grid_P(ngrid))
      ! Evaluate expectation values of <x^i>.
      do iexp = 0, 2*ngrid-1
        xpower_expval(iexp) = eval_xpower_expval (norder, orbcoeff, iexp)
      enddo ! iexp
      ! Solve for the quadrature grid.
      write(6,*)'Grid size '//trim(i2s(ngrid))//':'
      !call solve_grid_nl2sol (ngrid, xpower_expval, norder, orbcoeff, grid_x, &
      !   &grid_P)
      call solve_grid_newton (ngrid, xpower_expval, norder, orbcoeff, grid_x, &
         &grid_P, ierr)
      if (ierr/=0) then
        ! Report failure.
        write(6,*) '  Grid did not converge.'
      else
        ! Report success.
        write(12,*)grid_x(1)-1.d0, dble(ngrid-2)
        do igrid = 1, ngrid
          write(6,*)'  Grid point #'//trim(i2s(igrid))//': x, P = ', &
             &grid_x(igrid), grid_P(igrid)
          write(12,*) grid_x(igrid)-0.1d0, dble(ngrid-2)
          write(12,*) grid_x(igrid)-0.1d0, dble(ngrid-2) + grid_P(igrid)
          write(12,*) grid_x(igrid)+0.1d0, dble(ngrid-2) + grid_P(igrid)
          write(12,*) grid_x(igrid)+0.1d0, dble(ngrid-2)
        enddo ! igrid
        write(12,*)grid_x(ngrid)+1.d0, dble(ngrid-2)
        write(12,'(a)')'&'
        ! Evaluate expectation value of potential energy using this grid.
        fexpval = 0.d0
        do igrid = 1, ngrid
          x = grid_x(igrid)
          !fx = cos(x/sqrt(omega))
          !fx = 1.d0 / (1.d0 + x**2/omega)
          fx = exp(-(x/sqrt(omega)-0.5d0)**2)
          fexpval = fexpval + grid_P(igrid)*fx
        enddo ! igrid
        write(6,*)'  <f> = ', fexpval
        write(13,*) ngrid, fexpval
      endif
      write(6,*)
      ! Clean up.
      deallocate(xpower_expval, grid_x, grid_P)
    enddo ! ngrid

  END SUBROUTINE main


  SUBROUTINE obtain_ucentre_omega (norder_v, vcoeff_nat, norder, ucentre, omega)
    !----------------------------------------------------------!
    ! Given the natural-polynomial coefficients of a potential !
    ! VCOEFF_NAT(0:NORDER_V), obtain the values of ucentre and !
    ! omega that minimizes the variational energy or virial    !
    ! ratio error for a trial wave function of order NORDER.   !
    !----------------------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: norder_v, norder
    DOUBLE PRECISION, INTENT(in) :: vcoeff_nat(0:norder_v)
    DOUBLE PRECISION, INTENT(inout) :: ucentre, omega
    DOUBLE PRECISION, PARAMETER :: OMEGA_TOL = 1.d-9
    LOGICAL, PARAMETER :: SET_OMEGA_BY_X2 = .false.
    LOGICAL, PARAMETER :: MINIMIZE_E = .true.
    INTEGER, PARAMETER :: MAX_ITER_CENTRE = 10
    DOUBLE PRECISION vcoeff(0:norder_v), orbcoeff(0:norder), e0, vratio, &
       &ucentre_prev, omega_prev, omega_init
    DOUBLE PRECISION xu, xv, xw, fu, fv, fw, x, f
    INTEGER iter
    LOGICAL rejected

    ! Obtain initial guess for omega so that:
    ! * Scaling a potential V(u) by a multiplicative constant results in the
    !   same dimensionless potential v(x).
    ! * omega is the frequency for a harmonic potential.
    omega_init = sqrt(2.d0)*vcoeff_nat(norder_v)**(2.d0/dble(norder_v+2))

    ! Initialize ucentre.
    ucentre = 0.d0

    if (SET_OMEGA_BY_X2) then

      ! Loop over self-consistence cycles.
      omega = omega_init
      do iter = 1, MAX_ITER_CENTRE

        ! Evaluate <x> at this omega.
        call transform_potential (norder_v, vcoeff_nat, ucentre, omega, vcoeff)
        call get_ground_state (norder, norder_v, vcoeff, e0, orbcoeff, vratio)
        ucentre_prev = ucentre
        omega_prev = omega
        ucentre = ucentre_prev - &
           &eval_xpower_expval (norder, orbcoeff, 1)/sqrt(omega_prev)
        omega = 0.5d0 * omega_prev/eval_xpower_expval (norder, orbcoeff, 2)
        if (abs(ucentre-ucentre_prev)<1.d-10.and.abs(omega-omega_prev)<1.d-10) &
           &exit

      enddo

    else ! .not. SET_OMEGA_BY_X2

      ! Loop over self-consistence cycles.
      do iter = 1, MAX_ITER_CENTRE

        ! Obtain initial guess for omega so that:
        ! * Scaling a potential V(u) by a multiplicative constant results in the
        !   same dimensionless potential v(x).
        ! * omega is the frequency for a harmonic potential.
        xv = omega_init
        call transform_potential (norder_v, vcoeff_nat, ucentre, xv, vcoeff)
        call get_ground_state (norder, norder_v, vcoeff, e0, orbcoeff, vratio)
        if (MINIMIZE_E) then
          fv = e0*xv
        else
          fv = abs(vratio-1.d0)
        endif

        ! Perform line minimization for omega.

        ! Bracket to the left.
        xw = 0.d0
        fw = fv-1.d0
        do
          xu = 0.9d0*xv
          call transform_potential (norder_v, vcoeff_nat, ucentre, xu, vcoeff)
          call get_ground_state (norder, norder_v, vcoeff, e0, orbcoeff, vratio)
          if (MINIMIZE_E) then
            fu = e0*xu
          else
            fu = abs(vratio-1.d0)
          endif
          if (fu>fv) exit
          xw = xv
          fw = fv
          xv = xu
          fv = fu
        enddo

        ! Bracket to the right.
        if (fw<fu) then
          do
            xw = 1.1d0*xv
            call transform_potential (norder_v, vcoeff_nat, ucentre, xw, &
               &vcoeff)
            call get_ground_state (norder, norder_v, vcoeff, e0, orbcoeff, &
               &vratio)
            if (MINIMIZE_E) then
              fw = e0*xw
            else
              fw = abs(vratio-1.d0)
            endif
            if (fw>fv) exit
            xu = xv
            fu = fv
            xv = xw
            fv = fw
          enddo
        endif

        ! Zone in on minimum.
        do
          call parabolic_min (xu, xv, xw, fu, fv, fw, x, f, rejected)
          if (rejected) exit
          if (x<=xu.or.x>=xw) exit
          call transform_potential (norder_v, vcoeff_nat, ucentre, x, vcoeff)
          call get_ground_state (norder, norder_v, vcoeff, e0, orbcoeff, &
             &vratio)
          if (MINIMIZE_E) then
            f = e0*x
          else
            f = abs(vratio-1.d0)
          endif
          if (f<fv) then
            if (x<xv) then
              xw = xv
              fw = fv
            elseif (x>xv) then
              xu = xv
              fu = fv
            else
              exit
            endif
            xv = x
            fv = f
          elseif (f>fv) then
            if (x<xv) then
              xu = x
              fu = f
            elseif (x>xv) then
              xw = x
              fw = f
            else
              exit
            endif
          else
            exit
          endif
          if (xw-xu<OMEGA_TOL) exit
        enddo

        ! Return position of minimum.
        omega = xv

        ! Evaluate <x> at this omega and shift ucentre by it.
        call transform_potential (norder_v, vcoeff_nat, ucentre, omega, vcoeff)
        call get_ground_state (norder, norder_v, vcoeff, e0, orbcoeff, vratio)
        ucentre_prev = ucentre
        ucentre = ucentre_prev - &
           &eval_xpower_expval (norder, orbcoeff, 1)/sqrt(omega)
        if (abs(ucentre-ucentre_prev)<1.d-10) exit

      enddo

    endif ! SET_OMEGA_BY_X2 or not

  END SUBROUTINE obtain_ucentre_omega


  SUBROUTINE transform_potential (norder_v, vcoeff_nat, ucentre, omega, vcoeff)
    !---------------------------------------------------------!
    ! Given a function represented as a natural polynomial of !
    ! coefficients VCOEFF_NAT(0:NORDER_V), apply a horizontal !
    ! shift UCENTRE and rescale factor sqrt(OMEGA), and       !
    ! re-represent it in normalized Hermite polynomials of    !
    ! coefficients VCOEFF(0:NORDER_V).                        !
    !---------------------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: norder_v
    DOUBLE PRECISION, INTENT(in) :: vcoeff_nat(0:norder_v), ucentre, omega
    DOUBLE PRECISION, INTENT(inout) :: vcoeff(0:norder_v)
    DOUBLE PRECISION lu_hmatrix(0:norder_v,0:norder_v), &
       &vcoeff_nat1(0:norder_v), vcoeff_nat2(0:norder_v)
    INTEGER i, j, piv_hmatrix(0:norder_v)

    ! Prepare Hermite transformation matrix.
    call lu_decom_hermite_matrix (norder_v, lu_hmatrix, piv_hmatrix)

    ! Apply shift.
    do i = 0, norder_v
      vcoeff_nat1(i) = sum( (/ ( rchoose(i+j,i)*(-ucentre)**j*vcoeff_nat(i+j), &
          &                      j=0,norder_v-i ) /) )
    enddo ! i

    ! Apply change of variable.
    vcoeff_nat2 = (/ ( vcoeff_nat1(i)*omega**(-0.5d0*dble(i+2)), &
       &               i=0,norder_v ) /)

    ! Perform conversion.
    call convert_natpoly_to_hermite (norder_v, lu_hmatrix, piv_hmatrix, &
       &vcoeff_nat2, vcoeff)

  END SUBROUTINE transform_potential


  SUBROUTINE get_ground_state (norder, norder_v, vcoeff, e0, orbcoeff, &
     &vratio, nplot)
    !---------------------------------------------------------!
    ! Given a one-dimensional (anharmonic) potential,         !
    !                                                         !
    !   v(x) = Sum_i VCOEFF(i)*N_i(x) ,                       !
    !                                                         !
    ! where H_i is the i-th Hermite polynomial, construct the !
    ! matrix elements of the Hamiltonian                      !
    !                                                         !
    !   H(x) = -1/2 d/dx^2 + v(x) ,                           !
    !                                                         !
    ! in the basis of the eigenfunctions of the harmonic      !
    ! oscillator of unit frequency,                           !
    !                                                         !
    !   phi_i(x) = exp(-x^2/2) * N_i(x) ,                     !
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
    DOUBLE PRECISION, INTENT(in) :: vcoeff(0:norder_v)
    DOUBLE PRECISION, INTENT(inout) :: e0, orbcoeff(0:norder), vratio
    INTEGER, INTENT(in), OPTIONAL :: nplot
    ! Eigenproblem arrays.
    DOUBLE PRECISION alpha(0:norder), hmatrix(0:norder,0:norder), &
       &cmatrix(0:norder,0:norder)
    ! Buffer for logarithms of factorials.
    DOUBLE PRECISION log_fact(0:norder)
    ! LAPACK work arrays.
    DOUBLE PRECISION, ALLOCATABLE :: lapack_work(:)
    INTEGER, ALLOCATABLE :: lapack_iwork(:)
    INTEGER lapack_lwork, lapack_liwork
    ! Variables for plotting.
    DOUBLE PRECISION, ALLOCATABLE :: hbasis(:)
    DOUBLE PRECISION x, vx, psix
    ! Parameters.
    DOUBLE PRECISION, PARAMETER :: PLOT_ORB_SCALE_FACTOR = 1.d0
    DOUBLE PRECISION, PARAMETER :: PLOT_MAX_X = 5.d0 ! on either side of zero
    INTEGER, PARAMETER :: PLOT_NPOINT = 200 ! on either side of zero
    DOUBLE PRECISION, PARAMETER :: TOL_ZERO = 1.d3*epsilon(1.d0)
    ! Numerical constants.
    DOUBLE PRECISION, PARAMETER :: pi = 4.d0*atan(1.d0)
    DOUBLE PRECISION, PARAMETER :: fourth_root_pi_over_sqrt8 = &
       &                           pi**0.25d0*sqrt(0.125d0)
    ! Virial ratio evaluation.
    DOUBLE PRECISION vircoeff(0:norder_v), xdv_expval, t_expval
    ! Misc local variables.
    INTEGER i, j, k, iplot, ierr
    DOUBLE PRECISION t1, t2

    ! Get numerical constants to speed up operations.
    log_fact(0:norder) = eval_log_fact( (/ (i, i=0,norder) /) )

    ! Populate Hamiltonian matrix in the basis of harmonic-oscillator
    ! eigenfunctions.
    do i = 0, norder
      hmatrix(i,i) = eval_comb_Gamma (norder_v, i, i, vcoeff, log_fact) + &
         &dble(i)+0.25d0 - fourth_root_pi_over_sqrt8*eval_Gamma(i,i,2,log_fact)
      do j = i+1, norder
        hmatrix(i,j) = eval_comb_Gamma (norder_v, i, j, vcoeff, log_fact) - &
           &fourth_root_pi_over_sqrt8*eval_Gamma(i,j,2,log_fact)
        hmatrix(j,i) = hmatrix(i,j)
      enddo ! j
    enddo ! i
    cmatrix = hmatrix

    ! Diagonalize Hamiltonian.
    lapack_lwork = 1
    lapack_liwork = 1
    allocate(lapack_work(lapack_lwork), lapack_iwork(lapack_liwork))
    lapack_lwork = -1
    lapack_liwork = -1
    call dsyevd ('V', 'U', norder+1, cmatrix, norder+1, alpha, &
       &lapack_work, lapack_lwork, lapack_iwork, lapack_liwork, ierr)
    if (ierr/=0) call quit ('DSYEVD error '//trim(i2s(ierr))//'.')
    lapack_lwork = nint(lapack_work(1))
    lapack_liwork = lapack_iwork(1)
    deallocate(lapack_work, lapack_iwork)
    allocate(lapack_work(lapack_lwork), lapack_iwork(lapack_liwork))
    call dsyevd ('V', 'U', norder+1, cmatrix, norder+1, alpha, &
       &lapack_work, lapack_lwork, lapack_iwork, lapack_liwork, ierr)
    if (ierr/=0) call quit ('DSYEVD error '//trim(i2s(ierr))//'.')
    deallocate(lapack_work, lapack_iwork)

    ! Loop over all eigenstates to tidy them up.
    do i = 0, norder
      ! Normalize the wave function.
      t1 = 1.d0/sqrt(sum(cmatrix(0:norder,i)**2))
      cmatrix(0:norder,i) = t1*cmatrix(0:norder,i)
      ! Flush small coefficients to zero.
      where (abs(cmatrix(0:norder,i)) < TOL_ZERO) cmatrix(0:norder,i) = 0.d0
      ! Normalize the wave function again, and make first coefficient positive
      ! while at it.
      t1 = 1.d0/sqrt(sum(cmatrix(0:norder,i)**2))
      do j = 0, norder
        if (abs(cmatrix(j,i)) > TOL_ZERO) then
          if (cmatrix(j,i)<0.d0) t1 = -t1
          exit
        endif
      enddo ! j
      cmatrix(0:norder,i) = t1*cmatrix(0:norder,i)
      ! Recalculate the variational energy.
      alpha(i) = 0.d0
      do j = 0, norder
        t1 = 0.d0
        do k = j+1, norder
          t1 = t1 + cmatrix(k,i)*hmatrix(k,j)
        enddo ! k
        t1 = cmatrix(j,i)*(cmatrix(j,i)*hmatrix(j,j) + 2.d0*t1)
        alpha(i) = alpha(i)+t1
      enddo ! j
    enddo ! i

    ! Evaluate the virial ratio for the ground state.
    vircoeff = 0.d0
    do i = 0, norder_v
      if (i<norder_v-1) vircoeff(i) = vircoeff(i) + &
         &vcoeff(i+2)*sqrt(dble((i+1)*(i+2)))
      if (i>1) vircoeff(i) = vircoeff(i) + vcoeff(i)*dble(i)
    enddo ! i
    xdv_expval = 0.d0
    t_expval = 0.d0
    do i = 0, norder
      t1 = eval_comb_Gamma (norder_v, i, i, vircoeff, log_fact)
      t2 = dble(i)+0.25d0 - fourth_root_pi_over_sqrt8*eval_Gamma(i,i,2,log_fact)
      xdv_expval = xdv_expval + t1*cmatrix(i,0)**2
      t_expval = t_expval + t2*cmatrix(i,0)**2
      do j = i+1, norder
        t1 = eval_comb_Gamma (norder_v, i, j, vircoeff, log_fact)
        t2 = -fourth_root_pi_over_sqrt8*eval_Gamma(i,j,2,log_fact)
        xdv_expval = xdv_expval + 2.d0*t1*cmatrix(i,0)*cmatrix(j,0)
        t_expval = t_expval + 2.d0*t2*cmatrix(i,0)*cmatrix(j,0)
      enddo ! j
    enddo ! i
    vratio = 2.d0*t_expval/xdv_expval

    ! Return ground-state components.
    ! FIXME - if ground state is degenerate, how do we choose which
    ! eigenstate to return?
    e0 = alpha(0)
    orbcoeff(0:norder) = cmatrix(0:norder,0)

    ! Optionally, plot potential and first NPLOT eigenstates in xmgrace format
    ! to Fortran unit 11.
    if (present(nplot)) then
      allocate (hbasis(0:norder))
      ! Plot V(x).
      do i = -PLOT_NPOINT, PLOT_NPOINT
        x = dble(i)/dble(PLOT_NPOINT) * PLOT_MAX_X
        call eval_hermite_poly_norm (norder_v, x, hbasis)
        vx = sum(vcoeff(0:norder_v)*hbasis(0:norder_v))
        write(11,*) x, vx
      enddo ! i
      ! Plot Psi_n(x).
      do iplot = 0, min(norder,nplot)
        write(11,'(a)') '&'
        write(11,*) -PLOT_MAX_X, alpha(iplot)
        write(11,*) PLOT_MAX_X, alpha(iplot)
        write(11,'(a)') '&'
        do i = -PLOT_NPOINT, PLOT_NPOINT
          x = dble(i)/dble(PLOT_NPOINT) * PLOT_MAX_X
          call eval_hermite_poly_norm (norder, x, hbasis)
          t1 = exp(-0.5d0*x*x)
          psix = t1*sum(cmatrix(0:norder,iplot)*hbasis(0:norder))
          write(11,*) x, alpha(iplot)+PLOT_ORB_SCALE_FACTOR*psix
        enddo ! i
      enddo ! iplot
      deallocate (hbasis)
    endif

  END SUBROUTINE get_ground_state


  DOUBLE PRECISION FUNCTION eval_xpower_expval (norder, orbcoeff, iexp)
    !---------------------------------------------------------------------!
    ! Given the coefficients ORBCOEFF(0:NORDER) of a trial wave function, !
    ! evaluate the expectation value of the natural power <x^IEXP> .      !
    !---------------------------------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: norder, iexp
    DOUBLE PRECISION, INTENT(in) :: orbcoeff(0:norder)
    ! Misc local variables.
    INTEGER i, j, piv_hmatrix(0:iexp)
    DOUBLE PRECISION lu_hmatrix(0:iexp,0:iexp), xcoeff_nat(0:iexp), &
       &xcoeff(0:iexp), log_fact(0:max(iexp,norder))

    ! Re-represent x^iexp in the basis of Hermite polynomials.
    xcoeff_nat(0:iexp) = 0.d0
    xcoeff_nat(iexp) = 1.d0
    call lu_decom_hermite_matrix (iexp, lu_hmatrix, piv_hmatrix)
    call convert_natpoly_to_hermite (iexp, lu_hmatrix, piv_hmatrix, &
       &xcoeff_nat, xcoeff)

    ! Get numerical constants to speed up operations.
    log_fact(0:max(iexp,norder)) = &
       &eval_log_fact( (/ (i, i=0,max(iexp,norder)) /) )

    ! Evaluate expectation value.
    eval_xpower_expval = 0.d0
    do i = 0, norder
      eval_xpower_expval = eval_xpower_expval + orbcoeff(i) * orbcoeff(i) * &
         &eval_comb_Gamma (iexp, i, i, xcoeff, log_fact)
      do j = i+1, norder
        eval_xpower_expval = eval_xpower_expval + &
           &2.d0 * orbcoeff(j) * orbcoeff(i) * eval_comb_Gamma &
           &(iexp, i, j, xcoeff, log_fact)
      enddo ! j
    enddo ! i
    if (abs(eval_xpower_expval)<1.d3*epsilon(1.d0)) eval_xpower_expval = 0.d0

  END FUNCTION eval_xpower_expval


  SUBROUTINE solve_grid_newton (ngrid, xpower_expval, norder, orbcoeff, &
     &grid_x, grid_P, ierr)
    !------------------------------------------------------------!
    ! Given <x^n> for n=0:2*NGRID-1, solve for the parameters of !
    ! the NGRID-point quadrature grid using Newton's method.     !
    !------------------------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: ngrid, norder
    DOUBLE PRECISION, INTENT(in) :: xpower_expval(0:2*ngrid-1), &
       &orbcoeff(0:norder)
    DOUBLE PRECISION, INTENT(inout) :: grid_x(ngrid), grid_P(ngrid)
    INTEGER, INTENT(inout) :: ierr
    ! Maximum number of Nweton's method iterations to attempt.
    INTEGER, PARAMETER :: MAX_ITER = 50
    LOGICAL, PARAMETER :: VERBOSE = .false.
    ! Misc local variables.
    DOUBLE PRECISION fvec(2*ngrid), Jmat(2*ngrid,2*ngrid), t1, &
       &grid_P_test(ngrid), grid_x_test(ngrid), xscale, x0max, x, &
       &hbasis(0:norder), dorbcoeff(0:norder-1), grid_orb(ngrid), &
       &xpower_renorm(0:2*ngrid-1), xu, xv, xw, fu, fv, fw, f, f0, &
       &fvec_test(2*ngrid)
    INTEGER i, iter, piv(2*ngrid)
    LOGICAL rejected

    ierr = 0

    ! Define x rescaling factor.
    xscale = xpower_expval(2*ngrid-2)**(1.d0/dble(2*ngrid-2))
    xpower_renorm = (/ ( xpower_expval(i)/xscale**dble(i), i=0,2*ngrid-1 ) /)
    x0max = 1.d0 + 0.423d0*log(dble(ngrid-1))
    if (VERBOSE) then
      write(6,*) '  Using xscale = ', xscale
      write(6,*) '  Using x0max  = ', x0max
    endif

    ! Define wave function derivative.
    do i = 0, norder-1
      dorbcoeff(i) = orbcoeff(i+1) * sqrt(dble(2*i+2))
    enddo ! i

    ! Initialize to something sensible.
    grid_x(1:ngrid) = (/ ( xpower_renorm(1) - x0max + &
       &                   2*x0max*dble(i-1)/dble(ngrid-1), i=1,ngrid) /)
    grid_P = 1.d0

    ! Loop over Newton's method iterations.
    do iter = 1, MAX_ITER

      ! Obtain Newton's method step.
      call eval_fvec_Jmat (ngrid, grid_x, grid_P, xscale, xpower_renorm, &
         &norder, orbcoeff, dorbcoeff, fvec, Jmat)
      f0 = sum(fvec*fvec)
      call dgetrf (2*ngrid, 2*ngrid, Jmat, 2*ngrid, piv, ierr)
      if (ierr/=0) return
      call dgetrs ('N', 2*ngrid, 1, Jmat, 2*ngrid, piv, fvec, 2*ngrid, ierr)
      if (ierr/=0) return

      ! Test Newton's step.
      xu = 0.d0
      fu = f0
      xw = 1.d0
      do
        grid_P_test = grid_P - xw*fvec(ngrid+1:2*ngrid)
        grid_x_test = grid_x - xw*fvec(1:ngrid)
        if ( all(grid_P_test>0.d0) .and. &
           & all(grid_x_test(1:ngrid-1)<grid_x_test(2:ngrid)) ) then
          call eval_fvec (ngrid, grid_x_test, grid_P_test, xscale, &
             &xpower_renorm, norder, orbcoeff, fvec_test)
          fw = sum(fvec_test*fvec_test)
          exit
        endif
        xw = 0.8d0*xw
      enddo

      if (fw>fu) then
        ! Find middle point.
        xv = 0.5d0*xw
        do
          grid_P_test = grid_P - xv*fvec(ngrid+1:2*ngrid)
          grid_x_test = grid_x - xv*fvec(1:ngrid)
          call eval_fvec (ngrid, grid_x_test, grid_P_test, xscale, &
             &xpower_renorm, norder, orbcoeff, fvec_test)
          fw = sum(fvec_test*fvec_test)
          if (fv<=fu) exit
          xw = xv
          fw = fv
          xv = 0.5d0*xw
        enddo
      elseif (fw<fu) then
        ! Find right-side point.
        xv = xw
        fv = fw
        xw = 2.d0*xv
        do
          grid_P_test = grid_P - xw*fvec(ngrid+1:2*ngrid)
          grid_x_test = grid_x - xw*fvec(1:ngrid)
          if ( any(grid_P_test<=0.d0) .or. &
             & any(grid_x_test(1:ngrid-1)>=grid_x_test(2:ngrid)) ) then
            xw = xv + 0.9d0*(xw-xv)
            cycle
          endif
          call eval_fvec (ngrid, grid_x_test, grid_P_test, xscale, &
             &xpower_renorm, norder, orbcoeff, fvec_test)
          fw = sum(fvec_test*fvec_test)
          if (fv<=fw) exit
          xv = xw
          fv = fw
          xw = xv + 2.d0*(xw-xv)
        enddo
      else
        xv = 0.d0
        fv = fu
      endif

      ! Zone in on minimum.
      do
        call parabolic_min (xu, xv, xw, fu, fv, fw, x, f, rejected)
        if (rejected) exit
        if (x<=xu.or.x>=xw) exit
        grid_P_test = grid_P - x*fvec(ngrid+1:2*ngrid)
        grid_x_test = grid_x - x*fvec(1:ngrid)
        call eval_fvec (ngrid, grid_x_test, grid_P_test, xscale, &
           &xpower_renorm, norder, orbcoeff, fvec_test)
        f = sum(fvec_test*fvec_test)
        if (f<fv) then
          if (x<xv) then
            xw = xv
            fw = fv
          elseif (x>xv) then
            xu = xv
            fu = fv
          else
            exit
          endif
          xv = x
          fv = f
        elseif (f>fv) then
          if (x<xv) then
            xu = x
            fu = f
          elseif (x>xv) then
            xw = x
            fw = f
          else
            exit
          endif
        else
          exit
        endif
        if (xw-xu<1.d-3) exit
      enddo

      ! Apply step.
      fvec = xv*fvec
      grid_x = grid_x - fvec(1:ngrid)
      grid_P = grid_P - fvec(ngrid+1:2*ngrid)

      ! Check function and parameter convergence.
      t1 = sqrt(sum(fvec**2)/dble(2*ngrid))
      if (VERBOSE) then
        write(6,*) '  Iteration '//trim(i2s(iter))//':'
        write(6,*) '    delta f    = ', f0-fv
        write(6,*) '    Step scale = ', xv
        write(6,*) '    Disp norm  = ', t1
        write(6,*) '    Unscaled parameters:'
        do i = 1, ngrid
          write(6,*) '      x,P #'//trim(i2s(i))//': ', grid_x(i), grid_P(i)
        enddo ! i
      endif ! VERBOSE
      if (abs(f0)<1.d-12) exit
      if (abs(f0-fv)<1.d-8*abs(f0)) then
        ierr = 1
        exit
      endif
      if (t1<1.d-10) then
        ierr = 2
        exit
      endif

    enddo
    if (iter>MAX_ITER) ierr = 3

    ! Rescale grid parameters.
    grid_x = grid_x*xscale
    do i = 1, ngrid
      x = grid_x(i)
      call eval_hermite_poly_norm (norder, x, hbasis)
      grid_orb(i) = exp(-0.5d0*x*x) * sum(orbcoeff(0:norder)*hbasis(0:norder))
    enddo ! i
    grid_P = grid_P*grid_orb**2
    grid_P = grid_P/sum(grid_P)

  END SUBROUTINE solve_grid_newton


  SUBROUTINE eval_fvec (ngrid, grid_x, grid_P, xscale, xpower_renorm, norder, &
     &orbcoeff, fvec)
    !-----------------------------------------------!
    ! Evaluate target function for Newton's method. !
    !-----------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: ngrid, norder
    DOUBLE PRECISION, INTENT(in) :: grid_x(ngrid), grid_P(ngrid), xscale, &
       &xpower_renorm(0:2*ngrid-1), orbcoeff(0:norder)
    DOUBLE PRECISION, INTENT(inout) :: fvec(2*ngrid)
    DOUBLE PRECISION grid_orb(ngrid), x, hbasis(0:norder)
    INTEGER i

    ! Evaluate required objects.
    do i = 1, ngrid
      x = grid_x(i)*xscale
      call eval_hermite_poly_norm (norder, x, hbasis)
      grid_orb(i) = exp(-0.5d0*x*x) * sum(orbcoeff(0:norder)*hbasis(0:norder))
    enddo ! i

    ! Evaluate target function.
    fvec(1) = sum( grid_P(1:ngrid)*grid_orb(1:ngrid)**2 ) - 1.d0
    fvec(2:2*ngrid) = (/ ( sum( grid_P(1:ngrid)*grid_orb(1:ngrid)**2*&
       &                        ( grid_x(1:ngrid)**dble(i-1) - &
       &                          xpower_renorm(i-1) ) ), i=2,2*ngrid ) /)

  END SUBROUTINE eval_fvec


  SUBROUTINE eval_fvec_Jmat (ngrid, grid_x, grid_P, xscale, xpower_renorm, &
     &norder, orbcoeff, dorbcoeff, fvec, Jmat)
    !------------------------------------------------------------!
    ! Evaluate target function and Jacobian for Newton's method. !
    !------------------------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: ngrid, norder
    DOUBLE PRECISION, INTENT(in) :: grid_x(ngrid), grid_P(ngrid), xscale, &
       &xpower_renorm(0:2*ngrid-1), orbcoeff(0:norder), dorbcoeff(0:norder-1)
    DOUBLE PRECISION, INTENT(inout) :: fvec(2*ngrid), Jmat(2*ngrid,2*ngrid)
    DOUBLE PRECISION grid_orb(ngrid), grid_dorb(ngrid), x, hbasis(0:norder)
    INTEGER i

    ! Evaluate required objects.
    do i = 1, ngrid
      x = grid_x(i)*xscale
      call eval_hermite_poly_norm (norder, x, hbasis)
      grid_orb(i) = exp(-0.5d0*x*x) * sum(orbcoeff(0:norder)*hbasis(0:norder))
      grid_dorb(i) = exp(-0.5d0*x*x) * &
         &sum(dorbcoeff(0:norder-1)*hbasis(0:norder-1)) - x*grid_orb(i)
    enddo ! i

    ! Evaluate target function.
    fvec(1) = sum( grid_P(1:ngrid)*grid_orb(1:ngrid)**2 ) - 1.d0
    fvec(2:2*ngrid) = (/ ( sum( grid_P(1:ngrid)*grid_orb(1:ngrid)**2*&
       &                        ( grid_x(1:ngrid)**dble(i-1) - &
       &                          xpower_renorm(i-1) ) ), i=2,2*ngrid ) /)

    ! Evaluate Jacobian.
    Jmat = 0.d0
    Jmat(1,1:ngrid) = grid_P(1:ngrid)*2.d0*grid_orb(1:ngrid)*&
       &grid_dorb(1:ngrid)*xscale
    Jmat(1,ngrid+1:2*ngrid) = grid_orb(1:ngrid)**2
    do i = 2, 2*ngrid
      Jmat(i,1:ngrid) = grid_P(1:ngrid) * &
         &( dble(i-1)*grid_orb(1:ngrid)**2*grid_x(1:ngrid)**dble(i-2) + &
         &  2.d0*grid_orb(1:ngrid)*grid_dorb(1:ngrid)*xscale* &
         &  ( grid_x(1:ngrid)**dble(i-1) - xpower_renorm(i-1) ) )
      Jmat(i,ngrid+1:2*ngrid) = grid_orb(1:ngrid)**2 * &
         &( grid_x(1:ngrid)**dble(i-1) - xpower_renorm(i-1) )
    enddo ! i

  END SUBROUTINE eval_fvec_Jmat


  ! Hermite polynomial tools.  NB, we use "normalized" Hermite polynomials,
  ! N_n(x) = pi^-1/4 2^-n/2 (n!)^-1/2 H_n(x).


  DOUBLE PRECISION FUNCTION eval_Gamma (i, j, k, log_fact)
    !--------------------------------------------------------------!
    ! Evaluate the Gamma_i,j,k symbol,                             !
    !                                                              !
    !   Gamma_i,j,k = integral exp(-x^2) N_i(x) N_j(x) N_k(x) dx . !
    !--------------------------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: i, j, k
    DOUBLE PRECISION, INTENT(in) :: log_fact(0:)
    ! Numerical constants.
    DOUBLE PRECISION, PARAMETER :: pi = 4.d0*atan(1.d0)
    DOUBLE PRECISION, PARAMETER :: fourth_log_pi = 0.25d0*log(pi)
    ! Local variables.
    DOUBLE PRECISION t1
    eval_Gamma = 0.d0
    if (k<abs(i-j) .or. k>i+j .or. mod(i+j+k,2)/=0) return
    t1 = -fourth_log_pi + 0.5d0*(log_fact(i) + log_fact(j) + log_fact(k)) -&
      &log_fact((j+k-i)/2) - log_fact((k+i-j)/2) - log_fact((i+j-k)/2)
    eval_Gamma = exp(t1)
  END FUNCTION eval_Gamma


  DOUBLE PRECISION FUNCTION eval_comb_Gamma (n, i, j, coeff, log_fact)
    !--------------------------------------------------!
    ! Evaluate the linear combination of Gamma symbols !
    !                                                  !
    !   Sum_k=0^n coeff_k * Gamma_i,j,k .              !
    !                                                  !
    ! NB, Gamma evaluator inlined to skip the logic.   !
    !--------------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: n, i, j
    DOUBLE PRECISION, INTENT(in) :: coeff(0:n), log_fact(0:)
    ! Numerical constants.
    DOUBLE PRECISION, PARAMETER :: pi = 4.d0*atan(1.d0)
    DOUBLE PRECISION, PARAMETER :: fourth_log_pi = 0.25d0*log(pi)
    ! Local variables.
    DOUBLE PRECISION t1
    INTEGER k
    eval_comb_Gamma = 0.d0
    do k = i+j, abs(i-j), -2
      if (k>n) cycle
      t1 = -fourth_log_pi + 0.5d0*(log_fact(i) + log_fact(j) + log_fact(k)) -&
        &log_fact((j+k-i)/2) - log_fact((k+i-j)/2) - log_fact((i+j-k)/2)
      eval_comb_Gamma = eval_comb_Gamma + coeff(k)*exp(t1)
    enddo ! k
  END FUNCTION eval_comb_Gamma


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


  SUBROUTINE lu_decom_hermite_matrix (norder, lu_hmatrix, piv_hmatrix)
    !-----------------------------------------------------------------!
    ! Returns the LU decomposition of the matrix of coefficients of   !
    ! natural powers in Hermite polynomials up to order NORDER, which !
    ! can then be used to express natural polynomials in Hermite      !
    ! polynomials.                                                    !
    !-----------------------------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: norder
    DOUBLE PRECISION, INTENT(inout) :: lu_hmatrix(0:norder,0:norder)
    INTEGER, INTENT(inout) :: piv_hmatrix(0:norder)
    ! Numerical constants.
    DOUBLE PRECISION, PARAMETER :: pi = 4.d0*atan(1.d0)
    DOUBLE PRECISION, PARAMETER :: fourth_log_pi = 0.25d0*log(pi)
    DOUBLE PRECISION, PARAMETER :: log2 = log(2.d0)
    ! Misc local variables.
    DOUBLE PRECISION t1
    INTEGER n, i, isgn, ierr

    ! Build matrix of coefficients of N_n as a linear combination of natural
    ! powers.
    lu_hmatrix = 0.d0
    do n = 0, norder
      isgn = 1
      do i = n, 0, -2
        t1 = 0.5d0*(eval_log_fact(n) + dble(2*i-n)*log2) - &
           &eval_log_fact(i) - eval_log_fact(n/2-i/2) - fourth_log_pi
        lu_hmatrix(i,n) = dble(isgn) * exp(t1)
        isgn = -isgn
      enddo ! i
    enddo ! n

    ! LU-decompose the matrix.
    call dgetrf (norder+1, norder+1, lu_hmatrix, norder+1, piv_hmatrix, ierr)
    if (ierr/=0) call quit ('DGETRF error '//trim(i2s(ierr))//'.')

  END SUBROUTINE lu_decom_hermite_matrix


  SUBROUTINE convert_natpoly_to_hermite (norder, lu_hmatrix, piv_hmatrix, &
     &pcoeff, hcoeff)
    !----------------------------------------------------------!
    ! Convert a polynomial of coefficients PCOEFF(0:NORDER) to !
    ! a Hermite polynomial of coefficients HCOEFF(0:NORDER).   !
    !----------------------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: norder, piv_hmatrix(0:norder)
    DOUBLE PRECISION, INTENT(in) :: lu_hmatrix(0:norder,0:norder), &
       &pcoeff(0:norder)
    DOUBLE PRECISION, INTENT(inout) :: hcoeff(0:norder)
    INTEGER ierr
    hcoeff(0:norder) = pcoeff(0:norder)
    call dgetrs ('N', norder+1, 1, lu_hmatrix, norder+1, piv_hmatrix, hcoeff, &
       &norder+1, ierr)
    if (ierr/=0) call quit ('DGETRS error '//trim(i2s(ierr))//'.')
  END SUBROUTINE convert_natpoly_to_hermite


  ! Numerical tools.


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


  SUBROUTINE parabolic_min (x1, x2, x3, y1, y2, y3, x0, y0, rejected)
    !-----------------------------------------------------------------!
    ! Fit three points to a parabola and return (x,y) of the min/max. !
    !-----------------------------------------------------------------!
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(in) :: x1, x2, x3, y1, y2, y3
    DOUBLE PRECISION, INTENT(out) :: x0, y0
    LOGICAL, INTENT(out) :: rejected
    DOUBLE PRECISION a, b, c, numa, numb, numc, den, x1_sq, x2_sq, x3_sq,&
       &x21, x32, x31, z1, z2, z3, invden

    ! Initialize.
    x0 = x2
    y0 = y2
    rejected = .false.

    ! Check that x and y values are distinguishable and in correct order.
    if (x1>=x2 .or. x2>=x3 .or. (y2>=y1.eqv.y3>=y2))then
      rejected = .true.
      return
    endif

    ! Compute squares.
    x1_sq = x1*x1
    x2_sq = x2*x2
    x3_sq = x3*x3

    ! Renormalize for better numerics.
    x31 = x3-x1
    x21 = (x2-x1)/x31
    x32 = (x3-x2)/x31
    z1 = y1*x32
    z2 = y2
    z3 = y3*x21

    ! Solve linear system.
    den = -x1_sq*x32 + x2_sq - x3_sq*x21
    numa = -z1 + z2 - z3
    numb = z1*(x2+x3) - z2*(x1+x3) + z3*(x1+x2)
    numc = -z1*x2*x3 + z2*x3*x1 - z3*x1*x2

    ! Find x0 and y0.
    invden = 1.d0/den
    a = numa*invden
    b = numb*invden
    c = numc*invden
    x0 = -0.5d0*numb/numa
    y0 = (a*x0+b)*x0 + c

  END SUBROUTINE parabolic_min


  DOUBLE PRECISION FUNCTION rchoose(a,b)
    !-------------------------------------!
    ! This function returns a choose b as !
    ! a floating-point real number.       !
    !-------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: a,b
    INTEGER i
    rchoose=1.d0
    do i=1,b
      rchoose=rchoose*(dble(a+1-i)/dble(i))
    enddo ! i
  END FUNCTION rchoose


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
