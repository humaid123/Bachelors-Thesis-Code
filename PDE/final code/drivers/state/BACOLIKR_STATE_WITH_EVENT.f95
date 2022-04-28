module rootfinding
    
    use bacoli95_mod, only: bacoli95_sol, bacoli95_vals

    implicit none

!   Global variables and data structures
    integer, parameter          :: dp = kind(0d0)
    type(bacoli95_sol), pointer :: sol      ! pointer to sol object
    real(dp), pointer           :: work(:)  ! work array

    contains

        subroutine rt1(t, x, nint, ub, utb, neq, rval, nrt)
            ! Input arguments
            integer ::  nint, neq, nrt 
            real(dp) :: t, x(nint+1), ub(neq), utb(neq)
           
            ! Output arguments
            real(dp) :: rval(nrt)
    
            ! Local variables
            real(dp)            :: approx   ! Integral approximation
            real(dp)            :: int      ! Target value of the integral

            real(dp)            :: u(sol%npde, nint*(sol%kcol+3), 1) !u(6500)   ! Solution values at x values
                                            ! Breaks if nint_max > 500
            real(dp)            :: h(500)   ! h values; assumes nint_max <= 500
            real(dp)            :: rho(13)  ! Gauss points per subinterval
            real(dp)            :: wts(13)  ! Gauss weights per subinterval
            real(dp)            :: excol(nint*(sol%kcol+3)), ewts(nint*(sol%kcol+3)) ! excol(6500), ewts(6500) 
                                ! Total set of Gauss points and weights over 
                                ! entire domain; assumes nint_max <= 500.
                    
            integer             :: npts     ! Number of evaluation points
            integer, parameter  :: nderiv = 0 
                                   ! Number of spatial derivatives of
                                   ! solution that are required
            integer             :: i, ii, j        ! loop indexes

            real(dp) xplot(1), yplot(sol%npde, 1)
            integer, parameter :: nplot = 1
            integer n_times_called
            common /n_times_called/ n_times_called

            n_times_called = n_times_called + 1
            !-----------------------------------------------------------

            ! For efficiency, the computation of the Gauss pts and Wts should
            ! be done once beforehand and then the array can be available
            ! globally. For now we leave these computations here.
            ! For the kcol+3 Gaussian Quadrature rule approximation to the
            ! integral we will need to evaluate the approximate solution
            ! at nint*(kcol+3) points, i.e., the kcol+3 Gauss points on 
            ! each subinterval. Therefore
            npts = nint*(sol%kcol+3)
            
            ! Next we compute the entire set of Gauss points and weights 
            ! for the entire spatial domain. The code below is taken from bacoli.f.
            
            ! Calculate the subinterval size sequence.
            do i = 1, nint
                h(i) = x(i+1)-x(i)
            end do
    
            ! Compute the kcol+3 Gaussian points and Gaussian weights on [0,1].
            ! (In the call below to the F77 routine, gauleg, the work array
            ! needs to be of dimension at least (kcol+3)**2, which has a max
            ! of 169 (since the max kcol value is 10). The size of work is
            ! set in the next routine below, called setsol. As long as nint_max
            ! is large enough, work will be large enough for the call to gauleg.)
            ! The points and weights are returned in rho and wts resepctively.

            ! HUMAID - this routine is defined in bacoli-aux 
            ! HUMAID - it returns the GAUSS-LEGENDRE points 
            ! HUMAID - the 4 means it is defined from [0, 1]
             call gauleg(sol%kcol+3, (sol%kcol+3)*(sol%kcol+3), rho, wts, work, 4)

            ! HUMAID - this scales the gauss points over the whole interval.
            ! HUMAID - rho are the points, wts are the weights
            ! HUMAID the scaled points and weights are excol and ewts
            do i = 1, nint
                ii = (i - 1) * (sol%kcol + 3)
                do j = 1, sol%kcol+3
                    excol(ii + j) = x(i) + h(i) * rho(j)
                    ewts(ii + j) = h(i) * wts(j)
                end do
            end do

            ! The 'values' routine is an F77 routine that can be
            ! accessed globally.

            call values(sol%kcol, excol, nint, x, sol%npde, npts, nderiv,  &
                                                             u, ub, work)
    
            ! Apply the Gauss rule to obtain the integral approximation
            ! integral = sigma (wts * points)
            approx = 0.0d0
            do i = 1, nint
                ii = (i - 1) * (sol%kcol + 3)
                do j = 1, sol%kcol+3
                    approx = approx + ewts(ii+j)*u(2, ii+j, 1)
                end do
            end do

            ! to plot:
            xplot(1) = 0
            call values(sol%kcol, xplot, nint, x, sol%npde, nplot, nderiv,  &
                    yplot, ub, work)
    
            ! write(*, *) "plotting", xplot(1), yplot(2, 1)

            ! target value for the root
            ! old => int = 5d0
            int = 30000d0

            ! We next assign the gstop value to be difference between the
            ! approximate integral value and the target integral value. 
            rval(1) = approx - int

            ! we plot from the root function and clean up in Python.
            write(11, *) t,  yplot(2, 1), approx
    
        end subroutine rt1

        subroutine rt2(t, x, nint, ub, utb, neq, rval, nrt)
            ! Input arguments
            integer ::  nint, neq, nrt 
            real(dp) :: t, x(nint+1), ub(neq), utb(neq)
           
            ! Output arguments
            real(dp) :: rval(nrt)
    
            ! Local variables
            real(dp)            :: approx   ! Integral approximation
            real(dp)            :: int      ! Target value of the integral

            real(dp)            :: u(sol%npde, nint*(sol%kcol+3), 1) !u(6500)   ! Solution values at x values
                                            ! Breaks if nint_max > 500
            real(dp)            :: h(500)   ! h values; assumes nint_max <= 500
            real(dp)            :: rho(13)  ! Gauss points per subinterval
            real(dp)            :: wts(13)  ! Gauss weights per subinterval
            real(dp)            :: excol(nint*(sol%kcol+3)), ewts(nint*(sol%kcol+3)) ! excol(6500), ewts(6500) 
                                ! Total set of Gauss points and weights over 
                                ! entire domain; assumes nint_max <= 500.
                    
            integer             :: npts     ! Number of evaluation points
            integer, parameter  :: nderiv = 0 
                                   ! Number of spatial derivatives of
                                   ! solution that are required
            integer             :: i, ii, j        ! loop indexes

            real(dp) xplot(1), yplot(sol%npde, 1)
            integer, parameter :: nplot = 1
            integer n_times_called
            common /n_times_called/ n_times_called

            n_times_called = n_times_called + 1
            !-----------------------------------------------------------

            ! For efficiency, the computation of the Gauss pts and Wts should
            ! be done once beforehand and then the array can be available
            ! globally. For now we leave these computations here.
            ! For the kcol+3 Gaussian Quadrature rule approximation to the
            ! integral we will need to evaluate the approximate solution
            ! at nint*(kcol+3) points, i.e., the kcol+3 Gauss points on 
            ! each subinterval. Therefore
            npts = nint*(sol%kcol+3)
            
            ! Next we compute the entire set of Gauss points and weights 
            ! for the entire spatial domain. The code below is taken from bacoli.f.
            
            ! Calculate the subinterval size sequence.
            do i = 1, nint
                h(i) = x(i+1)-x(i)
            end do
    
            ! Compute the kcol+3 Gaussian points and Gaussian weights on [0,1].
            ! (In the call below to the F77 routine, gauleg, the work array
            ! needs to be of dimension at least (kcol+3)**2, which has a max
            ! of 169 (since the max kcol value is 10). The size of work is
            ! set in the next routine below, called setsol. As long as nint_max
            ! is large enough, work will be large enough for the call to gauleg.)
            ! The points and weights are returned in rho and wts resepctively.

            ! this routine is defined in bacali-aux 
            ! it returns the GAUSS-LEGENDRE points 
            ! the 4 means it is defined from [0, 1]
             call gauleg(sol%kcol+3, (sol%kcol+3)*(sol%kcol+3), rho, wts, work, 4)

            ! this scales the gauss points over the whole interval.
            ! rho are the points, wts are the weights
            ! the scaled points and weights are returned in excol and ewts
            do i = 1, nint
                ii = (i - 1) * (sol%kcol + 3)
                do j = 1, sol%kcol+3
                    excol(ii + j) = x(i) + h(i) * rho(j)
                    ewts(ii + j) = h(i) * wts(j)
                end do
            end do

            ! uses the BPSLINE coefficients to get the values at excol and 
            ! places them in u
            call values(sol%kcol, excol, nint, x, sol%npde, npts, nderiv,  &
                                                             u, ub, work)
    
            ! Apply the Gauss rule to obtain the integral approximation
            ! integral = sigma (wts * points)
            approx = 0.0d0
            do i = 1, nint
                ii = (i - 1) * (sol%kcol + 3)
                do j = 1, sol%kcol+3
                    approx = approx + ewts(ii+j)*u(2, ii+j, 1)
                end do
            end do

            ! to plot:
            xplot(1) = 0
            call values(sol%kcol, xplot, nint, x, sol%npde, nplot, nderiv,  &
                    yplot, ub, work)
    
            ! target value for the root
            ! old int = 1d0
            int = 10000d0

            ! We next assign the gstop value to be difference between the
            ! approximate integral value and the target integral value. 
            rval(1) = approx - int

            ! we plot from the root function and clean up in Python
            write(11, *) t,  yplot(2, 1), approx
    
        end subroutine rt2


    ! allocate data and keep a pointer to the sol
    subroutine setsol(in_sol)
        type(bacoli95_sol), target :: in_sol

        ! Local variables used in defining size of work array
        integer, parameter :: nconti = 2
        integer, parameter :: nderiv = 0

        ! Assign the pointer in the module scope to be equal to in_sol
        sol => in_sol
        
        ! Allocate work memory for use in 'rt'
        allocate(work((sol%kcol+nconti)*(nderiv+1)+(sol%kcol+nconti+1)   &
           *(sol%kcol+nconti+2)/2+2*sol%kcol*(sol%nint_max+1)+2*nconti))

    end subroutine setsol

end module rootfinding



program trimesh_example_driver

    use bacoli95_mod, only: bacoli95_init, bacoli95, bacoli95_vals
    use bacoli95_mod, only: bacoli95_sol, bacoli95_sol_teardown
    use rootfinding, only: rt1, rt2, setsol

    implicit none
    integer, parameter     :: dp = kind(0d0)

    type(bacoli95_sol)     :: sol

    ! Choose npde to be consistent with the problem specific source 
    integer,  parameter    :: npde = 4
    integer,  parameter    :: num_points = 10       ! controls the x spacing
    real(dp), parameter    :: xa = -5.0, xb = 5.0
    real(dp)               :: tout, tstart, tstop, atol(npde), rtol(npde)
    integer, parameter     :: nrt = 1           ! number of root functions 
    integer, parameter     :: ntout = 10        ! controls the time spacing

    real(dp) :: jroot

    external f1, f2, bndxa, bndxb, uinit, derivf
    character(len=100) atol_str, rtol_str, fname

    integer nfev
    common /nfev/ nfev
    integer n_times_called
    common /n_times_called/ n_times_called

    integer :: direction ! 1 when we are going up and 0 when we are going down
    direction = 1 ! go up

    ! tolerances
    atol = 1.d-11
    rtol = atol

    ! Initialization: Allocate storage and set problem parameters. - call bacoli95_init(sol, npde, (/xa,xb/), atol=atol, rtol=rtol, dirichlet=1)
    call bacoli95_init(sol, npde, (/xa,xb/), atol=atol, rtol=rtol, estimator=1, kcol=5, nrt=nrt)
    ! Allocate memory for the work array to be used during the root finding
    call setsol(sol)

    tout = sol%t0     ! starts at 0
    tstop = 200d0        ! final time

    ! initialise number of function evaluations
    nfev = 0
    n_times_called = 0

    ! opening files
    write(atol_str, *) floor(log10(atol(1)))
    write(rtol_str, *) floor(log10(rtol(1)))
    fname = "BACOLIKR_EVENT_STATE_a_" // trim(adjustl(atol_str)) // "_r_" // trim(adjustl(rtol_str))
    open(unit=11, file=fname)

    ! Print solution output to file, modified to increase resolution of points
    do while (tout < tstop)

        if (direction == 1) then
            call bacoli95(sol, tstop, f1, bndxa, bndxb, uinit, rt=rt1)    ! integrate until the time step
        else
            call bacoli95(sol, tstop, f2, bndxa, bndxb, uinit, rt=rt2)    ! integrate until the time step
        end if
        if (sol%idid <= 0) goto 800

        tout = sol%t0
        
        if (sol%idid == 5) then
            ! we now need to set up bacoli95 for another call
            sol%idid = 1     ! make it restart integration
            sol%mflag(1) = 2 ! cold start

            ! root detected
            jroot = sol%jroot(1)
            
            !jroot = 1 or -1, means that rt1 or rt2 gave a root function
            if (jroot == 1 .or. jroot == -1) then
                ! switch direction to integrate with other root and model
                if (direction == 1) then
                    direction = 0
                else 
                    direction = 1
                end if
                ! write (*, *) "switched direction to", direction
            end if

            ! need to make a marker so that we can clean the plots ...
            write(11,*) "root" , tout
        end if
    end do

    ! close file
    close(11)
    
    write(*, *)

    print '("IDID       = ",i10)', sol%idid
    print '("nsteps     = ",i10)', sol%num_accepted_time_steps
    write (*, *) fname, "nfev = ", nfev, "n_times_called", n_times_called
    
   
    !-------------------------------------------------------------------

    call bacoli95_sol_teardown(sol) ; stop
600 print '("Error: Improperly formatted input")' ; stop
700 print '("Error: Could not allocate storage")' ; stop
800 print '("Error: Was not able to integrate to tsop")' ; stop

    !-------------------------------------------------------------------
    ! Formats!
900 format("kcol = ",i2,", nint0 = ",i4,", npde = ",i3,", tout = ",es7.1)
901 format("atol = ",es7.1,", rtol = ",es7.1,",",17x,a3)
902 format("Number of subintervals in the current mesh:",i8)
903 format("X(1) = ",f2.0)

end program

