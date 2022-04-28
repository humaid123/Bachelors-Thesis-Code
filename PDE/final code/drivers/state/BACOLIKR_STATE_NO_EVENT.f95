program trimesh_example_driver

    use bacoli95_mod, only: bacoli95_init, bacoli95, bacoli95_vals
    use bacoli95_mod, only: bacoli95_sol, bacoli95_sol_teardown

    implicit none
    integer, parameter     :: dp = kind(0d0)

    type(bacoli95_sol)     :: sol

    ! Choose npde to be consistent with the problem specific source 
    integer,  parameter    :: npde = 4
    integer,  parameter    :: num_points = 100       ! controls the x spacing
    real(dp), parameter    :: xa = -5.0, xb = 5.0
    real(dp)               :: u(npde, num_points + 1, 1)
    real(dp)               :: tout, tstart, tstop, atol(npde), rtol(npde)
    real(dp)               :: xsol(num_points + 1)
    integer                :: i, j, k, ier
    integer, parameter     :: ntout = 10000        ! controls the time spacing
    integer m                                   ! the integral loop index
    real(dp) :: approx, h   
    external F_state, bndxa, bndxb, uinit, derivf
    character(len=100)      :: str, fname, atol_str, rtol_str ! required to create descriptive file names

    ! variables to plot the cross-section
    real(dp), parameter :: xplot(1)  = (/ 0 /)
    real(dp) yplot(npde, 1, 1) 

    ! count the number of function evaluations
    integer nfev
    common /nfev/ nfev

    !NEW => changed from = real(dp), parameter :: min_val = 1, max_val = 5
    real(dp), parameter :: min_val = 10000, max_val = 30000
    integer measures_implemented        ! 1 for true, 0 for false
    ! pipe around beta to the model
    common /measures_implemented/ measures_implemented
    measures_implemented = 0

    ! tolerances
    atol = 1.d-11
    rtol = atol

    ! init nfev count
    nfev = 0

    ! Initialization: Allocate storage and set problem parameters. - call bacoli95_init(sol, npde, (/xa,xb/), atol=atol, rtol=rtol, dirichlet=1)
    call bacoli95_init(sol, npde, (/xa,xb/), atol=atol, rtol=rtol, estimator=1, kcol=5)

    tstart = sol%t0     ! starts at 0
    tstop = 200d0        ! final time

    ! calculate the uniform x-values where the integral are going to be approximated
    xsol(1) = xa
    h = (xb - xa) / num_points
    do i = 2, num_points + 1
        xsol(i) = xsol(i - 1) + h
    end do

    ! opening file
    write(str,*) ntout
    write(atol_str,*) floor(log10(atol(1)))
    write(rtol_str,*) floor(log10(rtol(1)))
    fname = 'BACOLIKR_STATE_' // trim(adjustl(str)) // "_a_" // trim(adjustl(atol_str)) // "_rtol_" // trim(adjustl(rtol_str))
    open(unit=11, file=fname)

    k = 2 ! we are integrating E
    ! Print solution output to file, modified to increase resolution of points
    do j = 2, ntout
        ! at what time to stop based on the time spacing
        tout = tstart + (j-1)*(tstop-tstart)/(ntout-1)

        call bacoli95(sol, tout, F_state, bndxa, bndxb, uinit)    ! integrate until the time step
        if (sol%idid <= 0) goto 800

        ! this call is to calculate for the integration
        call bacoli95_vals(sol, xsol, u)    ! calculate the values based on the xsol uniform mesh and place in u
        
        ! integrate
        approx = 0d0
        do m = 1, num_points
            approx = approx + (xsol(m+1) - xsol(m)) * (u(k, m+1, 1) + u(k, m, 1))
        end do
        approx = (approx/2.0d0)

        ! this call to _vals() is for the cross section
        call bacoli95_vals(sol, xplot, yplot)
        ! write(*, *) "[", tout, ",", yplot(k, 1, 1), "],"
        write(11, *) tout, yplot(k, 1, 1), approx

        if (measures_implemented == 1) then
            if (approx < min_val) then
                measures_implemented = 0
            end if
        else
            if (approx > max_val) then
                measures_implemented = 1
            end if
        end if
    end do

    close(11)
    write(*, *)

    print '("IDID       = ",i10)', sol%idid
    print '("nsteps     = ",i10)', sol%num_accepted_time_steps
    write (*, *) fname, "nfev = ", nfev
    
   
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
