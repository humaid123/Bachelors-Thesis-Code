
program time_disc_cross

    use bacoli95_mod, only: bacoli95_init, bacoli95, bacoli95_vals
    use bacoli95_mod, only: bacoli95_sol, bacoli95_sol_teardown

    implicit none
    integer, parameter     :: dp = kind(0d0)
    type(bacoli95_sol)     :: sol

    integer,  parameter    :: npde = 4
    real(dp), parameter    :: xa = -5.0, xb = 5.0
    real(dp)               :: tout, tstart, tstop, atol(npde), rtol(npde)
    real(dp), parameter    :: xplot = 0d0
    real(dp)               :: new_uout(npde)

    integer                :: j, k, ntout
    character(len=100)     :: fname 
    character(len=50)      :: npde_str, atol_str, rtol_str

    external F_time, bndxa, bndxb, uinit, derivf
    
    ! count the number of function evaluations...
    integer :: nfev
    common /nfev/ nfev
    nfev = 0

    !-------------------------------------------------------------------
    write(6,*) 'The number of PDEs is assumed to be ', npde
    write(6,*)
    
    !"the end of the temporal domain."
    tstop = 70d0

    !"At how many equally-spaced points along the time domain is output desired?"
    ntout = 1000

    !"Please choose an error tolerance"
    atol = 1.d-11
    rtol = atol

    ! Initialization: Allocate storage and set problem parameters.
    call bacoli95_init(sol, npde, (/xa,xb/), atol=atol, rtol=rtol, estimator=1, kcol=5)

    tstart = sol%t0


    k = 2    ! only print E
    write(npde_str,*) k
    write(atol_str,*) floor(log10(atol(1)))
    write(rtol_str,*) floor(log10(rtol(1)))

    fname = "BACOLIKR_TIME_"//trim(adjustl(npde_str))//"_a_"//trim(adjustl(atol_str))//"_r_"//trim(adjustl(rtol_str))
    write(*, *) "fname", fname
    open(unit=11,file=fname)

    !-------------------------------------------------------------------
    ! Integrate solution from t=0 to t=tout.
    print '(/"THE INPUT IS")'
    print 900, sol%kcol, sol%nint, sol%npde, tstop
    print 901, sol%atol(1), sol%rtol(1), "LOI"
    print 903, -sol%X(1)
    print 903, sol%x(sol%nint+1)

    ! Print solution output to file, modified to increase resolution of points
    do j = 2, ntout

        tout = tstart + (j-1)*(tstop-tstart)/(ntout-1)

        call bacoli95(sol, tout, F_time, bndxa, bndxb, uinit)
        if (sol%idid <= 0) goto 800

        if (j == 2) then
            call uinit(xplot, new_uout, npde)
            write(11, *) xplot, sol%t0, new_uout(2)
        end if

        call bacoli95_vals(sol, (/xplot/), new_uout)

        write(11, *) xplot, sol%t0, new_uout(2)
    end do

    print '("IDID       = ",i10)', sol%idid
    print '("nsteps     = ",i10)', sol%num_accepted_time_steps
    write (*, *) fname, "nfev = ", nfev 
    
    !-------------------------------------------------------------------

    call bacoli95_sol_teardown(sol) ; stop
800 print '("Error: Was not able to integrate to tsop")' ; stop

    !-------------------------------------------------------------------
    ! Formats!
900 format("kcol = ",i2,", nint0 = ",i4,", npde = ",i3,", tout = ",es7.1)
901 format("atol = ",es7.1,", rtol = ",es7.1,",",17x,a3)
903 format("X(1) = ",f2.0)

end program
