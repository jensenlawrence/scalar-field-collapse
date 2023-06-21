! Performs a simulation of spherically symmetric scalar field collapse
! in Einstein-aether theory with spatial compactification

! Original code by David Garfinkle
! Updated to Fortran 90 by Jensen Lawrence

! Program execution
program main
    implicit none

    ! Initialize parameters and variables
    integer n
    parameter(n=2000)
    real*8 k(n),ar(n),phi(n)
    real*8 dtk(n),dtar(n),dtphi(n)
    real*8 knew(n),arnew(n),phinew(n)
    real*8 dtknew(n),dtarnew(n),dtphinew(n)
    real*8 psi(n),p(n),dtpsi(n),dtp(n)
    real*8 psinew(n),pnew(n),dtpsinew(n),dtpnew(n)
    real*8 grr(n),dtgrr(n),grrnew(n),dtgrrnew(n)
    real*8 du,u(n),rs,r,dt
    integer itime,ntime,iter,i
    real*8 amp1,amp2,sigma,r0
    real*8 cvec(4)
    real*8 cnstr(n),cnstr2(n)
    real*8 alpha(n)
    real*8 w,lam,s1,s2,s3,wt
    real*8 trap(n),trap2(n)
    real*8 term1,term2
    integer ichk
    real*8 curv(n)
    real*8 mass(n),q(n)
    real*8 time
    real*8 cfl,cflmin,tmax
    real*8 tll(n),tnn(n)
    real*8 asq,trapmin,trapmin2
    integer iteramp
    real*8 c1,c2,c3,c4,c13,c14,c123,ck,v

    ! Open files for saving results
    open(unit=15, file='arfin.dat', status='unknown')
    open(unit=16, file='kfin.dat', status='unknown')
    open(unit=17, file='cnstrfin.dat', status='unknown')
    open(unit=18, file='Phifin.dat', status='unknown')
    open(unit=19, file='alphafin.dat', status='unknown')
    open(unit=21, file='trapfin.dat', status='unknown')
    open(unit=22, file='psifin.dat', status='unknown')
    open(unit=23, file='pfin.dat', status='unknown')
    open(unit=24, file='curvfin.dat', status='unknown')
    open(unit=25, file='massfin.dat', status='unknown')
    open(unit=26, file='grrfin.dat', status='unknown')
    open(unit=27, file='kfin2.dat', status='unknown')
    open(unit=28, file='arfin2.dat', status='unknown')
    open(unit=29, file='qfin2.dat', status='unknown')
    open(unit=30, file='trapfin2.dat', status='unknown')
    open(unit=31, file='cnstrfin2.dat', status='unknown')
    open(unit=32, file='asqfin.dat', status='unknown')
    open(unit=33, file='trapmin.dat', status='unknown')
    open(unit=34, file='trapmin2.dat', status='unknown')
    open(unit=35, file='tnnfin.dat', status='unknown')
    open(unit=36, file='attime.dat', status='unknown')
    open(unit=37, file='ttime.dat', status='unknown')
    open(unit=38, file='Psi.dat', status='unknown')
    open(unit=39, file='P.dat', status='unknown')
    open(unit=40, file='ar.dat', status='unknown')
    open(unit=41, file='K.dat', status='unknown')
    open(unit=42, file='alpha.dat', status='unknown')
    open(unit=43, file='Tll.dat', status='unknown')
    open(unit=44, file='Tnn.dat', status='unknown')
    open(unit=45, file='qfin.dat', status='unknown')

    ! Set simulation parameters
    time  = 0.0
    tmax  = 30.0
    ntime = 1000000
    rs    = 10.0
    du    = 1.0/float(n-1)
    dt    = 0.5*rs*du

    amp1  = 0.15
    amp2  = 0.0
    r0    = 10.0
    sigma = 4.0

    ! Calculate initial conditions
    call init(n, psi, p, ar, k, phi, grr, du, u, amp1, amp2, sigma, r0, cvec, rs)

    ! Set c values
    c1 = cvec(1)
    c2 = cvec(2)
    c3 = cvec(3)
    c4 = cvec(4)

    c13  = c1 + c3
    c14  = c1 + c4
    c123 = c1 + c2 + c3
    ck   = -1.0*c14*(1.0 - c13)/c123
    v    = sqrt((c14 - 2.0)/(ck*(2.0 + c13 + 3.0*c2)))

    ! Loop through time and evolve equations
    do itime = 1, ntime
        ichk = itime/1000
        ichk = ichk*1000
        if(ichk.eq.itime) then
            do i = 1, n/8
                r = rs*u(i)/(1.0 - u(i))
                write(38,77) r,time,psi(i)
                write(39,77) r,time,p(i)
                write(40,77) time,r,ar(i)
                write(41,77) time,r,k(i)
                write(42,77) r,time,alpha(i)
                write(43,77) r,time,tll(i)
                write(44,77) r,time,tnn(i)
                77 format(1x,3f15.7)
            end do

            write(38,99) 
            write(39,99) 
            write(40,99) 
            write(41,99) 
            write(42,99) 
            write(43,99) 
            write(44,99) 
            99 format(1x,'')

            write(*,*) itime,time
        end if

        ! Find the time derivatives for the current values of the variables
        call evolve(n, psi, p, ar, k, phi, grr, dtpsi, dtp, dtar, dtk, &
        & dtphi, dtgrr, du, dt, cvec, cnstr, alpha, trap, trap2, curv, &
        & mass, q, cnstr2, tll, tnn, rs)

        do i = 1, n
            psinew(i) = psi(i)
            pnew(i)   = p(i)
            arnew(i)  = ar(i)
            knew(i)   = k(i)
            phinew(i) = phi(i)
            grrnew(i) = grr(i)
        end do
    
        ! Crank-Nicholson iteration
        do iter = 1, 3

            ! Find derivatives at next time step
            call evolve(n, psinew, pnew, arnew, knew, phinew, grrnew, &
            & dtpsinew, dtpnew, dtarnew, dtknew, dtphinew, dtgrrnew, &
            & du, dt, cvec, cnstr, alpha, trap, trap2, curv, mass, q, &
            & cnstr2, tll, tnn, rs)

            ! Find the new values of the variables using the average
            ! of the time derivatives
            do i = 2, n-1
                psinew(i) = psi(i) + dt*0.5*(dtpsi(i) + dtpsinew(i))
                pnew(i)   = p(i) + dt*0.5*(dtp(i) + dtpnew(i))
                arnew(i)  = ar(i) + dt*0.5*(dtar(i) + dtarnew(i))
                knew(i)   = k(i) + dt*0.5*(dtk(i) + dtknew(i))
                phinew(i) = phi(i) + dt*0.5*(dtphi(i) + dtphinew(i))
                grrnew(i) = grr(i) + dt*0.5*(dtgrr(i) + dtgrrnew(i))
            end do

            ! Apply boundary conditions at u = 0
            psinew(1) = (4.0*psinew(2) - psinew(3))/3.0
            pnew(1)   = (4.0*pnew(2) - pnew(3))/3.0
            arnew(1)  = 0.0
            knew(1)   = (4.0*knew(2) - knew(3))/3.0
            grrnew(1) = (4.0*grrnew(2) - grrnew(3))/3.0
            phinew(1) = sqrt(grrnew(1))

            ! Apply boundary conditions at u = 1
            psinew(n) = 0.0
            pnew(n)   = 0.0
            arnew(n)  = 0.0
            knew(n)   = 0.0
            grrnew(n) = 1.0
            phinew(n) = 1.0
            ! radnew(n) = 2.0*radnew(n-1) - radnew(n-2)
            ! grrnew(n) = grrnew(n-1)
            ! grrnew(n) = grr(n) + dt*0.5*(dtgrr(n) + dtgrrnew(n))
            ! phinew(n) = phi(n) + dt*0.5*(dtphi(n) + dtphinew(n))

        end do

        do i = 1, n
            psi(i) = psinew(i)
            p(i)   = pnew(i)
            ar(i)  = arnew(i)
            k(i)   = knew(i)
            phi(i) = phinew(i)
            grr(i) = grrnew(i)
        end do

        time = time + dt

        trapmin  = trap(1)
        trapmin2 = trap2(1)

        do i = 2, n
            if(trap(i).lt.trapmin) then
                trapmin = trap(i)
            end if
            if(trap2(i).lt.trapmin2) then
                trapmin2 = trap2(i)
            end if
        end do

        write(33,30) time,trapmin
        write(34,30) time,trapmin2

        ! if(trapmin.lt.0.) then
        !     write(37,30) amp2,time
        !     write(*,45) amp2,time
        !     goto 99
        ! end if

        ! if(trapmin2.lt.0.) then
        !     write(36,30) amp2,time
        !     write(*,55) amp2,time
        !     goto 99
        ! end if

        ! 45 format(1x,'amp',3x,f12.6,4x,'found trapped surface at time',
        ! & 3x,f12.6)
        ! 55 format(1x,'amp',3x,f12.6,4x,'found anti-trapped surface at &
        ! & time',3x,f12.6)

        if(time.gt.tmax) goto 88

        ! implement the Courant condition
        cflmin = 1.0
        do i = 1, n
            cfl = sqrt(grr(i))/alpha(i)
            if(cfl.lt.cflmin) then
                cflmin = cfl
            end if
        end do

        cflmin = cflmin/v
        dt     = 0.5*rs*du*cflmin
    end do

    88 continue
    write(*,65)
    ! stop
    65 format(1x,'reached final time')

    ! we have reached the final time and will output the final values of the variables
    do i = 1, 3*n/4
        r=rs*u(i)/(1.0 - u(i))
        write(17,30) r,cnstr(i)
        write(31,30) r,cnstr2(i)
        write(15,30) r,ar(i)
        write(16,30) r,k(i)
        write(18,30) r,phi(i)
        write(19,30) r,alpha(i)
        write(21,30) r,trap(i)
        write(30,30) r,trap2(i)
        write(22,30) r,psi(i)
        write(23,30) r,p(i)
        asq=ar(i)*ar(i)/grr(i)
        write(32,30) r,asq
        write(35,30) r,tnn(i)
        write(45,30) r,q(i)
        ! if((rad(i).gt.1.5).and.(curv(i).gt.0.)) then
        ! write(24,30) log(rad(i)),log(curv(i))
        ! end if
        ! write(25,30) rad(i),mass(i)
        30 format(1x,2f16.8)
    end do

    ! do i=1138,n
    !     write(26,30) rad(i),grr(i)
    !     write(27,30) rad(i),k(i)
    !     write(28,30) rad(i),ar(i)
    !     write(29,30) rad(i),q(i)
    !     end do
    !     99 continue
    !     amp2=amp2+0.0125
    ! end do

end program main


! Subroutine to find the initial values of the variables
subroutine init(n, psi, p, ar, k, phi, grr, du, u, amp1, amp2, sigma, &
    & r0, cvec, rs)
    implicit none

    ! Initialize parameters and variables
    integer n
    real*8 psi(n), p(n), ar(n), k(n), phi(n), grr(n), du, u(n), rs
    real*8 rad(n), r
    integer i
    real*8 amp1, amp2, sigma, r0, rsq
    real*8 cvec(4), c1, c2, c3, c4
    real*8 c13, c14, c123, v, ck, exp1, exp2, term1, term2

    ! Put in the values of the c parameters
    c1 = 0.001
    c2 = 0.1
    c3 = 0.0
    c4 = 0.0

    cvec(1) = c1
    cvec(2) = c2
    cvec(3) = c3
    cvec(4) = c4

    c13  = c1 + c3
    c14  = c1 + c4
    c123 = c1 + c2 + c3
    ck   = -1.0*c14*(1.0 - c13)/c123
    v    = sqrt((c14 - 2.0)/(ck*(2.0 + c13 + 3.0*c2)))

    ! Set variables at t = 0 on the domain interior
    do i = 2, n-1
        u(i)   = du*float(i - 1)
        r      = rs*u(i)/(1.0 - u(i))
        psi(i) = amp1*exp(-1.0*((r**2 - r0**2)/sigma**2)**2)
        p(i)   = 0.0
        ar(i)  = 0.0
        k(i)   = 0.0
    end do

    ! Set variables at t = 0 at u = 0
    u(1)   = 0.0
    psi(1) = (4.0*psi(2) - psi(3))/3.0
    p(1)   = (4.0*p(2) - p(3))/3.0
    ar(1)  = 0.0
    k(1)   = (4.0*k(2) - k(3))/3.0

    ! Set variables at t = 0 at u = 1
    psi(n) = 0.0
    ar(n)  = 0.0
    k(n)   = 0.0

    ! Solve the Hamiltonian constraint to find the metric component
    call calcphi(n, phi, du, psi, rs)
    do i = 1, n
        grr(i) = phi(i)*phi(i)
    end do

    return
end subroutine init


! Subroutine to calculate the time derivatives of the variables
subroutine evolve(n, psi, p, ar, k, phi, grr, dtpsi, dtp, dtar, dtk, &
    & dtphi, dtgrr, du, dt, cvec, cnstr, alpha, trap, trap2, curv, mass, &
    & q, cnstr2, tll, tnn, rs)
    implicit none

    ! Initialize parameters and variables
    integer n
    real*8 psi(n), p(n), ar(n), k(n), phi(n), grr(n)
    real*8 dtpsi(n), dtp(n), dtar(n), dtk(n), dtphi(n), dtgrr(n)
    real*8 du, dt, rs
    integer i
    real*8 duk, duar, duphi, duuphi, dupsi, dup, dugrr
    real*8 drk, drar, drphi, drrphi, r, rad, drrad, drrrad
    real*8 drpsi, drp, drgrr, lappsi, diva
    real*8 up, um, dupsip, dupsim
    real*8 eps, arko(n), kko(n), psiko(n), pko(n), phiko(n), grrko(n)
    real*8 alpha(n), trap(n), trap2(n), q(n), cvec(4)
    real*8 c1, c2, c3, c4, c13, c14
    real*8 term1, term2, term3, term4, term5, term6, term7, term8, term9
    real*8 cnstr(n)
    real*8 cof1, cof2, cof3
    real*8 curv(n), mass(n)
    real*8 cnstr2(n), w(n), intw
    real*8 rho, js, s, mss, tll(n), tnn(n)
    real*8 u(n)

    ! Unpack the values of the c parameters
    c1 = cvec(1)
    c2 = cvec(2)
    c3 = cvec(3)
    c4 = cvec(4)

    c13 = c1 + c3
    c14 = c1 + c4

    ! Set u values
    do i = 1, n
        u(i) = du*float(i - 1)
    end do

    ! Compute the Kreiss-Oliger dissipation terms for extra stability
    eps = 0.5
    do i = 1, n
        arko(i)  = 0.0
        kko(i)   = 0.0
        psiko(i) = 0.0
        pko(i)   = 0.0
        phiko(i) = 0.0
        grrko(i) = 0.0
    end do

    do i = 3, n-2
        arko(i) = (-1.0/16.0)*(eps/dt)*(ar(i+2) + ar(i-2) + 6.0*ar(i) &
        & - 4.0*(ar(i+1) + ar(i-1)))
        kko(i) = (-1.0/16.0)*(eps/dt)*(k(i+2) + k(i-2) + 6.0*k(i) &
        & - 4.0*(k(i+1) + k(i-1)))
        psiko(i) = (-1.0/16.0)*(eps/dt)*(psi(i+2) + psi(i-2) + 6.0*psi(i) &
        & - 4.0*(psi(i+1) + psi(i-1)))
        pko(i) = (-1.0/16.0)*(eps/dt)*(p(i+2) + p(i-2) + 6.0*p(i)&
        & - 4.0*(p(i+1) + p(i-1))) 
        phiko(i) = (-1.0/16.0)*(eps/dt)*(phi(i+2) + phi(i-2) + 6.0*phi(i)&
        & - 4.0*(phi(i+1) + phi(i-1))) 
        grrko(i) = (-1.0/16.0)*(eps/dt)*(grr(i+2) + grr(i-2) + 6.0*grr(i)&
        & - 4.0*(grr(i+1) + grr(i-1)))
    end do

    ! Find the lapse, the trace-free part of the extrinsic curvature,
    ! and the shift
    call calcalpha(n, alpha, ar, du, u, rs)
    call calcq(n, q, phi, k, du, cvec, psi, p, u)

    ! Time evolution calculations
    do i = 2, n-1
        ! Calculate the spatial derivatives of the variables
        duk    = (k(i+1) - k(i-1))/(2.0*du)
        duar   = (ar(i+1) - ar(i-1))/(2.0*du)
        duphi  = (phi(i+1) - phi(i-1))/(2.0*du)
        drk    = duk*(1.0 - u(i))*(1.0 - u(i))/rs
        drar   = duar*(1.0 - u(i))*(1.0 - u(i))/rs
        drphi  = duphi*(1.0 - u(i))*(1.0 - u(i))/rs
        duuphi = (phi(i+1) + phi(i-1) - 2.0*phi(i))/du**2

        up    = u(i+1)
        um    = u(i-1)
        term1 = (1.5/(du*(up*up + up*um + um*um))) &
        & * (up*up*ar(i+1) - um*um*ar(i-1))*((1.0 - u(i))**2)/rs
        term2 = 2.0*(1.0 - u(i))*ar(i)/rs
        term3 = 2.0*ar(i)*drphi/phi(i) 
        diva  = term1 + term2 + term3

        r      = rs*u(i)/(1.0 - u(i))
        rad    = r*phi(i)
        drrad  = phi(i) + u(i)*(1.0 - u(i))*duphi
        dupsi  = (psi(i+1) - psi(i-1))/(2.0*du)
        dup    = (p(i+1) - p(i-1))/(2.0*du)
        drpsi  = dupsi*(1.0 - u(i))*(1.0 - u(i))/rs
        drp    = dup*(1.0 - u(i))*(1.0 - u(i))/rs       
        dugrr  = (grr(i+1) - grr(i-1))/(2.0*du)
        drgrr  = dugrr*(1.0 - u(i))*(1.0 - u(i))/rs
        drrrad = ((1.0 - u(i))**3)*(u(i)*duuphi + 2.0*duphi)/rs

        up     = u(i) + 0.5*du
        um     = u(i) - 0.5*du
        dupsip = (psi(i+1) - psi(i))/du
        dupsim = (psi(i) - psi(i-1))/du
        term1  = (3.0/(du*(up*up + up*um + um*um))) &
        & * (up*up*dupsip - um*um*dupsim)*((1.0 - u(i))**4)/(rs*rs)
        term2  = 2.0*drpsi*drphi/phi(i) 
        lappsi = term1 + term2

        ! Calculate the time derivatives of the variables
        term1    = alpha(i)*p(i)
        term2    = psiko(i)
        dtpsi(i) = term1 + term2

        term1  = alpha(i)*(p(i)*k(i) + ar(i)*drpsi/grr(i))
        term2  = alpha(i)*lappsi/grr(i)
        term3  = -0.5*alpha(i)*drgrr*drpsi/(grr(i)*grr(i))
        term4  = pko(i)
        dtp(i) = term1 + term2 + term3 + term4

        cof1    = c13/(c14*(1. - c13))
        cof2    = -1.0*(c2 + c13)/(c14*(1.0 - c13))
        term1   = alpha(i)*(k(i)/3.0 - 2.0*q(i))*ar(i)
        term2   = alpha(i)*cof1*p(i)*drpsi
        term3   = alpha(i)*cof2*drk
        term4   = arko(i)
        dtar(i) = term1 + term2 + term3 + term4

        cof1   = (c14 - 2.0)/(2.0 + c13 + 3.0*c2)
        cof2   = 2.0/(2.0 + c13 + 3.0*c2)
        cof3   = 3.0*(1.0 - c13)/(2.0 + c13 + 3.0*c2)
        term1  = (alpha(i)/3.0)*k(i)*k(i)
        term2  = alpha(i)*cof1*(diva + ar(i)*ar(i))/grr(i)
        term3  = -0.5*alpha(i)*cof1*drgrr*ar(i)/(grr(i)*grr(i))   
        term4  = alpha(i)*cof2*p(i)*p(i)
        term5  = alpha(i)*cof3*q(i)*q(i)
        term6  = kko(i)
        dtk(i) = term1 + term2 + term3 + term4 + term5 + term6

        term1    = alpha(i)*phi(i)*(q(i)/2.0 - k(i)/3.0)
        term2    = phiko(i)
        dtphi(i) = term1 + term2

        term1    = -2.0*alpha(i)*grr(i)*(q(i) + k(i)/3.0)
        term2    = grrko(i)
        dtgrr(i) = term1 + term2

        ! Calculate the constraint to use as a code check
        term1    = drrrad
        term2    = (0.5/rad)*(drrad**2 - grr(i))
        term3    = (-0.5/grr(i))*drgrr*drrad
        term4    = c14*rad*(0.5*drar + 0.25*ar(i)**2)
        term5    = c14*ar(i)*drrad
        term6    = -0.25*(rad*c14/grr(i))*drgrr*ar(i)
        term7    = (3.0/8.0)*rad*grr(i)*(1.0 - c13)*(q(i)**2)  
        term8    = (-1.0/12.0)*rad*grr(i)*(2.0 + c13 + 3.0*c2)*(k(i)**2)  
        term9    = 0.25*rad*(grr(i)*(p(i)**2) + drpsi**2)
        cnstr(i) = term1 + term2 + term3 + term4 + term5 + term6 + term7 &
        & + term8 + term9
        w(i)     = (2.0*rad*drrad/grr(i))*(term4 + term5 + term6 + term7 &
        & + term8 + term9)
        cnstr(i) = cnstr(i)*rad/(10.0+rad)
 
        ! Calculate the mass and check for the presence of two kinds
        ! of trapped surfaces
        trap(i)  = drrad/sqrt(grr(i)) + rad*(q(i)/2.0 - k(i)/3.0) 
        trap2(i) = drrad/sqrt(grr(i)) - rad*(q(i)/2.0 - k(i)/3.0) 
        ! trap2(i) = rad(i)*(q(i)/2.0 - k(i)/3.0) + drrad*sqrt(2.0)
        ! grr(i)   = drrad**2 - (rad(i)**2)*((q(i)/2.0 - k(i)/3.0)**2)
        mass(i) = (rad/2.0)*(1.0 - rad)

        ! Calculate the components of the stress-energy
        term1 = (0.5*c14/grr(i))*(2.0*drar + ar(i)*(4.0*drrad/rad + ar(i)))
        term2 = 0.5*(p(i)*p(i) + drpsi*drpsi/grr(i))
        term3 = (-1.0/6.0)*(c13 + 3.0*c2)*k(i)*k(i)
        term4 = -0.75*c13*q(i)*q(i)
        term5 = -0.5*c14*ar(i)*drgrr/(grr(i)*grr(i))
        rho   = term1 + term2 + term3 + term4 + term5
        js    = (1.0/(1.0 - c13))*((c13 + c2)*drk - p(i)*drpsi)/(sqrt(grr(i)))

        term1 = 2.0*p(i)*p(i) - 4.5*(c13 + c2)*q(i)*q(i)
        term2 = (c14 + c13 + 2.*c2)*(drar + ar(i)*(2.0*drrad/rad + ar(i))) &
        & / grr(i)
        term3 = -0.5*(c14 + c13 + 2.0*c2)*ar(i)*drgrr/(grr(i)*grr(i))
        s     = -1.0*rho + (2.0/(2.0 + c13 + 3.0*c2))*(term1 + term2 + term3)

        term1 = drpsi*drpsi + (c3 - c4)*ar(i)*ar(i)
        term2 = c13*(drar + ar(i)*drrad/rad)
        term3 = (c13/rad)*(drrrad + (grr(i) - drrad*drrad)/rad)
        term4 = (-0.5*drgrr/grr(i))*(c13*ar(i) + drrad/rad)
        mss   = (2.0*grr(i)/(3.0*(1.0 - c13)))*(term1 + term2 + term3 + term4)

        tll(i) = rho + (s/3.0) + mss - 2.0*js
        tnn(i) = rho + (s/3.0) + mss + 2.0*js
    end do

    ! First constraint at boundaries
    cnstr(1)   = cnstr(2)
    cnstr(n-1) = cnstr(n-2)
    cnstr(n)   = cnstr(n-1)

    ! Trapped surfaces at boundaries
    trap(1)  = 1.0
    trap(n)  = trap(n-1)
    trap2(1) = sqrt(2.0)
    trap2(n) = trap2(n-1)

    ! Mass at boundaries
    mass(1) = 0.0
    mass(n) = mass(n-1)

    ! dtphi(n) = alpha(n)*phi(n)*(q(n)/2.0 - k(n)/3.0)
    ! dtgrr(n) = -2.0*alpha(n)*grr(n)*(q(n) + k(n)/3.0)
    ! grr(1)   = 1.0
    ! grr(n)   = grr(n-1)

    ! Calculate the second constraint
    w(1) = 0.0
    intw = 0.0
    do i = 2, n-2
        r         = rs*u(i)/(1.0 - u(i))
        rad       = r*phi(i)
        duphi     = (phi(i+1) - phi(i-1))/(2.0*du)
        drrad     = phi(i) + u(i)*(1.0-u(i))*duphi
        intw      = intw + 0.5*(du*rs/((1.0 - u(i))**2))*(w(i) + w(i-1))
        cnstr2(i) = rad*(1.0 - drrad*drrad/grr(i)) - intw
        ! cnstr2(i) = 4.0*cnstr2(i)
    end do

    ! Second constraint at boundaries
    cnstr2(1)   = 0.0
    cnstr2(n-1) = 0.0
    cnstr2(n)   = 0.0

    ! Find the squared curvature (this blows up at a singularity)
    ! call calccurv(n, curv, ar, k, q, rad, psi, p, dr, cvec)

    return
end subroutine evolve


! Subroutine to calculate the squared curvature
subroutine calccurv(n, curv, ar, k, q, rad, psi, p, dr, cvec)
    implicit none

    ! Initialize parameters and variables
    integer n
    real*8 curv(n), ar(n), k(n), q(n), rad(n), psi(n), p(n), dr, cvec(4)
    integer i
    real*8 c1, c2, c3, c4, c13, c14
    real*8 y, r1, r2, k1, k2, l1, l2, m1, m2
    real*8 term1, term2, term3, term4, term5
    real*8 drk, drar, drrad, drpsi, drrrad
    real*8 luk

    ! Unpack the values of the c parameters
    c1 = cvec(1)
    c2 = cvec(2)
    c3 = cvec(3)
    c4 = cvec(4)

    c13 = c1 + c3
    c14 = c1 + c4

    ! Calculate squared curvature
    do i = 2, n-1
        drk    = (k(i+1) - k(i-1))/(2.0*dr)
        drar   = (ar(i+1) - ar(i-1))/(2.0*dr)
        drrad  = (rad(i+1) - rad(i-1))/(2.0*dr)
        drpsi  = (psi(i+1) - psi(i-1))/(2.0*dr)
        drrrad = (rad(i+1) + rad(i-1) - 2.0*rad(i))/(dr*dr)

        y  = (p(i)*drpsi - (c13 + c2)*drk)/(1.0 - c13)
        r1 = -2.0*drrrad/rad(i)
        r2 = (1.0 - (rad(i)*drrrad + drrad*drrad))/(rad(i)*rad(i))
        k1 = q(i) + (k(i)/3.0)
        k2 = (k(i)/3.0) - (q(i)/2.0)
        l1 = r1 + 2.0*k1*k2
        l2 = r2 + k2*(k1 + k2)

        term1 = (k(i)*k(i))/3.0
        term2 = (c14 - 2.0)*(drar + ar(i)*(ar(i) + 2.0*drrad/rad(i)))
        term3 = 2.0*p(i)*p(i)
        term4 = 3.0*(1.0 - c13)*q(i)*q(i)
        luk   = term1 + (term2 + term3 + term4)/(2.0 + c13 + 3.0*c2)

        term1 = -1.0*l2
        term2 = 0.5*(c1 + c4)*(drar + ar(i)*(ar(i) + 2.0*drrad/rad(i)))
        term3 = 0.5*(c13 + c2)*(luk - k(i)*k(i))
        term4 = c13*(k(i)*k2 + ar(i)*drrad/rad(i))
        m2    = (term1 + term2 + term3 + term4)/(c13 - 1.0)

        term1 = (drar + ar(i)*(ar(i) + 2.0*drrad/rad(i)))
        term2 = luk - 2.0*m2
        m1    = term1 + term2

        term1   = 3.0*l1*l1 + 4.0*l2*(l2 - l1)
        term2   = 4.0*m1*m1 - 8.0*m2*m2
        term3   = -4.0*y*y
        curv(i) = term1 + term2 + term3
    end do

    ! Squared curvature at boundaries
    curv(1) = curv(2)
    curv(n) = curv(n-1)

    return
end subroutine calccurv


! Subroutine to calculate the lapse
subroutine calcalpha(n, alpha, ar, du, u, rs)
    implicit none

    ! Initialize parameters and variables
    integer n
    real*8 alpha(n), ar(n), du, u(n), rs
    real*8 lnalpha(n), uavg, aravg
    integer i

    ! Calculate alpha by integrating the acceleration of the aether field
    lnalpha(n) = 0.0
    do i = n, 2, -1
        aravg        = 0.5*(ar(i) + ar(i-1))
        uavg         = 0.5*(u(i) + u(i-1))
        lnalpha(i-1) = lnalpha(i) - du*rs*aravg/(1.0 - uavg)**2
    end do

    ! Convert from ln(alpha) to alpha
    do i = 1, n
        alpha(i) = exp(lnalpha(i))
    end do

    return
end subroutine calcalpha


! Subroutine to calculate the trace-free part of the extrinsic curvature
subroutine calcq(n, q, phi, k, du, cvec, psi, p, u)
    implicit none

    ! Initialize parameters and variables
    integer n
    real*8 q(n),phi(n),k(n),du,cvec(4),psi(n),p(n),u(n)
    real*8 c2,c13
    real*8 uavg,phiavg,duphi,duk,qint
    integer i
    real*8 cof,kavg
    real*8 dupsi,pavg
    real*8 term1,term2

    ! Unpack c values
    c2  = cvec(2)
    c13 = cvec(1) + cvec(3)
    cof = (2.0 + c13 + 3.0*c2)/(3.0 - 3.0*c13)

    ! Calculate trace-free part of the extrinsic curvature by 
    ! integrating the momentum constraint
    q(1) = 0.0
    do i = 1, n-1
        uavg   = du*(float(i) - 0.5)
        phiavg = 0.5*(phi(i) + phi(i+1))
        kavg   = 0.5*(k(i) + k(i+1))
        pavg   = 0.5*(p(i) + p(i+1))
        duphi  = (phi(i+1) - phi(i))/du
        duk    = (k(i+1) - k(i))/du
        dupsi  = (psi(i+1) - psi(i))/du

        term1 = -3.0*q(i)*(1.0/(uavg*(1.0 - uavg)) + duphi/phiavg)
        term2 = cof*duk - pavg*dupsi/(1.0 - c13) 
        qint  = q(i) + 0.5*du*(term1 + term2)

        term1  = -3.0*qint*(1.0/(uavg*(1.0 - uavg)) + duphi/phiavg)
        q(i+1) = q(i) + du*(term1 + term2)
    end do

    return
end subroutine calcq


! Subroutine to calculate phi
subroutine calcphi(n, phi, du, psi, rs)
    implicit none

    ! Initialize parameters and variables
    integer n, i
    real*8 phi(n), du, psi(n), rs
    real*8 w, wint, phiint, u, dupsi, term1, term2, term3

    ! Calculate phi
    phi(1) = 1.0
    w = 1.0
    do i = 1, n-1
        u      = du*(float(i) - 0.5)
        dupsi  = (psi(i+1) - psi(i))/du
        phiint = phi(i) + 0.5*du*(w - phi(i))/(u*(1.0 - u))

        term1    = (w - phi(i))*(w - phi(i))/(2.0*phi(i)*u*(1.0 - u))
        term2    = -0.25*phi(i)*u*(1.0 - u)*dupsi*dupsi
        wint     = w + 0.5*du*(term1 + term2)
        phi(i+1) = phi(i) + du*(wint - phiint)/(u*(1.0 - u))

        term1 = (wint - phiint)*(wint - phiint)/(2.0*phiint*u*(1.0 - u))
        term2 = -0.25*phiint*u*(1.0 - u)*dupsi*dupsi
        w     = w + du*(term1 + term2)
    end do

    return
end subroutine calcphi