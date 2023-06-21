c  program to evolve the Einstein-Aether equations
c  with additional scalar field matter
c
       implicit none
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
       open(unit=15,file='arfin.dat',status='unknown')
       open(unit=16,file='kfin.dat',status='unknown')
       open(unit=17,file='cnstrfin.dat',status='unknown')
       open(unit=18,file='Phifin.dat',status='unknown')
       open(unit=19,file='alphafin.dat',status='unknown')
       open(unit=21,file='trapfin.dat',status='unknown')
       open(unit=22,file='psifin.dat',status='unknown')
       open(unit=23,file='pfin.dat',status='unknown')
       open(unit=24,file='curvfin.dat',status='unknown')
       open(unit=25,file='massfin.dat',status='unknown')
       open(unit=26,file='grrfin.dat',status='unknown')
       open(unit=27,file='kfin2.dat',status='unknown')
       open(unit=28,file='arfin2.dat',status='unknown')
       open(unit=29,file='qfin2.dat',status='unknown')
       open(unit=30,file='trapfin2.dat',status='unknown')
       open(unit=31,file='cnstrfin2.dat',status='unknown')
       open(unit=32,file='asqfin.dat',status='unknown')
       open(unit=33,file='trapmin.dat',status='unknown')
       open(unit=34,file='trapmin2.dat',status='unknown')
       open(unit=35,file='tnnfin.dat',status='unknown')
       open(unit=36,file='attime.dat',status='unknown')
       open(unit=37,file='ttime.dat',status='unknown')
       open(unit=38,file='Psi.dat',status='unknown')
       open(unit=39,file='P.dat',status='unknown')
       open(unit=40,file='ar.dat',status='unknown')
       open(unit=41,file='K.dat',status='unknown')
       open(unit=42,file='alpha.dat',status='unknown')
       open(unit=43,file='Tll.dat',status='unknown')
       open(unit=44,file='Tnn.dat',status='unknown')
       open(unit=45,file='qfin.dat',status='unknown')
c
       ntime=1000000
       rs=10.
       amp1=0.15
       amp2=0.
       r0=10.0
       sigma=4.
       du=1./float(n-1)
       dt=0.5*rs*du
c       do iteramp=1,81
       time=0.
c       tmax=0.1
       tmax=20.
c
c   find the initial data
c
       call init(n,psi,p,ar,k,phi,grr,du,u,amp1,amp2,sigma,r0,cvec,rs)
       c1=cvec(1)
       c2=cvec(2)
       c3=cvec(3)
       c4=cvec(4)
       c13=c1+c3
       c14=c1+c4
       c123=c1+c2+c3
       ck=-1.*c14*(1.-c13)/c123
       v=sqrt((c14-2.)/(ck*(2.+c13+3.*c2)))
c
c   loop through time
c
       do itime=1,ntime
       ichk=itime/1000
       ichk=ichk*1000
       if(ichk.eq.itime) then
       do i=1,n/8
       r=rs*u(i)/(1.-u(i))
       write(38,77) r,time,psi(i)
       write(39,77) r,time,p(i)
       write(40,77) time,r,ar(i)
       write(41,77) time,r,k(i)
       write(42,77) r,time,alpha(i)
       write(43,77) r,time,tll(i)
       write(44,77) r,time,tnn(i)
  77   format(1x,3f15.7)
       enddo
       write(38,99) 
       write(39,99) 
       write(40,99) 
       write(41,99) 
       write(42,99) 
       write(43,99) 
       write(44,99) 
  99   format(1x,'')
       write(*,*) itime,time
       call xvs('Psi',time,u,psi,n)
       call xvs('P',time,u,p,n)
       call xvs('ar',time,u,ar,n)
       call xvs('K',time,u,k,n)
       call xvs('alpha',time,u,alpha,n)
       call xvs('grr',time,u,grr,n)
       call xvs('trap',time,u,trap,n)
       call xvs('trap2',time,u,trap2,n)
       call xvs('rad',time,u,Phi,n)
       call xvs('Tll',time,u,tll,n)
       call xvs('Tnn',time,u,tnn,n)
       endif
c
c   find the time derivatives for the current value of the variables
c
       call evolve(n,psi,p,ar,k,phi,grr,dtpsi,dtp,dtar,dtk,dtphi,
     & dtgrr,du,dt,cvec,cnstr,alpha,trap,trap2,curv,mass,q,cnstr2,
     & tll,tnn,rs)
       do i=1,n
       psinew(i)=psi(i)
       pnew(i)=p(i)
       arnew(i)=ar(i)
       knew(i)=k(i)
       phinew(i)=phi(i)
       grrnew(i)=grr(i)
       enddo
c
c   do the Crank-Nicholson iteration
c
       do iter=1,3
c
c   find the time derivatives at the next time step
c
       call evolve(n,psinew,pnew,arnew,knew,phinew,grrnew,dtpsinew,
     & dtpnew,dtarnew,dtknew,dtphinew,dtgrrnew,du,dt,cvec,cnstr,alpha,
     & trap,trap2,curv,mass,q,cnstr2,tll,tnn,rs)
c
c   find the new values of the variables using the 
c   average of the time derivatives
c
       do i=2,n-1
       psinew(i)=psi(i)+dt*0.5*(dtpsi(i)+dtpsinew(i))
       pnew(i)=p(i)+dt*0.5*(dtp(i)+dtpnew(i))
       arnew(i)=ar(i)+dt*0.5*(dtar(i)+dtarnew(i))
       knew(i)=k(i)+dt*0.5*(dtk(i)+dtknew(i))
       phinew(i)=phi(i)+dt*0.5*(dtphi(i)+dtphinew(i))
       grrnew(i)=grr(i)+dt*0.5*(dtgrr(i)+dtgrrnew(i))
       enddo
c
c   apply boundary conditions at r=0 and r=rmax
c
       psinew(1)=(4.*psinew(2)-psinew(3))/3.
       pnew(1)=(4.*pnew(2)-pnew(3))/3.
       arnew(1)=0.
       knew(1)=(4.*knew(2)-knew(3))/3.
       grrnew(1)=(4.*grrnew(2)-grrnew(3))/3.
       phinew(1)=sqrt(grrnew(1))
c
       psinew(n)=0.
       pnew(n)=0.
       arnew(n)=0.
       knew(n)=0.
       grrnew(n)=1.
       phinew(n)=1.
c       radnew(n)=2.*radnew(n-1)-radnew(n-2)
c       grrnew(n)=grrnew(n-1)
c       grrnew(n)=grr(n)+dt*0.5*(dtgrr(n)+dtgrrnew(n))
c       phinew(n)=phi(n)+dt*0.5*(dtphi(n)+dtphinew(n))
       enddo
       do i=1,n
       psi(i)=psinew(i)
       p(i)=pnew(i)
       ar(i)=arnew(i)
       k(i)=knew(i)
       phi(i)=phinew(i)
       grr(i)=grrnew(i)
       enddo
       time=time+dt
       trapmin=trap(1)
       trapmin2=trap2(1)
       do i=2,n
       if(trap(i).lt.trapmin) then
       trapmin=trap(i)
       endif
       if(trap2(i).lt.trapmin2) then
       trapmin2=trap2(i)
       endif
       enddo
       write(33,30) time,trapmin
       write(34,30) time,trapmin2
c       if(trapmin.lt.0.) then
c       write(37,30) amp2,time
c       write(*,45) amp2,time
c       goto 99
c       endif
c       if(trapmin2.lt.0.) then
c       write(36,30) amp2,time
c       write(*,55) amp2,time
c       goto 99
c       endif
c 45    format(1x,'amp',3x,f12.6,4x,'found trapped surface at time',
c     & 3x,f12.6)
c 55    format(1x,'amp',3x,f12.6,4x,'found anti-trapped surface at 
c     & time',3x,f12.6)
       if(time.gt.tmax) goto 88
c
c      implement the Courant condition
c
       cflmin=1.
       do i=1,n
       cfl=sqrt(grr(i))/alpha(i)
       if(cfl.lt.cflmin) then
       cflmin=cfl
       endif
       enddo
       cflmin=cflmin/v
       dt=0.5*rs*du*cflmin
       enddo
  88   continue
       write(*,65)
c       stop
  65   format(1x,'reached final time')
c
c   we have reached the final time and will output 
c   the final values of the variables
c
       do i=1,3*n/4
       r=rs*u(i)/(1.-u(i))
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
c       if((rad(i).gt.1.5).and.(curv(i).gt.0.)) then
c       write(24,30) log(rad(i)),log(curv(i))
c       endif
c       write(25,30) rad(i),mass(i)
  30   format(1x,2f16.8)
       enddo
c       do i=1138,n
c       write(26,30) rad(i),grr(i)
c       write(27,30) rad(i),k(i)
c       write(28,30) rad(i),ar(i)
c       write(29,30) rad(i),q(i)
c       enddo
c  99   continue
c       amp2=amp2+0.0125
c       enddo
       end
c
       subroutine init(n,psi,p,ar,k,phi,grr,du,u,amp1,amp2,sigma,
     & r0,cvec,rs)
       implicit none
c
c   subroutine to find the initial data
c
       integer n
       real*8 psi(n),p(n),ar(n),k(n),phi(n),grr(n),du,u(n),rs
       real*8 rad(n),r
       integer i
       real*8 amp1,amp2,sigma,r0,rsq
       real*8 cvec(4),c1,c2,c3,c4
       real*8 c13,c14,c123,v,ck,exp1,exp2,term1,term2
c
c   put in the values of the c parameters
c
c       c1=0.7
       c1=0.001
       c2=0.1
       c3=0.
       c4=0.
       cvec(1)=c1
       cvec(2)=c2
       cvec(3)=c3
       cvec(4)=c4
       c13=c1+c3
       c14=c1+c4
       c123=c1+c2+c3
       ck=-1.*c14*(1.-c13)/c123
       v=sqrt((c14-2.)/(ck*(2.+c13+3.*c2)))
       u(1)=0.
       do i=2,n-1
       u(i)=du*float(i-1)
       r=rs*u(i)/(1.-u(i))
c
c   the scalar field and the aether field are ingoing a gaussians
c   
c      
       psi(i)=amp1*exp(-1.*((r**2-r0**2)/sigma**2)**2)
       p(i)=0.
       ar(i)=0.
       k(i)=0.
       enddo
       psi(1)=(4.*psi(2)-psi(3))/3.
       p(1)=(4.*p(2)-p(3))/3.
       ar(1)=0.
       k(1)=(4.*k(2)-k(3))/3.
       psi(n)=0.
       ar(n)=0.
       k(n)=0.
c
c   solve the Hamiltonian constraint to find the metric component
c
       call calcphi(n,phi,du,psi,rs)
       do i=1,n
       grr(i)=phi(i)*phi(i)
       enddo
       return
       end
c
       subroutine evolve(n,psi,p,ar,k,phi,grr,dtpsi,dtp,dtar,dtk,
     & dtphi,dtgrr,du,dt,cvec,cnstr,alpha,trap,trap2,curv,mass,q,
     & cnstr2,tll,tnn,rs)
       implicit none
c
c   subroutine to calculate the time derivatives of the variables
c
       integer n
       real*8 psi(n),p(n),ar(n),k(n),phi(n),grr(n)
       real*8 dtpsi(n),dtp(n),dtar(n),dtk(n),dtphi(n),dtgrr(n)
       real*8 du,dt,rs
       integer i
       real*8 duk,duar,duphi,duuphi,dupsi,dup,dugrr
       real*8 drk,drar,drphi,drrphi,r,rad,drrad,drrrad
       real*8 drpsi,drp,drgrr,lappsi,diva
       real*8 up,um,dupsip,dupsim
       real*8 eps,arko(n),kko(n),psiko(n),pko(n),phiko(n),grrko(n)
       real*8 alpha(n),trap(n),trap2(n),q(n),cvec(4)
       real*8 c1,c2,c3,c4,c13,c14
       real*8 term1,term2,term3,term4,term5,term6,term7,term8,term9
       real*8 cnstr(n)
       real*8 cof1,cof2,cof3
       real*8 curv(n),mass(n)
       real*8 cnstr2(n),w(n),intw
       real*8 rho,js,s,mss,tll(n),tnn(n)
       real*8 u(n)
c
c   unpack the values of the c parameters
c
       c1=cvec(1)
       c2=cvec(2)
       c3=cvec(3)
       c4=cvec(4)
       c13=c1+c3
       c14=c1+c4 
       eps=0.5
       do i=1,n
       u(i)=du*float(i-1)
       enddo
c 
c   compute the Kreiss-Oliger dissipation terms
c   for extra stability
c
       do i=1,n
       arko(i)=0.
       kko(i)=0.
       psiko(i)=0.
       pko(i)=0.
       phiko(i)=0.
       grrko(i)=0.
       enddo
       do i=3,n-2
       arko(i)=(-1./16.)*(eps/dt)*(ar(i+2)+ar(i-2)+6.*ar(i)
     & -4.*(ar(i+1)+ar(i-1)))
       kko(i)=(-1./16.)*(eps/dt)*(k(i+2)+k(i-2)+6.*k(i)
     & -4.*(k(i+1)+k(i-1)))
       psiko(i)=(-1./16.)*(eps/dt)*(psi(i+2)+psi(i-2)+6.*psi(i)
     & -4.*(psi(i+1)+psi(i-1)))
       pko(i)=(-1./16.)*(eps/dt)*(p(i+2)+p(i-2)+6.*p(i)
     & -4.*(p(i+1)+p(i-1))) 
       phiko(i)=(-1./16.)*(eps/dt)*(phi(i+2)+phi(i-2)+6.*phi(i)
     & -4.*(phi(i+1)+phi(i-1))) 
       grrko(i)=(-1./16.)*(eps/dt)*(grr(i+2)+grr(i-2)+6.*grr(i)
     & -4.*(grr(i+1)+grr(i-1))) 
       enddo
c 
c   find the lapse, the trace-free part 
c   of the extrinsic curvature, and the shift
c
       call calcalpha(n,alpha,ar,du,u,rs)
       call calcq(n,q,phi,k,du,cvec,psi,p,u)
c
c   calculate the spatial derivatives of the variables
c
       do i=2,n-1
       duk=(k(i+1)-k(i-1))/(2.*du)
       duar=(ar(i+1)-ar(i-1))/(2.*du)
       duphi=(phi(i+1)-phi(i-1))/(2.*du)
       drk=duk*(1.-u(i))*(1.-u(i))/rs
       drar=duar*(1.-u(i))*(1.-u(i))/rs
       drphi=duphi*(1.-u(i))*(1.-u(i))/rs
       duuphi=(phi(i+1)+phi(i-1)-2.*phi(i))/du**2      
       up=u(i+1)
       um=u(i-1)
       term1=(1.5/(du*(up*up+up*um+um*um)))*
     & (up*up*ar(i+1)-um*um*ar(i-1))*((1.-u(i))**2)/rs
       term2=2.*(1.-u(i))*ar(i)/rs
       term3=2.*ar(i)*drphi/phi(i) 
       diva=term1+term2+term3
       r=rs*u(i)/(1.-u(i))
       rad=r*phi(i)
       drrad=phi(i)+u(i)*(1.-u(i))*duphi
       dupsi=(psi(i+1)-psi(i-1))/(2.*du)
       dup=(p(i+1)-p(i-1))/(2.*du)
       drpsi=dupsi*(1.-u(i))*(1.-u(i))/rs
       drp=dup*(1.-u(i))*(1.-u(i))/rs       
       dugrr=(grr(i+1)-grr(i-1))/(2.*du)
       drgrr=dugrr*(1.-u(i))*(1.-u(i))/rs
       drrrad=((1.-u(i))**3)*(u(i)*duuphi+2.*duphi)/rs
       up=u(i)+0.5*du
       um=u(i)-0.5*du
       dupsip=(psi(i+1)-psi(i))/du
       dupsim=(psi(i)-psi(i-1))/du
       term1=(3./(du*(up*up+up*um+um*um)))*
     & (up*up*dupsip-um*um*dupsim)*((1.-u(i))**4)/(rs*rs)
       term2=2.*drpsi*drphi/phi(i) 
       lappsi=term1+term2
c
c   calculate the time derivatives of the variables
c
       term1=alpha(i)*p(i)
       term2=psiko(i)
       dtpsi(i)=term1+term2
       term1=alpha(i)*(p(i)*k(i)+ar(i)*drpsi/grr(i))
       term2=alpha(i)*lappsi/grr(i)
       term3=-0.5*alpha(i)*drgrr*drpsi/(grr(i)*grr(i))
       term4=pko(i)
       dtp(i)=term1+term2+term3+term4
       cof1=c13/(c14*(1.-c13))
       cof2=-1.*(c2+c13)/(c14*(1.-c13))
       term1=alpha(i)*(k(i)/3.-2.*q(i))*ar(i)
       term2=alpha(i)*cof1*p(i)*drpsi
       term3=alpha(i)*cof2*drk
       term4=arko(i)
       dtar(i)=term1+term2+term3+term4
       cof1=(c14-2.)/(2.+c13+3.*c2)
       cof2=2./(2.+c13+3.*c2)
       cof3=3.*(1.-c13)/(2.+c13+3.*c2)
       term1=(alpha(i)/3.)*k(i)*k(i)
       term2=alpha(i)*cof1*(diva+ar(i)*ar(i))/grr(i)
       term3=-0.5*alpha(i)*cof1*drgrr*ar(i)/(grr(i)*grr(i))   
       term4=alpha(i)*cof2*p(i)*p(i)
       term5=alpha(i)*cof3*q(i)*q(i)
       term6=kko(i)
       dtk(i)=term1+term2+term3+term4+term5+term6
       term1=alpha(i)*phi(i)*(q(i)/2.-k(i)/3.)
       term2=phiko(i)
       dtphi(i)=term1+term2
       term1=-2.*alpha(i)*grr(i)*(q(i)+k(i)/3.)
       term2=grrko(i)
       dtgrr(i)=term1+term2
c
c   calculate the constraint to use as a code check
c
       term1=drrrad
       term2=(0.5/rad)*(drrad**2-grr(i))
       term3=(-0.5/grr(i))*drgrr*drrad
       term4=c14*rad*(0.5*drar+0.25*ar(i)**2)
       term5=c14*ar(i)*drrad
       term6=-0.25*(rad*c14/grr(i))*drgrr*ar(i)
       term7=(3./8.)*rad*grr(i)*(1.-c13)*(q(i)**2)  
       term8=(-1./12.)*rad*grr(i)*(2.+c13+3.*c2)*(k(i)**2)  
       term9=0.25*rad*(grr(i)*(p(i)**2)+drpsi**2)
       cnstr(i)=term1+term2+term3+term4+term5+term6+term7
     & +term8+term9
       w(i)=(2.*rad*drrad/grr(i))*(term4+term5+term6+term7
     & +term8+term9)
       cnstr(i)=cnstr(i)*rad/(10.+rad)
c 
c   calculate the mass and check for the presence 
c   of two kinds of trapped surface
c
       trap(i)=rad*(q(i)/2.-k(i)/3.)+drrad/sqrt(grr(i))
       trap2(i)=-1.*rad*(q(i)/2.-k(i)/3.)+drrad/sqrt(grr(i))
c       trap2(i)=rad(i)*(q(i)/2.-k(i)/3.)+drrad*sqrt(2.)
c       grr(i)=drrad**2-(rad(i)**2)*((q(i)/2.-k(i)/3.)**2)
       mass(i)=(rad/2.)*(1.-rad)
c
c   calculate components of the stress-energy
c
       term1=(0.5*c14/grr(i))*(2.*drar+ar(i)*(4.*drrad/rad+ar(i)))
       term2=0.5*(p(i)*p(i)+drpsi*drpsi/grr(i))
       term3=(-1./6.)*(c13+3.*c2)*k(i)*k(i)
       term4=-0.75*c13*q(i)*q(i)
       term5=-0.5*c14*ar(i)*drgrr/(grr(i)*grr(i))
       rho=term1+term2+term3+term4+term5
       js=(1./(1.-c13))*((c13+c2)*drk-p(i)*drpsi)/(sqrt(grr(i)))
       term1=2.*p(i)*p(i)-4.5*(c13+c2)*q(i)*q(i)
       term2=(c14+c13+2.*c2)*(drar+ar(i)*(2.*drrad/rad+ar(i)))
     &/grr(i)
       term3=-0.5*(c14+c13+2.*c2)*ar(i)*drgrr/(grr(i)*grr(i))
       s=-1.*rho+(2./(2.+c13+3.*c2))*(term1+term2+term3)
       term1=drpsi*drpsi+(c3-c4)*ar(i)*ar(i)
       term2=c13*(drar+ar(i)*drrad/rad)
       term3=(c13/rad)*(drrrad+(grr(i)-drrad*drrad)/rad)
       term4=(-0.5*drgrr/grr(i))*(c13*ar(i)+drrad/rad)
       mss=(2.*grr(i)/(3.*(1.-c13)))*(term1+term2+term3+term4)
       tll(i)=rho+(s/3.)+mss-2.*js
       tnn(i)=rho+(s/3.)+mss+2.*js
       enddo
       cnstr(1)=cnstr(2)
       cnstr(n-1)=cnstr(n-2)
       cnstr(n)=cnstr(n-1)
       trap(1)=1.
       trap(n)=trap(n-1)
       trap2(1)=sqrt(2.)
       trap2(n)=trap2(n-1)
       mass(1)=0.
       mass(n)=mass(n-1)
c       dtphi(n)=alpha(n)*phi(n)*(q(n)/2.-k(n)/3.)
c       dtgrr(n)=-2.*alpha(n)*grr(n)*(q(n)+k(n)/3.)
c       grr(1)=1.
c       grr(n)=grr(n-1)
c
c      calculate the second constraint
c
       cnstr2(1)=0.
       w(1)=0.
       intw=0.
       do i=2,n-2
       r=rs*u(i)/(1.-u(i))
       rad=r*phi(i)
       duphi=(phi(i+1)-phi(i-1))/(2.*du)
       drrad=phi(i)+u(i)*(1.-u(i))*duphi
       intw=intw+0.5*(du*rs/((1.-u(i))**2))*(w(i)+w(i-1))
       cnstr2(i)=rad*(1.-drrad*drrad/grr(i))-intw
c       cnstr2(i)=4.*cnstr2(i)
       enddo
       cnstr2(n-1)=0.
       cnstr2(n)=0.
c
c   find the squared curvature
c   (this blows up at a singularity)
c
c       call calccurv(n,curv,ar,k,q,rad,psi,p,dr,cvec)
       return
       end
c
       subroutine calccurv(n,curv,ar,k,q,rad,psi,p,dr,cvec)
       implicit none
c
c   subroutine to find the squared curvature
c
       integer n
       real*8 curv(n),ar(n),k(n),q(n),rad(n),psi(n),p(n),dr,cvec(4)
       integer i
       real*8 c1,c2,c3,c4,c13,c14
       real*8 y,r1,r2,k1,k2,l1,l2,m1,m2
       real*8 term1,term2,term3,term4,term5
       real*8 drk,drar,drrad,drpsi,drrrad 
       real*8 luk
       c1=cvec(1)
       c2=cvec(2)
       c3=cvec(3)
       c4=cvec(4)
       c13=c1+c3
       c14=c1+c4
       do i=2,n-1
       drk=(k(i+1)-k(i-1))/(2.*dr)
       drar=(ar(i+1)-ar(i-1))/(2.*dr)
       drrad=(rad(i+1)-rad(i-1))/(2.*dr)
       drpsi=(psi(i+1)-psi(i-1))/(2.*dr)
       drrrad=(rad(i+1)+rad(i-1)-2.*rad(i))
       y=(p(i)*drpsi-(c13+c2)*drk)/(1.-c13)
       r1=-2.*drrrad/rad(i)
       r2=(1.-(rad(i)*drrrad+drrad*drrad))/(rad(i)*rad(i))
       k1=q(i)+(k(i)/3.)
       k2=(k(i)/3.)-(q(i)/2.)
       l1=r1+2.*k1*k2
       l2=r2+k2*(k1+k2)
       term1=(k(i)*k(i))/3.
       term2=(c14-2.)*(drar+ar(i)*(ar(i)+2.*drrad/rad(i)))
       term3=2.*p(i)*p(i)
       term4=3.*(1.-c13)*q(i)*q(i)
       luk=term1+(term2+term3+term4)/(2.+c13+3.*c2) 
       term1=-1.*l2
       term2=0.5*(c1+c4)*(drar+ar(i)*(ar(i)+2.*drrad/rad(i)))  
       term3=0.5*(c13+c2)*(luk-k(i)*k(i))
       term4=c13*(k(i)*k2+ar(i)*drrad/rad(i))
       m2=(term1+term2+term3+term4)/(c13-1.)
       term1=(drar+ar(i)*(ar(i)+2.*drrad/rad(i)))
       term2=luk-2.*m2
       m1=term1+term2
       term1=3.*l1*l1+4.*l2*(l2-l1)
       term2=4.*m1*m1-8.*m2*m2
       term3=-4.*y*y
       curv(i)=term1+term2+term3
       enddo
       curv(1)=curv(2)
       curv(n)=curv(n-1)
       return
       end
c
       subroutine calcalpha(n,alpha,ar,du,u,rs)
       implicit none
c
c   subroutine to find the lapse by integrating
c   the acceleration of the aether field
c
       integer n
       real*8 alpha(n),ar(n),du,u(n),rs
       real*8 lnalpha(n),uavg,aravg
       integer i
       lnalpha(n)=0.
       do i=n,2,-1
       aravg=0.5*(ar(i)+ar(i-1))
       uavg=0.5*(u(i)+u(i-1))
       lnalpha(i-1)=lnalpha(i)-du*rs*aravg/(1.-uavg)**2
       enddo
       do i=1,n
       alpha(i)=exp(lnalpha(i))
       enddo 
       return
       end
c
       subroutine calcq(n,q,phi,k,du,cvec,psi,p,u)
       implicit none
c
c   subroutine to find the trace-free part of the 
c   extrinsic curvature by integrating the momentum constraint
c
       integer n
       real*8 q(n),phi(n),k(n),du,cvec(4),psi(n),p(n),u(n)
       real*8 c2,c13
       real*8 uavg,phiavg,duphi,duk,qint
       integer i
       real*8 cof,kavg
       real*8 dupsi,pavg
       real*8 term1,term2
       c2=cvec(2)
       c13=cvec(1)+cvec(3)
       cof=(2.+c13+3.*c2)/(3.-3.*c13)
       q(1)=0.
       do i=1,n-1
       uavg=du*(float(i)-0.5)
       phiavg=0.5*(phi(i)+phi(i+1))
       kavg=0.5*(k(i)+k(i+1))
       pavg=0.5*(p(i)+p(i+1))
       duphi=(phi(i+1)-phi(i))/du
       duk=(k(i+1)-k(i))/du
       dupsi=(psi(i+1)-psi(i))/du
       term1=-3.*q(i)*(1./(uavg*(1.-uavg))+duphi/phiavg)
       term2=cof*duk-pavg*dupsi/(1.-c13) 
       qint=q(i)+0.5*du*(term1+term2)
       term1=-3.*qint*(1./(uavg*(1.-uavg))+duphi/phiavg)
       q(i+1)=q(i)+du*(term1+term2)
       enddo
       return
       end
c
       subroutine calcphi(n,phi,du,psi,rs)
       implicit none
       integer n,i
       real*8 phi(n),du,psi(n),rs
       real*8 w,wint,phiint,u,dupsi,term1,term2,term3
c
       phi(1)=1.
       w=1.
       do i=1,n-1
       u=du*(float(i)-0.5)
       dupsi=(psi(i+1)-psi(i))/du
       phiint=phi(i)+0.5*du*(w-phi(i))/(u*(1.-u))
       term1=(w-phi(i))*(w-phi(i))/(2.*phi(i)*u*(1.-u))
       term2=-0.25*phi(i)*u*(1.-u)*dupsi*dupsi
       wint=w+0.5*du*(term1+term2)
       phi(i+1)=phi(i)+du*(wint-phiint)/(u*(1.-u))
       term1=(wint-phiint)*(wint-phiint)/(2.*phiint*u*(1.-u))
       term2=-0.25*phiint*u*(1.-u)*dupsi*dupsi
       w=w+du*(term1+term2)
       enddo
       return
       end
