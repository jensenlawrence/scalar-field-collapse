c      program to evolve the Einstein-scalar equations
c      in spherical symmetry
c      with maximal slicing and zero shift
c
       implicit none
       integer n,nvars,nvaux
       parameter(n=9601,nvars=3,nvaux=3)
       real*8 vars(n,nvars),newvars(n,nvars)
       real*8 dtvars(n,nvars),dtnewvars(n,nvars)
       real*8 r(n),dr
       integer i,j,itime,ntime
       real*8 time,dt
       integer iter
       integer ich,ipr
       real*8 phi(n),p(n),s(n)
       real*8 amp,sigma
       real*8 rmax,r0
       real*8 auxvars(n,nvaux),newauxvars(n,nvaux)
       real*8 alpha(n),a(n),krr(n)
       real*8 alphanew(n),anew(n),krrnew(n)
       real*8 cnstr(n),adot,avgalpha,avgkrr
       real*8 cnstrl2
       real*8 cnstr2(n)
       real*8 tmax
       real*8 trap,trapmin,mass,gam
       open(unit=14,file='lncnstrl2.dat',status='unknown')
       open(unit=15,file='cnstrfin.dat',status='unknown')
       open(unit=16,file='cnstr2fin.dat',status='unknown')
       open(unit=18,file='phifin.dat',status='unknown')
       open(unit=19,file='trapmin.dat',status='unknown')
       open(unit=20,file='phi0.dat',status='unknown')
       open(unit=21,file='alphafin.dat',status='unknown')
       open(unit=23,file='pfin.dat',status='unknown')
       open(unit=24,file='gammafin.dat',status='unknown')
c       ntime=2000
       ntime=20000
       amp=0.15
       sigma=4.
       rmax=80.
       r0=10. 
       tmax=19.
       call init(n,nvars,r,vars,dr,dt,amp,sigma,rmax,r0)
       time=0.
       do itime=1,ntime
       call evolve(n,nvars,vars,dtvars,dr,r,dt,nvaux,auxvars)
       do i=1,n
       do j=1,nvars
       newvars(i,j)=vars(i,j)
       enddo
       enddo
       do iter=1,3
       call evolve(n,nvars,newvars,dtnewvars,dr,r,dt,
     & nvaux,newauxvars)
       do i=2,n-1
       do j=1,nvars
       newvars(i,j)=vars(i,j)+dt*0.5*(dtvars(i,j)+dtnewvars(i,j))
       enddo
       enddo
       newvars(1,1)=(4.*newvars(2,1)-newvars(3,1))/3.
       newvars(1,2)=(4.*newvars(2,2)-newvars(3,2))/3.
       newvars(1,3)=0.
       newvars(n,1)=0.
       newvars(n,2)=0.       
       newvars(n,3)=0.
       enddo
c
c      unpack variables
c
       do i=1,n
       phi(i)=vars(i,1)
       p(i)=vars(i,2)
       s(i)=vars(i,3)
       alpha(i)=auxvars(i,1)
       a(i)=auxvars(i,2)
       krr(i)=auxvars(i,3)
       alphanew(i)=newauxvars(i,1)
       anew(i)=newauxvars(i,2)
       krrnew(i)=newauxvars(i,3)
       enddo
c
c      write out the data
c
       ipr=itime/100
       ich=ipr*100
       if(ich.eq.itime) then
       write(*,*) itime,time
       call xvs('Phi',time,r,phi,n)
       call xvs('alpha',time,r,alpha,n)
       call xvs('P',time,r,p,n)
       call xvs('S',time,r,s,n)
       call xvs('A',time,r,a,n)
       call xvs('Krr',time,r,krr,n)
       endif
c      calculate constraint
       do i=2,n-1
       adot=(anew(i)-a(i))/dt
       avgalpha=0.5*(alpha(i)+alphanew(i))
       avgkrr=0.5*(krr(i)+krrnew(i))
       cnstr(i)=adot-avgalpha*avgkrr
       cnstr2(i)=(phi(i+1)-phi(i-1))/(2.*dr)-s(i)
       cnstr(i)=4.*cnstr(i)
c       cnstr2(i)=4.*cnstr2(i)
       enddo
       cnstr(1)=cnstr(2)
       cnstr(n)=cnstr(n-1)
       cnstr2(1)=cnstr2(2)
       cnstr2(n)=cnstr2(n-1)
       cnstrl2=0.5*dr*cnstr(1)*cnstr(1)
       do i=2,n-1
       cnstrl2=cnstrl2+dr*cnstr(i)*cnstr(i)
       enddo
       cnstrl2=cnstrl2+0.5*dr*cnstr(n)*cnstr(n)
       cnstrl2=sqrt(cnstrl2)
       write(14,22) time,log(cnstrl2)
  22   format(1x,2e18.6)
c      update variables
       do i=1,n
       do j=1,nvars
       vars(i,j)=newvars(i,j)
       enddo
       enddo  
       time=time+dt
       write(20,22) time,phi(1)
       if(time.gt.tmax) goto 88
       trapmin=1.
       do i=2,n-1
       trap=1.+0.5*r(i)*((0.5/dr)*(a(i+1)-a(i-1))
     & +krr(i)*exp(-1.*a(i)))
       if(trap.lt.trapmin) then
       trapmin=trap
       endif
c       if(trap.le.0.) then
c       mass=0.5*r(i)*exp(0.5*a(i))
c       write(*,29) mass,time
c  29   format(1x,'a black hole with mass',2x,f12.7,2x,
c     & 'forms at time',2x,f12.7)
c       stop
c       endif
       enddo
       write(19,22) time,trapmin
       enddo
  88   continue
       write(*,31) 
  31   format(1x,'reached maximum time')
       do i=1,n/4
       write(15,22) r(i),cnstr(i)
       write(16,22) r(i),cnstr2(i)
       write(18,22) r(i),phi(i)
       write(21,22) r(i),alpha(i)
       write(23,22) r(i),p(i)
       gam=exp(-2.*a(i))
       write(24,22) r(i),gam
       enddo
       end
c
       subroutine init(n,nvars,r,vars,dr,dt,amp,sigma,rmax,r0)
       implicit none
       integer n,nvars,i
       real*8 r(n),vars(n,nvars),dr,dt,amp,sigma,rmax,r0
       real*8 rarg1,rarg2,exp1,exp2,term1,term2,phi,p,s,exparg
       dr=rmax/float(n-1)
       dt=0.5*dr
       do i=1,n
       r(i)=dr*float(i-1)
       enddo
       do i=2,n
       phi=0.
       p=0.
       exparg=((r(i)**2-r0**2)/sigma**2)**2
       if(exparg.lt.50.) then    
       phi=amp*exp(-1.*exparg)
       s=(phi/sigma**4)*(-4.*r(i)*(r(i)**2-r0**2))
       endif
       vars(i,1)=phi
       vars(i,2)=p
       vars(i,3)=s 
       enddo
       vars(1,1)=(4.*vars(2,1)-vars(3,1))/3
       vars(1,2)=(4.*vars(2,2)-vars(3,2))/3.
       vars(1,3)=0.
       return
       end
c
       subroutine evolve(n,nvars,vars,dtvars,dr,r,dt,
     & nvaux,auxvars)
       implicit none
       integer n,nvars,nvaux
       real*8 vars(n,nvars),dtvars(n,nvars),dr,r(n),dt
       real*8 auxvars(n,nvaux)
       real*8 phi(n),p(n),s(n),dtphi,dtp,dts
       integer i
       real*8 a(n),krr(n),mass(n)
       real*8 drp,drs,dra,dralpha
       real*8 eps,phiko(n),pko(n),sko(n)
       real*8 alpha(n)
       real*8 term1,term2,term3,term4
       real*8 rp,rm
       real*8 trap1(n)
c       eps=0.
       do i=1,n
       phi(i)=vars(i,1)
       p(i)=vars(i,2)
       s(i)=vars(i,3)
       enddo
       eps=0.5
       do i=1,n
       phiko(i)=0.
       pko(i)=0.
       sko(i)=0.
       enddo
       do i=3,n-2
       phiko(i)=(-1./16.)*(eps/dt)*(phi(i+2)+phi(i-2)+6.*phi(i)
     & -4.*(phi(i+1)+phi(i-1)))
       pko(i)=(-1./16.)*(eps/dt)*(p(i+2)+p(i-2)+6.*p(i)
     & -4.*(p(i+1)+p(i-1)))
       sko(i)=(-1./16.)*(eps/dt)*(s(i+2)+s(i-2)+6.*s(i)
     & -4.*(s(i+1)+s(i-1)))
       enddo
       call calcka(n,a,krr,p,s,dr,r)
       call calcalpha(n,alpha,a,krr,p,s,dr,r)
       do i=2,n-1
       dtphi=alpha(i)*p(i)+phiko(i)
       drp=(p(i+1)-p(i-1))/(2.*dr)
       drs=(s(i+1)-s(i-1))/(2.*dr)
       dra=(a(i+1)-a(i-1))/(2.*dr)
       dralpha=(alpha(i+1)-alpha(i-1))/(2.*dr)
       term1=dralpha*p(i)*exp(a(i))
       term2=alpha(i)*(drp*exp(a(i))+krr(i)*s(i))
       term3=sko(i)
       dts=term1+term2+term3
       term1=dralpha*s(i)*exp(a(i))
       rp=r(i+1)
       rm=r(i-1)
       term2=alpha(i)*(1.5/(dr*(rp*rp+rp*rm+rm*rm)))*
     & (rp*rp*s(i+1)*exp(a(i+1))-rm*rm*s(i-1)*exp(a(i-1)))
       term3=pko(i)
       dtp=term1+term2+term3
       term1=1.+exp(a(i))*(0.5*r(i)*krr(i))**2
       term2=-1.*exp(3.*a(i))*(1.+0.5*r(i)*dra)**2
       mass(i)=0.5*r(i)*exp(0.5*a(i))*(term1+term2)
       trap1(i)=1.+0.5*r(i)*(dra+exp(-1.*a(i))*krr(i))
       dtvars(i,1)=dtphi
       dtvars(i,2)=dtp
       dtvars(i,3)=dts
       enddo
       mass(1)=0.
       mass(n)=mass(n-1)
       trap1(1)=trap1(2)
       trap1(n)=trap1(n-1)
c
c      pack up the auxiliary variables
c
       do i=1,n
       auxvars(i,1)=alpha(i)
       auxvars(i,2)=a(i)
       auxvars(i,3)=krr(i)
       enddo
       return
       end
c
       subroutine calcka(n,a,krr,p,s,dr,r)
       implicit none
       integer n
       real*8 a(n),krr(n),p(n),s(n),dr,r(n)
       integer i
       real*8 k1(3),k2(3)
       real*8 term1,term2,term3,term4
       real*8 w
       real*8 atmp,wtmp,krrtmp,ptmp,stmp,rtmp
       real*8 aa,bb
       aa=(-1./30.)*(p(1)*p(1))
       bb=-0.2*p(1)*(s(2)/dr)
       do i=1,4
       a(i)=aa*r(i)*r(i)
       krr(i)=bb*r(i)*r(i)
       w=2.*aa*r(i)
       enddo
       do i=4,n-1
       k1(1)=dr*w
       term1=(exp(-3.*a(i))-1.)/(r(i)*r(i))
       term2=-1.*w*((5./r(i))+(7./4.)*w)
       term3=-0.5*exp(-2.*a(i))*(1.5*krr(i)**2+p(i)**2+s(i)**2)
       k1(2)=dr*(term1+term2+term3)
       term1=-1.*krr(i)*(1.5*w+3./r(i))
       term2=-1.*exp(-1.*a(i))*p(i)*s(i)
       k1(3)=dr*(term1+term2)
       atmp=a(i)+0.5*k1(1)
       wtmp=w+0.5*k1(2)  
       krrtmp=krr(i)+0.5*k1(3)
       ptmp=0.5*(p(i)+p(i+1))
       stmp=0.5*(s(i)+s(i+1))
       rtmp=0.5*(r(i)+r(i+1))
       k2(1)=dr*wtmp
       term1=(exp(-3.*atmp)-1.)/(rtmp*rtmp)
       term2=-1.*wtmp*((5./rtmp)+(7./4.)*wtmp)
       term3=-0.5*exp(-2.*atmp)*(1.5*krrtmp**2+ptmp**2+stmp**2)
       k2(2)=dr*(term1+term2+term3)
       term1=-1.*krrtmp*(1.5*wtmp+3./rtmp)
       term2=-1.*exp(-1.*atmp)*ptmp*stmp
       k2(3)=dr*(term1+term2)
       a(i+1)=a(i)+k2(1)
       w=w+k2(2)
       krr(i+1)=krr(i)+k2(3)   
       enddo
       return
       end
c
       subroutine calcalpha(n,alpha,aa,krr,p,s,dr,rr)
       implicit none
       integer n
       real*8 alpha(n),aa(n),krr(n),p(n),s(n),dr,rr(n)
       integer i
       real*8 capa(n),capb(n),a(n-2),b(n-2),c(n-2),r(n-2),u(n-2)
       real*8 ap
       do i=2,n-1
       ap=(aa(i+1)-aa(i-1))/(2.*dr)
       capa(i)=2./rr(i)+2.*ap
       capb(i)=exp(-2.*aa(i))*(1.5*(krr(i)**2)+(p(i)**2))
       enddo
       do i=1,n-2
       a(i)=0.5*capa(i+1)*dr-1.
       b(i)=2.+capb(i+1)*dr*dr
       c(i)=-1.*(1.+0.5*capa(i+1)*dr)
       r(i)=0.
       enddo
       b(1)=(2./3.)*(1.+dr*capa(2))+capb(2)*dr*dr
       c(1)=(-2./3.)*(1.+capa(2)*dr)
       r(n-2)=1.+0.5*capa(n-1)*dr
       call tridag(a,b,c,r,u,n-2)
       do i=1,n-2
       alpha(i+1)=u(i)
       enddo
       alpha(1)=(4.*alpha(2)-alpha(3))/3.
       alpha(n)=1.
       return
       end
c
       subroutine tridag(a,b,c,r,u,n)
       implicit none
       integer n
       real*8 a(n),b(n),c(n),r(n),u(n)
       real*8 bet,gam(n)
       integer j
       bet=b(1)
       u(1)=r(1)/bet
       do j=2,n
       gam(j)=c(j-1)/bet
       bet=b(j)-a(j)*gam(j)
       u(j)=(r(j)-a(j)*u(j-1))/bet
       enddo
       do j=n-1,1,-1
       u(j)=u(j)-gam(j+1)*u(j+1)
       enddo
       return 
       end
