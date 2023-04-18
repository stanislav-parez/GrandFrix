module mod_principalstress
    use mycommons 
    implicit none

    private
    
    public :: PRINCIPALSTRESS, STRESS

    contains

!PRINCIPALSTRESS*******************************************
!This subroutine calculates the principal stresses at the *
!center of the grains                                     *
!                                                         *
!Variables used :                                         *
!u1,u2,u3 : arrays calculated in subroutine "stress"      *
!b,c : coefficients of the quadratic equation             *
!sigxx,sigyy,sigxy : coefficients of the stress matrix at *
!the center of the grains                                 *
!l1,l2 : eigenvalues of the stress tensor                 *
!cron : angle from x-axis to max principal stress sigmax  *
!sigmax,sigmin : max and min principal stresses at the    *
!center of the grains                                     *
!                                                         *
!Input : u1,u2,u3                                         *
!                                                         *
!Output : sigmax,sigmin,l1,l2,cron                        *
!                                                         *
!Programmer : Richard A. Lang                             *
!**********************************************************
!=======================================================================
SUBROUTINE PRINCIPALSTRESS()

    integer k
    
    real(8) q,f
    real(8) sigxx(MAXN),sigyy(MAXN),sigxy(MAXN)
    real(8) b(MAXN),c(MAXN)
    

    f=0.
    q=0.
    
    do k=1,n
        sigxx(k)=0.
        sigyy(k)=0.
        sigxy(k)=0.
        sigmax(k)=0.
        sigmin(k)=0.
        cron(k)=0.
        l1(k)=0.
        l2(k)=0.
        b(k)=0.
        c(k)=0.
    enddo


    do k=1,n
        sigxx(k)= 0.5*(u1(k)-u2(k))
        sigyy(k)= 0.5*(u1(k)+u2(k))
        sigxy(k)= 0.5*u3(k)
    enddo
    
    !    write(*,*) 'U1:'
    !       write(*,*) (u1(k),k=1,20)

    do k=1,n
        if (coordn(k) >= 2) then
            !         print*, coordn(k)
            b(k)= -(sigxx(k)+sigyy(k))
            c(k)= -((sigxy(k)*sigxy(k))-(sigxx(k)*sigyy(k)))
            f= sqrt((b(k)*b(k))-(4.0*c(k)))
            f=sign(f,b(k))
            q=-0.5*(b(k)+f)
            l1(k)=q
            l2(k)=c(k)/q

            if (l1(k) >= l2(k)) then
                sigmax(k)=l1(k)
                sigmin(k)=l2(k)
            else
                sigmax(k)=l2(k)
                sigmin(k)=l1(k)
            endif

            if (sigxy(k) == 0.) then
                cron(k)=0.5*PI
            else
                    cron(k)=0.5*acos((2.*sigxx(k)-(sigmin(k)+sigmax(k)))/ &
                            (sigmin(k)-sigmax(k)))
            endif
        else
            sigmax(k)= 0.
            sigmin(k)=0.
            cron(k)=0.
        endif

    enddo

end subroutine
!STRESS****************************************
!This subroutine calculates the quantities u1,*
!u2,u3 whenever two particles m,i interact    *
!Such quantites are used to calculate the     *
!stress tensor at the center of the particles *
!in the subroutine Principalstress            *
!                                             *
!Variables used:                              *
!m,i : particles indices                      *
!xm,ym,xi,yi : x and y components of forces   *
!radius : radius of particles                 *
!eta,beta,omega,zeta,mu : angles between      *
!interparticle forces and x-axis              *
!dm,em,fm,gm : sin and cos for particle m     *
!di,ei,fi,gi : sin and cos for particle i
!u1,u2,u3 : arrays used for the calculation   *
!of stress tensor                             *
!                                             *
!Input : m,i,xm,ym,xi,yi,radius,u1,u2,u3      *
!                                             *
!Output : Updated u1,u2,u3                    *
!                                             *  
!Programmer : Richard A. Lang                 *
!**********************************************
!=======================================================================
SUBROUTINE STRESS(m,i,rx1,rx2,xm,ym,xi,yi)

    integer m,i
    
    real(8) xm,ym,xi,yi,eta,beta,omega,zeta,mu
    real(8) dm,em,fm,gm,di,ei,fi,gi,rx1,rx2
    real(8) sino,coso
    real(8) d1
    
    d1=180./pi
    eta=0.
    beta=0.
    mu=0.
    omega=0.
    zeta=0.
    dm=0.
    em=0.
    fm=0.
    gm=0.
    di=0.
    ei=0.
    fi=0.
    gi=0. 

    if (rx1 == rx2) then
        if (ry(m) < ry(i)) then
            dm=1.
            em=0.
            di=-1.
            ei=0.
        else
            dm=-1.
            em=0.
            di=1.
            ei=0.
        endif
    elseif (ry(m) == ry(i)) then
        if (rx1 < rx2) then
            dm=0.
            em=1.
            di=0.
            ei=-1.
        else
            dm=0.
            em=-1.
            di=0.
            ei=1.
        endif
    else

        omega=atan(abs((ry(m)-ry(i))/(rx1-rx2)))
        sino=sin(omega)
        coso=cos(omega)

        ! where is i relative to m?
        if (rx1 < rx2) then
            if (ry(m) < ry(i)) then
                ! upper right quadrant
                dm=sino
                em=coso
                di=-sino
                ei=-coso
            else
                ! lower right
                dm=-sino
                em=coso
                di=sino
                ei=-coso
            endif
        else
            if (ry(m) < ry(i)) then
                ! upper left
                dm=sino
                em=-coso
                di=-sino
                ei=coso
            else
                ! lower left
                dm=-sino
                em=-coso
                di=sino
                ei=coso
            endif
        endif
    endif

    fm=dm*(3.-4.*dm*dm)
    gm=em*(4.*em*em-3.)
    fi=di*(3.-4.*di*di)
    gi=ei*(4.*ei*ei-3.)

    u1(m)=u1(m)+(xm*em+ym*dm)/(PI*radius(m))
    u2(m)=u2(m)-(xm*(em+gm)+ym*(fm-dm)) / (PI*radius(m))
    u3(m)=u3(m)+(xm*(dm+fm)+ym*(em-gm)) / (PI*radius(m))

    u1(i)=u1(i)+(xi*ei+yi*di)/(PI*radius(i))
    u2(i)=u2(i)-(xi*(ei+gi)+yi*(fi-di)) / (PI*radius(i))
    u3(i)=u3(i)+(xi*(di+fi)+yi*(ei-gi))/ (PI*radius(i))

end subroutine

end module
