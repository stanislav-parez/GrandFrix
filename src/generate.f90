module mod_generate
    use mycommons
    use mod_inputReader
                               
    implicit none

    private

    public :: addgrain2, gauss, shearstep, partgen, STRENGTH, rand1

    contains
!=======================================================================
subroutine partgen(phigoal,boxvol,grainarea,alap1,gfac)
    ! generate a randomly-filled box of particles, starting in the lower
    ! left corner and adding particles along the bottom, then starting
    ! a new row. Put the particles as close together as possible without
    ! overlapping them. By putting particles closer (reduce the factor
    ! gfac), reduces initial pore space but there will be overlap.)
    
    integer nlastrow,l,krow,nrow,nrowmax,np
    integer i
    integer dflag
    integer kl0,kl1
    
    real(8) rmaxl,rmaxr
    real(8) gfac
    real(8) yedge(maxn),yedgepos(maxn),ynewpos(maxn),ynewedge(maxn)
    real(8) rad0,rxg,grad0,ryg,xedge
    real(8) yt,solidvol,xr,alap,boxvol
    real(8) grainarea,phigoal,alap1,newarea
    real(8) highest
    real(8) dist,deviate,MyValue
    real(8) MySum,averad

    grainarea=0d0

    print*, 'generating random grid of grains...'
    print*, 'imposed overlap =', alap1

    !  add grains, until their volume reaches the desired porosity
    !  at the eventual box size
    solidvol = (1.-phigoal)*boxvol

    !nlastrow=nbound(1)

    ! create a surface that drapes over the last row of particles, defined
    !  by the points(yedgepos,yedge) in (x,y)
    ! add new grains on top of this surface, starting at left wall and working
    ! toward right wall. Then drape a new surface over the top of this and
    ! start over. 
    !  rad0=atom radius in this grain
    !  grad0=effective grain radius(some factor*rad0)
    !  nlastrow = number of points that make up the underlying surface
    !    (given by the number of grains placed into the previous row)
    !  krow = counter for number of grains in this row
    !  nrow = counter for number of rows

    if (ib(3) >= 0) then
        rmaxl=radius(nbound(3))
        do i=nbound(2)+1,nbound(3)-1
            rmaxl=max(rmaxl,radius(i))
        enddo
        xleft = xleft_out + w_offset + rmaxl
        rmaxr=radius(nbound(4))
        do i=nbound(3)+1,nbound(4)-1
            rmaxr=max(rmaxr,radius(i))
        enddo
        xright = xright_out - w_offset - rmaxr
        ytop = ytop_out - w_offset
        ybot = ybot_out + w_offset
    else
        xleft = xleft_out
        xright = xright_out
        ytop = ytop_out
        ybot = ybot_out
    endif

    xr = xright-0.5
    yt = ytop-0.5

!    do concurrent (l = 1:nlastrow)
!        yedge(l) = ry(l)+radius(l)
!        yedgepos(l) = rx(l)
!    enddo
    nlastrow=int(xr-xleft)
    do concurrent (l = 1:nlastrow)
        yedgepos(l) = xleft + l - 1
    enddo
    highest = maxval(radius(1:nbound(1)))
    yedge(1:nlastrow) = ybot + highest
    
    yedge(nlastrow+1)=yedge(1)
    yedgepos(nlastrow+1)=xr
    xedge=xleft

    krow=0
    nrow=1
    nrowmax=300
    kl0=1

    do np=1,maxn

        ! find grain size from set of gaussian distributions
        dist=rand1(randum)
        
        if (dist <= frac1) then 
            do
                deviate=gauss(randum)
                deviate=std1*deviate
                MyValue=mean1 + deviate
                dflag=1 
                if (.not.(MyValue < logr .or. MyValue > higr)) exit
            enddo
        else 
            do
                deviate=gauss(randum)
                deviate=std2*deviate
                MyValue=mean2 + deviate 
                if(.not.(MyValue < logr2 .or. MyValue > higr2)) exit
            enddo
            dflag=2
        endif
        rad0=0.5*MyValue
       
        alap=alap1
        newarea = grainarea + rad0*rad0*(pi-alap)
        
        if (grainarea>solidvol) exit

        grad0=gfac*rad0
        rxg=xedge+grad0
        xedge=rxg+grad0          

        !  make a new row
        if (xedge > xright) then
            nlastrow=krow
            krow=0
            nrow=nrow+1
            kl0=1
            rxg=xleft+grad0
            xedge = rxg+grad0
            do concurrent (l = 1:nlastrow)
                yedgepos(l) = ynewpos(l)
                yedge(l) = ynewedge(l)
            enddo
            yedge(nlastrow+1)=yedge(nlastrow)
            yedgepos(nlastrow+1)=xr
            if (nrow > nrowmax) exit
        endif

        !           loop over lastrow to find height of pile
        do
            krow=krow+1
            do l=2,nlastrow+1
                if ((rxg-grad0) < yedgepos(l)) then
                    kl0=l-1
                    exit
                endif
            enddo
            
            kl1=kl0-1
            do l=kl0,nlastrow+1
                if ((rxg+grad0) < yedgepos(l)) then
                    kl1=l
                    exit
                endif
            enddo
            
            highest=yedge(kl0)
            if(size(yedge(kl0+1:kl1))>0) highest = maxval(yedge(kl0+1:kl1))
            ryg = highest + grad0

            ynewedge(krow)=ryg+grad0
            ynewpos(krow)=rxg

            !  if this particle runs into the top boundary, quit adding particles
            if (ynewedge(krow) > (yt)) then
                rxg=rxg+grad0
                xedge=rxg+grad0
                
                if (xedge > xr) then
                    nrowmax = nrow
                    print*, 'box too short, increase boxy'
                    stop
                else
                    cycle
                endif
            endif
            
            exit
        enddo
        
        
        rxg=rxg+(ryg-ybot)*.001
        call addgrain2(rxg,ryg,rad0,0)
        grainarea=newarea
        gtype(n)=dflag
        
        if (dflag == 1) then 
            emod(n)=e1
        else
            emod(n)=e1
        endif

    enddo

    print*, 'grainarea, solidvol',grainarea, solidvol

    MySum = 0.
    do concurrent (i = 1:n)
        MySum=MySum+radius(i)
    enddo
    
    averad=MySum/dble(n)
    print*, "average particle diameter =", 2*averad

end subroutine
!=======================================================================
subroutine addgrain2(rxg,ryg,rad0,isbound)
    
    integer isbound
    
    real(8),parameter :: pi2 = 2.*pi
    real(8) rxg,ryg,rad0

    n=n+1

    if (n > maxn) then
        print*, 'TOO MANY GRAINS!'
        stop
    endif

    bdgrain(n) = isbound

    ! set the positions of the particles
    rx(n)=rxg
    ry(n)=ryg
    radius(n) = rad0
    radinv2(n) = .125/(radius(n)*radius(n)*radius(n))
    rnt(n)=0.
    call strength(n)
    vx(n)=0.
    vy(n)=0.
    w(n)=0.
    gtype(n)=1

end subroutine
!=======================================================================
function gauss(idum)

    integer idum
    integer,save :: iset=0
    
    real(8) gauss
    real(8) fac, rsq,v1,v2
    real(8),save :: gset
    
    if (iset == 0) then
        do
            v1 = 2*rand1(idum)-1.
            v2 = 2*rand1(idum)-1.
            rsq = v1*v1+v2*v2
            if(.not.(rsq >= 1. .or. rsq == 0.)) exit
        end do
        fac = sqrt(-2*log(rsq)/rsq)
        gset = v1*fac
        gauss = v2*fac
        iset = 1
    else
        gauss = gset
        iset = 0
    endif
    
end function
!=======================================================================
subroutine shearstep(gstep)        

    integer i
    
    real(8) dxs,gstep

    do concurrent (i = 1:n)
        dxs = gstep*(ry(i)-ybot) 
        rx(i) = rx(i)+dxs
    enddo

end subroutine
!======================================================================= 
!STRENGTH*********************************************
!This subroutine assignes strength to the grain i    *
!according to a switch. So, if                       *
!typest=0                                            *
!all the grains have the same strength regardless    *
!regardless thier size                               *
!                                                    *
!typest=1                                            *
!the grains have a size-dependent strength and a     *
!random component                                    *
!                                                    *
!typest=2                                            *
!the grains have only a size-dependent strength      *
!                                                    *
!typest=3                                            *
!the grains have only a random strength              *
!                                                    *
!Variables used:                                     *
!i : grain index                                     *
!randnm : random number between 0 and 1              *
!A,B : coefficients                                  *
!typest : switch for type of strength                *
!purest : min. possibile strength                    *
!esp : esponent used for the strength-size dependence*
!strgth : strength of the grains                     *
!                                                    *
!Input : i,randnm,typest,esp,purest,radius           *
!Output : strgth                                     * 
!Programmer: Richard A. Lang                         *
!*****************************************************
!SUBROUTINE STRENGTH(i,randnm)
subroutine STRENGTH(i)

    integer  i
    
    strgth(i)=purest/((2.0*radius(i))**esp)
       
end subroutine
!=======================================================================
real(8) function rand1(idum)

    integer,parameter :: ia=16807
    integer,parameter :: im=2147483647
    integer,parameter :: iq=127773
    integer,parameter :: ir=2836
    integer,parameter :: ntab=32
    integer,parameter :: ndiv=1 + int(dble(im-1) / dble(ntab))
    integer idum
    integer j,k
    integer,save :: iv(ntab) = 0
    integer,save :: iy = 0
    
    real(8),parameter    :: am=1./dble(im)
    real(8),parameter    :: eps=1.2e-7
    real(8),parameter    :: rnmx=1.-eps

    if (idum <= 0 .or. iy == 0) then
        idum=max(-idum,1)
        do j=ntab+8,1,-1
            k=idum/iq
            idum=ia*(idum-k*iq)-ir*k
            if (idum < 0) idum=idum+im
            
            if (j <= ntab) iv(j)=idum
        enddo
        iy=iv(1)
    endif

    k=idum/iq
    idum=ia*(idum-k*iq)-ir*k
    if (idum < 0) idum=idum+im
    
    j=1+iy/ndiv
    iy=iv(j)
    iv(j)=idum
    rand1=min(dble(am*iy),rnmx)
    
end function

end module
