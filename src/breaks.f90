module mod_breaks
    use mycommons
    use mod_generate, only : addgrain2, STRENGTH
    use mod_linklist, only : links
    use mod_principalstress, only : PRINCIPALSTRESS
    use mod_findforces, only : getneighbors
    use mod_inputReader
    implicit none

    private
    
    public :: breaking

    contains
    
!=======================================================================
subroutine breaking(rsmall,ibreak,tresume,grainarea,icrush,inext)

    
    integer bcontact(325),icrush
    integer ibreak
    integer i,ik,j,k,ibpart,nevent
    integer ibknt,ig,itf
    integer inext
    
    real(8) f2(13),breakmin
    real(8) strlimit,fmax(maxn),olap,olap2
    real(8) rsmall,tresume
    real(8) ff,masskept,grainarea,rnew,distint,tol,maxstr
    
    199  format(a,f6.3,a,e11.3)
    
    ibreak=0
    
! Routine to break grains. This has three steps
!  1. Check principal stresses against grain strengths
!  2. If grain is breakable, replace with new grains
!  3. Do series of fake time steps to find low-stress configuration of grains

! Step 1. 
    teindx(1:n)=0.
    stot(1:n)=0.
    fmax(1:n)=0.
    
    do ik=1,contactknt
    
        i=contacti(ik)
        j=contactj(ik)
        olap=sqrt(contfn(ik)*contfn(ik) + contft(ik)*contft(ik))
        olap2=olap*0.999
        stot(i)=stot(i)+contfn(ik)
        stot(j)=stot(j)+contfn(ik)
        
        if (radius(i) > radius(j)) then
            fmax(i)=max(fmax(i),olap)
            fmax(j)=max(fmax(j),olap2)
        else
            fmax(i)=max(fmax(i),olap2)
            fmax(j)=max(fmax(j),olap)
        endif
    enddo

    if (icrush == 1) call principalstress()

    maxstr=0.
    ibpart=0

    ! this checks that only interior grains that are under a total compressive
    !  force that is LARGE compared to net forces (accelerations) are broken

    ! this is the radius of the smallest particle that is allowed to break
    breakmin=0.075

    ! teindex = breakability (tensional stress - tensional strength)
    ! The grain with the largest teindex (if positive) will break.
    if (icrush == 1) then
        do i=nbound(2)+1,n
            if (radius(i) < breakmin) cycle
            
            ff=(fx(i)*fx(i)+fy(i)*fy(i))
                
            if (ff<(0.01*stot(i)*stot(i))) then
            teindx(i)= sigmax(i)-strgth(i)
                if (teindx(i) > maxstr)then
                    maxstr=teindx(i)
                    ibpart=i
                endif
            else
                teindx(i)=100.
            endif
        enddo
    endif

    strlimit=0.5e-2
    if (icrush == 2) then
        do i=nbound(2)+1,n
            
            if (radius(i) < breakmin) cycle
            
            ff=(fx(i)*fx(i)+fy(i)*fy(i))
            
            if (ff<(0.01*stot(i)*stot(i))) then
                teindx(i)= fmax(i)-strlimit
                if (teindx(i) > maxstr)then
                    maxstr=teindx(i)
                    ibpart=i
                endif
            else
                teindx(i)=100.
            endif
        end do
    endif

    if (ibpart > 0) then
        !  Step 2. Break the grain
        ibreak=ibpart
        tresume=tau+20.
        nevent=nevent+1

        inext=inext+1


        masskept=0.92
        grainarea=grainarea-(1.-masskept)*pi*(radius(ibreak)**2)
        call breakgrain (ibreak,rnew,masskept)    

        do i=n-(pieces-2),n
            teindx(i)=teindx(ibreak)
        enddo 
        ! Change time step dt
        if (radius(n) < rsmall) then
            dt=dt*((radius(n)/rsmall)**1.5)
            rsmall=radius(n)
            write(*,199) 'rsmall = ',rsmall,'  dt = ',dt
        endif
        

        do i=1,n
            bvx(i)=0.
            bvy(i)=0.
        enddo
        !
        !

        ! Step 3. Relax the stresses on the new grains (Fake time steps)
        !             call movean(pieces,newgrn)
        call links
        
        distint=1.5
        call getneighbors(distint)
        !     make smaller neighbor list for new grains
        ibknt=0
        do k=1,contactknt
            ig=ibreak
            if (contacti(k) == ig .or. contactj(k) == ig) then
                ibknt=ibknt+1
                bcontact(ibknt)=k               
            else 
                do ig=n-pieces+2,n
                    if (contacti(k) == ig .or. contactj(k) == ig) then
                        ibknt=ibknt+1
                        bcontact(ibknt)=k
                        exit            
                    endif 
                enddo
            endif 
        enddo
        
        call breakforce(bcontact,ibknt) 

        call movebn

        tol=1.0e-7

        !       fake time step loop
  outer:do itf=1,10000

            call movean()
            call breakforce(bcontact,ibknt)
            call movebn()

            do j=1,pieces
                k=newgrn(j)
                f2(j)=fx(k)*fx(k)+fy(k)*fy(k)
                if (f2(j) > 0.01*tol*tol) cycle outer
            enddo

            exit         

        enddo outer

        print*, itf-1,' fake time steps'

        do j=1,pieces
            k=newgrn(j)
            vx(k)=0.
            vy(k)=0.
            w(k)=0.
        enddo 


    endif

       
!***************************************************************
!                   END OF BREAKING GRAINS                     *
!***************************************************************
    
end subroutine


!BREAKGRAIN***************************************************************
!                                                                        *
!This program breaks a chosen particle in a prefixed number of pieces    *
!                                                                        *
!Variables used:                                                         *
!pieces : Number of pieces in which a particle breaks                    *
!m,i : Indices                                                           *
!rnew : Radius of the broken particles                                   *
!x,y,r : Coordinates, and radius the particles                           *
!MAXN : Maximum number of particles                                      *
!x1,y1,r1,x2,y2,r2,x3,y3,r3                                              *
!x4,y4,r4,x5,y5,r5,x6,y6,r6 : Dummy variables                            *
!newgrn : Array of index of newly-formed particles                       *      *                                                                        *
!Input : x,y,r (before breaking),pieces,n,m,MAXN                  *
!                                                                        *
!Output : x,y,r (after breaking),rnew,MAXN                        *
!                                                                        *
!                                                                        *
!Programmer : Richard A. Lang                                            *
!                                                                        *
!************************************************************************* 
!=======================================================================
subroutine BREAKGRAIN(ibreak,rnew,masskept)

    integer i,j,k
    integer ibreak

    real(8) rnew
    real(8) masskept,bratio,rnew2
    real(8) dl,dtheta,theta0,xl,yl

    newgrn(1:13)=0

    newgrn(1)=ibreak
    !
    !
    ! Break the particle in pieces, calculate new radii, and create
    ! an array with the newly-formed particles       
    print*

    print*,'BREAK particle ',ibreak 
    !  masskept=fraction of original particle mass that goes into the 
    !     new pieces (rest is lost)

    theta0 = 0.5*pi+cron(ibreak)
      

    if (pieces == 3) then
        do k=1,2
            dl=0.5*rnew/sqrt(3.)
            dtheta=theta0+pi/6.*dble(k-1)
            call addgrain2(xl,yl,rnew,0)
            newgrn(k+1)=n
        enddo
        
    else if (pieces == 9) then 
        bratio=0.5
        rnew=radius(ibreak)*sqrt(masskept/(1+8*bratio*bratio))
        rnew2=rnew*bratio
        !       print*, 'Old radius=',radius(ibreak),'NEW RADII = ',rnew,rnew2
        radius(ibreak)=rnew
        radinv2(ibreak)=.125/(radius(ibreak)**3)
        call strength(ibreak)

        do k=1,8
            dtheta=pi*dble(k-1)/4.
            dl=rnew+rnew2
            xl=rx(ibreak)+dl*cos(dtheta)
            yl=ry(ibreak)+dl*sin(dtheta)
            call addgrain2(xl,yl,rnew2,0)
            newgrn(k+1)=n
        enddo

    else if (pieces >= 7) then 
        rnew=radius(ibreak)*sqrt(masskept/real(pieces))
        radius(ibreak)=rnew
        radinv2(ibreak)=.125/(radius(ibreak)**3)
        call strength(ibreak)
        do k=1,6
            dtheta=theta0+pi*dble(k-1)/3.
            dl=2.*rnew
            xl=rx(ibreak)+dl*cos(dtheta)
            yl=ry(ibreak)+dl*sin(dtheta)
            call addgrain2(xl,yl,rnew,0)
            newgrn(k+1)=n
        enddo

        if (pieces == 13) then
            do k=1,6
                dtheta=theta0+(pi/6.)+pi/3.*dble(k-1)
                dl=2.*sqrt(3.)*rnew
                xl=rx(ibreak)+dl*cos(dtheta)
                yl=ry(ibreak)+dl*sin(dtheta)
                call addgrain2(xl,yl,rnew,0)
                newgrn(k+1)=n
            enddo
        endif

    endif 
    
    vx(newgrn(1))=0.
    vy(newgrn(1))=0.
    w(newgrn(1))=0.
    
    do j = 1,pieces
        i=newgrn(j)
        call strength(i)
    enddo

end subroutine

!=======================================================================
subroutine MOVEAN()

!    *******************************************************************
!    ** FIRST PART OF VELOCITY VERLET ALGORITHM                       **
!    **                                                               **
!    ** USAGE:                                                        **
!    **                                                               **
!    ** THE FIRST PART OF THE ALGORITHM IS A TAYLOR SERIES WHICH      **
!    ** ADVANCES POSITIONS FROM T TO T + DT AND VELOCITIES FROM       **
!    ** T TO T + DT/2.  AFTER THIS, THE FORCE ROUTINE IS CALLED.      **
!    *******************************************************************

    integer i
    integer k
    
    real(8)    dt2,dtsq2
    real(8) xl0,xr0,yb0,yt0
    real(8) rms,xl

    dt2   = 0.5*dt  
    dtsq2 = dt * dt2

    xl0 = 0.
    xr0 = 0.
    yb0 = 0.
    yt0 = 0.

    do k=1,pieces 
        i=newgrn(k) 
        ! interior grains
        rms=radinv2(i)
        rx(i) = rx(i) + dt * bvx(i) + dtsq2 * fx(i) *rms
        ry(i) = ry(i) + dt * bvy(i) + dtsq2 * fy(i) *rms
        bvx(i) = bvx(i) + dt2 * fx(i) *rms
        bvy(i) = bvy(i) + dt2 * fy(i) *rms
    end do

    !  ugly bit to take care of grain wraparound
    !  in case of periodic boundaries
    if (ib(3) == -1) then
        xl=xright-xleft

        do k=1,pieces
            i=newgrn(k) 
            !            wrap atoms that have moved out of the box           
            if (rx(i) >= xright) then 
                rx(i)=rx(i)-xl
            else if (rx(i) < xleft) then
                rx(i)=rx(i)+xl
            endif
        enddo               
    endif
        
end subroutine

!=======================================================================
subroutine MOVEBN   

!    *******************************************************************
!    ** SECOND PART OF VELOCITY VERLET ALGORITHM                      **
!    **                                                               **
!    ** USAGE:                                                        **
!    **                                                               **
!    ** THE SECOND PART OF THE ALGORITHM ADVANCES VELOCITIES FROM     **
!    ** T + DT/2 TO T + DT. THIS ASSUMES THAT FORCES HAVE BEEN        **
!    ** COMPUTED IN THE FORCE ROUTINE AND STORED IN FX, FY, FZ.       **
!    *******************************************************************

    integer k,i
    real(8) DT2
    real(8) rms

    dt2 = 0.5*dt  

    do k=1,pieces
        i=newgrn(k)

        rms=radinv2(i)

        bvx(i) = bvx(i) + dt2 * fx(i) *rms
        bvy(i) = bvy(i) + dt2 * fy(i) *rms
    enddo
        

end subroutine

!======================================================================= 
subroutine breakforce (bcontact,ibknt)
     
        

!    *******************************************************************
!    ** ROUTINE TO COMPUTE FORCES AND POTENTIAL USING A LINK LIST     **
!    **                                                               **
!    ** PRINCIPAL VARIABLES:                                          **
!    **                                                               **
!    ** INTEGER NGRAINS            NUMBER OF GRAINS                   **
!    ** INTEGER NGLIST(NGRAINS)    POINTER TO 1ST ATOM IN EACH GRAIN  **    
!    ** INTEGER MX,MY              NUMBER OF CELLS IN EACH DIRECTION  **
!    ** INTEGER NCELL              NUMBER OF SMALL CELLS (M**3)       **
!    ** INTEGER LIST(N)            THE LINKED LIST                    **
!    ** INTEGER HEAD(NCELL)        THE HEAD OF CHAIN ARRAY            **
!    ** INTEGER MAP(MAPSIZ)        LIST OF NEIGHBOURING CELLS         **
!    ** real(8)    RX(N),RY(N)        POSITIONS                          **
!    ** real(8)    FX(N),FY(N)        FORCES                             **
!    ** real(8)    SIGMA              THE LJ LENGTH PARAMETER            **
!    ** real(8)    RCUT               THE CUT-OFF DISTANCE               **
!    **                                                               **
!    **                                                               **
!    ** PROCEDURE:                              **
!    **   LOOP OVER THE CELLS                            **
!    **     WITHIN THIS CELL                              **
!    **         DO INTRAMOLECULAR FORCES IN GRAINS              **
!    **         LOOP THROUGH THE LINK LIST FOR GRAINS              **
!    **        IN NEIGHBORING CELLS                      **
!    **         LOOP THROUGH THE LINK LIST FOR GRAINS              **
!    *******************************************************************
       
    integer i,j,k 
    integer ncell,ixcell,jxcell
    integer ipf
    integer bcontact(325),ibknt,ik
    
    real(8) xl2,xl

    !    ** ZERO FORCES AND POTENTIAL **
    !       print*, 'IN FORCE'
    do k=1,pieces
        i=newgrn(k)
        fx(i)=0.
        fy(i)=0.
    enddo


    ncell=mx*my
    xl2=(xright-xleft)/2.
    xl=xright-xleft
    dxc=xl/dble(mx)

    !  loop over each particle
    !   loop over neighbor list for that particle
    do ik=1,ibknt
        k=bcontact(ik)
        i=contacti(k)
        j=contactj(k)
        ipf=0    
        ixcell=1+int((rx(i)-xleft)/dxc)
        if (ixcell == 1) then 
            jxcell=1+int((rx(j)-xleft)/dxc)
            if (jxcell == mx) ipf=1
        else if (ixcell == mx) then
            jxcell=1+int((rx(j)-xleft)/dxc)
            if (jxcell == 1) ipf=1
        endif           

        call binteractions(i,j,k,ipf)
    enddo

end subroutine

!=======================================================================
subroutine binteractions(i,j,knt,ipf)
!**************************************************************
!   calculate if two particles are interacting, and the
!   interaction forces
!c
!   fnij - normal force
!   ftij - tangential force
!
!   forces on paritcle imn are postive if outward or counterclockwise
!   relative velocities of paritcle imx wrt. imn are 
!     postive if outward or counterclockwise

    integer i,j,imx,imn
    integer ipf
    integer knt
    
    real(8) rij,rijsq,sk1,comp
    real(8) fyij,fxij,fnij,sk
    real(8) xl,cohesion,ryij,rxij,rmaxsq,eqlength
    real(8) tang
    real(8) vxij,vyij,vnij,g1,g2,gammax

    tang = 180./pi

    cohesion = 0.0
    sk = 1.
    xl = xright-xleft

    ! denote the two interacting grains as imn, imx
    imn = min(i,j)
    imx = max(i,j)

    eqlength = radius(imn)+radius(imx)
    rmaxsq = eqlength*eqlength

    rxij = rx(imx) - rx(imn)
    if (ipf == 1) then
        rxij=rxij-sign(xl,rxij)
    endif
    
    ryij  = ry(imx) - ry(imn)
    rijsq = (rxij*rxij + ryij*ryij)

    ! particles are interacting
    if (rijsq < rmaxsq) then
        coordn(imn)=coordn(imn)+1
        coordn(imx)=coordn(imx)+1
        rij=sqrt(rijsq)
        comp=eqlength-rij

        sk1=sk/eqlength

        !      get relative translational velocities
        vxij=bvx(imx)-bvx(imn)
        vyij=bvy(imx)-bvy(imn)

        !      resolve relative velocties into normal and tangential parts
        vnij = (vxij*rxij + vyij*ryij)/rij

        !      calculate the maximum allowable damping factor, gamma (size dependent)
        g1=radius(imn)+radius(imx)
        g2=radius(imn)*radius(imn)+radius(imx)*radius(imx)
        gammax=4*radius(imn)*radius(imx)/sqrt(g1+g2)

        !      calculate normal and tangential forces
        fnij = -sk1*comp + 0.25*gamma*gammax*vnij


        !      resolve forces back into x and y and add to sums for atom imn
        fxij= (rxij*fnij-ryij*ftij)/rij
        fyij= (ryij*fnij+rxij*ftij)/rij
        fx(imn)   = fx(imn) + fxij
        fy(imn)   = fy(imn) + fyij

        !   resolve opposite forces back into x and y and add to sums for imx
        fx(imx)   = fx(imx) - fxij
        fy(imx)   = fy(imx) - fyij

        contfn(knt)=fnij
        contft(knt)=0.

    else

        ! particles are not interacting
        contfn(knt)=0.
        contft(knt)=0.
    endif

end subroutine

 end module
