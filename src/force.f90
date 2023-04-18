module mod_force
    use mycommons 
    use mod_principalstress, only : STRESS
    use mod_inputReader
    implicit none

    real(8) fsin

    private
    
    public :: force, fsin

    contains
    
!===================================================================
subroutine force (gx,gy,ev,presx,presy,presy0,fytop,fybot, &
           fxtop,fxbot,ipfc,crush1,ifluid,fgrain,fbc)
     

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
!    ** PROCEDURE:                                          **
!    **   LOOP OVER THE CELLS                                      **
!    **     WITHIN THIS CELL                                      **
!    **         DO INTRAMOLECULAR FORCES IN GRAINS                  **
!    **         LOOP THROUGH THE LINK LIST FOR GRAINS                  **
!    **          IN NEIGHBORING CELLS                              **
!    **         LOOP THROUGH THE LINK LIST FOR GRAINS                  **
!    **                                                               **
!    **                                                               **
!    **  16/09/09 - Liran Goren.                                      **
!    **  Modified from force15h.f contains changes with the way the   **
!    **  the fluid force and the confining force are calculated on    **
!    **  the walls.                                                   **   
!    *******************************************************************
      
    integer i,j,k 
    integer ipf
    integer ipfc
    integer ixcell,jxcell
    integer crush1
    integer ifluid
    integer fbc

    real(8) gx,gy,ev,fytop,fybot
    real(8) fxbot,fxtop,presx,presy,presy0
    real(8) xl
    real(8) fgrain

    !    ** ZERO FORCES AND POTENTIAL **

    ! gravity forces proportional to mass of particle, 

    !print*, 'XX f 1-: fx,fy = ',fx(nbound(1)+2),fy(nbound(1)+2)
    call initf(presx,presy,presy0,gx,gy,fgrain)
    !print*, 'XX f 1+: fx,fy = ',fx(nbound(1)+2),fy(nbound(1)+2)

    ev = 0.0

    coordn(1:n) =0

    if (crush1 == 1) then 
      u1(1:n)=0.
      u2(1:n)=0.
      u3(1:n)=0.
    endif

    xl=xright_out-xleft_out
    dxc=xl/dble(mx)

    if (ipfc == 1) then
      ncontacts=0
      sumcomp=0d0
      maxcomp=0d0
      sumvn=0d0
      maxvn=0d0
      refi=0
      refj=0
    endif

    !  loop over each particle
    !   loop over neighbor list for that particle
    do k=1,contactknt
      i=contacti(k)
      j=contactj(k)
      ipf=0
      if (mx == 2) then 
          ipf=2
      else
          ixcell=1+int((rx(i)-xleft_out)/dxc)
          if (ixcell == 1) then 
              jxcell=1+int((rx(j)-xleft_out)/dxc)
              if (jxcell == mx) ipf=1
          else if (ixcell == mx) then
              jxcell=1+int((rx(j)-xleft_out)/dxc)
              if (jxcell == 1) ipf=1
          endif           
      endif
        
      if(intruder==1.)then
          if((i==n).or.(j==n)) ipf=2
      endif
     
      !print*, 'XX f 2-: fx,fy = ',fx(nbound(1)+2),fy(nbound(1)+2)
      call interactions(i,j,k,ev,ipf,ipfc,crush1)
      !print*, 'XX f 2+: fx,fy = ',fx(nbound(1)+2),fy(nbound(1)+2)
    enddo
    !     add the force due to the fluid gradient (liran)
    !print*, 'XX f 3-: fx,fy = ',fx(nbound(1)+2),fy(nbound(1)+2)
    if (ifluid == 0) call fluidf(fbc)
    !print*, 'XX f 3+: fx,fy = ',fx(nbound(1)+2),fy(nbound(1)+2)

    !print*, 'XX f 4-: fx,fy = ',fx(nbound(1)+2),fy(nbound(1)+2)
    call bsum(fytop,fybot,fxtop,fxbot)
    !print*, 'XX f 4+: fx,fy = ',fx(nbound(1)+2),fy(nbound(1)+2)
          
end subroutine     
                       
 
!===================================================================       
subroutine initf(presx, presy,presy0,gx,gy,fgrain)
!*****************************************************
! set up initial forces on particles
!  * body forces on internal particles
!  * applied forces on boundary particles  
      
    integer i

    real(8) presy0,gx,gy,presx,presy
    real(8) fgrain,fgrain0,fgrainx,fgrain0x
    real(8) yb1,xb1

    fgrain = 0.
    fgrainx = 0d0
    fgrain0 = 0.
    fgrain0x = 0d0
    ! presy, presy0 are applied stresses
    ! convert to force on each grain

    do concurrent (i = 1:nbound(1))
     fgrain0 = fgrain0 + presy0*pi*radius(i)**2
     fgrain0x = fgrain0x + presx*pi*radius(i)**2
    enddo

    fgrain0 = fgrain0/dble(nbound(1))
    fgrain0x = fgrain0x/dble(nbound(1))

    do concurrent (i = nbound(1)+1:nbound(2))
     fgrain = fgrain + presy*pi*radius(i)**2
     fgrainx = fgrainx + presx*pi*radius(i)**2
    enddo

    fgrain = fgrain/dble(nbound(2) - nbound(1))
    fgrainx = fgrainx/dble(nbound(2) - nbound(1))
       
    ! bottom boundary
    do concurrent (i = 1:nbound(1))
       fx(i) = 0.0
       if (ib(1)==1)                   fx(i)=fgrain0x
       if (ib(1) == 6 .or. ib(1) == 3) fx(i)=fb(1)
       if (ib(1) == 7)                 fx(i)=fsin
       
       if (ib(1)==1 .or. ib(1)==2) then
          fy(i)=fgrain0 
       else
          fy(i)=0.
       endif
    enddo

    ! top boundary
    do concurrent (i = nbound(1)+1:nbound(2))
      !fx(i)=0.0
      if (ib(2)==1) then
          fx(i)=fgrainx
      elseif (ib(2) == 5) then
          fspring = kspring*(dspring-dwall)
          fx(i) = fspring
      else
          fx(i)=0.
      endif
      
      if (any(ib(2)==[1,2,5])) then
          !fy(i) = -fgrain-wallv*vy(nbound(2))
          fy(i) = -fgrain
      else
          fy(i)=0.
      endif
    enddo

    ! internal particles
    do i = nbound(4)+1,n
        tq(i)=0.

        if((intruder==1).and.(i==n))then
            fx(i) = gx/radinv2(i)
            fy(i) = gy/radinv2(i)
            xb1=gx**2+gy**2
            
            if(gy/=0)then
                if(ry(i)+radius(i)<wl)then  
                    !fully immersed grain      
                    fy(i) = fy(i)-gy*fstress*6*radius(i)**2
                elseif(ry(i)-radius(i)>wl)then
                    continue
                elseif(ry(i)>wl) then !partially immersed grain
                    yb1=2d0*atan(sqrt(radius(i)**2-(ry(i)-wl)**2) / (ry(i)-wl))
                    fy(i) = fy(i) - gy*3d0/pi*fstress*radius(i)**2*(yb1-sin(yb1))                 
                else
                    yb1=2d0*atan(sqrt(radius(i)**2-(ry(i)-wl)**2) / (wl-ry(i)))
                    fy(i) = fy(i)-gy*fstress*6*radius(i)**2+ &
                           gy*3d0/pi*fstress*radius(i)**2*(yb1-sin(yb1))    
                endif
            endif           
        else
            fx(i) = gx*8*radius(i)**3
            fy(i) = gy*8*radius(i)**3
            xb1=gx**2+gy**2
            if(gy/=0)then
                if(ry(i)+radius(i)<wl)then  
                    !fully immersed grain      
                    fy(i) = fy(i)-gy*fstress*8*radius(i)**3
                 elseif(ry(i)-radius(i)>wl)then
                    continue
                elseif(ry(i)>wl)then !partially immersed grain
                    yb1=radius(i)-ry(i)+wl
                    fy(i) = fy(i) - gy*6d0*fstress*yb1**2*(radius(i)-yb1/3d0)       
              
                else
                    yb1=radius(i)-wl+ry(i)
                     fy(i) = fy(i)- gy*fstress*8*radius(i)**3 + &
                         gy*6d0*fstress*yb1**2*(radius(i)-yb1/3d0)    
                endif
            endif 
         endif
    enddo

    !external force to 
    if(use_ext_force) then
        if(ext_force_to <= nbound(4) .or. ext_force_to > n) then
            print*, "ERROR ext_force_to is part of bounds or > number of grains"
            stop
        else
            fx(ext_force_to) = fx(ext_force_to) + ext_force(1)
            fy(ext_force_to) = fy(ext_force_to) + ext_force(2)
        endif
    endif
      
end subroutine
      

!=======================================================================
subroutine bsum(fytop,fybot,fxtop,fxbot)
!*************************************************************        
! get the sums of stresses on the rigid boundaries
    integer i

    real(8) fytop,fybot,fxtop,fxbot

    ! find average normal and shear forces on boundaries
    fytop=0.
    fybot=0.
    fxtop=0.
    fxbot=0.

    do concurrent (i = 1:nbound(1))
       fybot=fybot + fy(i)
       fxbot=fxbot + fx(i)
    enddo

    do concurrent (i = nbound(1)+1:nbound(2))
       fytop=fytop + fy(i)
       fxtop=fxtop + fx(i)
    enddo

    fybot=fybot/dble(nbound(1))
    fxbot=fxbot/dble(nbound(1))
    fytop=fytop/dble(nbound(2)-nbound(1))
    fxtop=fxtop/dble(nbound(2)-nbound(1))

    ! to keep rigid, apply the average force to all boundary particles
    do concurrent (i = 1:nbound(1))
       fy(i)=fybot
       fx(i)=fxbot
    enddo

    do concurrent(i = nbound(1)+1:nbound(2))
       fy(i)=fytop
       fx(i)=fxtop
    enddo
    !      print*, 'new force on bottom particle=',fy(nbound(1))
    !      print*, 'new force on top particle=',fy(nbound(1)+1)

    
end subroutine


!=======================================================================
subroutine interactions(i,j,knt,ev,ipf,ipfc,crush1)
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
    integer ipfc
    integer knt
    integer iycelln, iycellx
    integer crush1
    
    real(8) vxij,vyij,vnij,vthij,sk2,rij,rijsq,sk1,comp,vtij0
    real(8) fyij,fxij,fnij,vtij,ftmax,sk
    real(8) xl,cohesion,ev,ryij,rxij,rmaxsq,eqlength
    real(8) slipij
    real(8) tang
    real(8) gammax,g1,g2
    real(8) xm,ym,xi,yi
    real(8) rxij2,rx1,rx2
    real(8) rbar, ar
    real(8) nu
    real(8) histdy

    tang=180./pi

    cohesion = 0.0
    sk = 1.
    xl = xright_out-xleft_out

    ! denote the two interacting grains as imn, imx
    imn = min(i,j)
    imx = max(i,j)

    eqlength = radius(imn)+radius(imx)
    rmaxsq = eqlength*eqlength

    ryij = ry(imx) - ry(imn)
    rxij = rx(imx) - rx(imn)
    
    if (ipf == 1) then
        rxij  =rxij-sign(xl,rxij)
    else if (ipf == 2) then
        rxij2 = rxij-sign(xl,rxij)
        if (abs(rxij2) < abs(rxij)) rxij=rxij2
    endif
    
    rijsq = (rxij*rxij + ryij*ryij)

    ! particles are interacting
    if (rijsq < rmaxsq) then
        coordn(i)=coordn(i)+1
        coordn(j)=coordn(j)+1

        rij=sqrt(rijsq)
        comp=eqlength-rij

        ! stiffness for two different grains, length dependent

        !  LINEAR SPRING CONSTANTS
        !  stiffness for E=1 for both grains, but length dependent
        sk1=skmol
        sk2=skshear

        !     potential energy from overlap
        ev = ev + 0.5 * skmol*comp*comp

        !      get relative translational velocities
        vxij=vx(imx)-vx(imn)
        vyij=vy(imx)-vy(imn)

        !      get relative rotations
        vthij=w(imx)*radius(imx)+w(imn)*radius(imn)

        !      resolve relative velocties into normal and tangential parts
        vnij = (vxij*rxij + vyij*ryij)/rij
        vtij0 = (-vxij*ryij + vyij*rxij)/rij
        vtij = vtij0 - vthij
        slipij=contft(knt)

        !      calculate the maximum allowable damping factor, gamma (size dependent)
        g1=radius(imn)+radius(imx)
        g2=radius(imn)*radius(imn)+radius(imx)*radius(imx)
        gammax=4*radius(imn)*radius(imx)/sqrt(g1+g2)

        if(contact_model==0)then
            !      calculate normal and tangential forces
            fnij = -sk1*comp + gamma*vnij
            slipij = slipij + sk2*dt*vtij
        elseif(contact_model==1)then
            !  HERTZ-MINDLIN CONTACT MODEL
            !     harmonic mean of grain radii, rbar
            !     contact radius,ar
            rbar=2*radius(imn)*radius(imx)/(radius(imn)+radius(imx))
            ar=sqrt(rbar*comp)

            ! calculate Poisson's ratio from skshear
            nu = (6.-2.*skshear)/(6.-skshear)

            !      calculate normal and tangential forces
            fnij = -(sqrt(2.)/3./(1.-nu*nu))*ar*comp + gamma*ar*vnij


            ! following Potopov, keep track of tangential force on contact

            !      update tangential force on contact
            !                 call getslip(imn,imx,slipij,jptr)
            !               print*, imn,imx,slipij,jptr

            ! FOR MINDLIN TANGENTIAL CONTACT
            ! sk2=ratio of tangential stiffness to normal stiffness
            !         =6*(1-nu)/(2-nu)  , where nu=Poisson's ratio
            !                  slipij = slipij + skshear*ar*dt*vtij
            slipij = slipij +  &
            (2.*sqrt(2.)/(2-nu)/(1+nu))*ar*dt*vtij
        elseif(contact_model==2)then
            ! HERTZ-MINDLIN CONTACT MODEL according to Di Renzo et al. Chemical Engineering Science 59, 2004
            rbar=2*radius(imn)*radius(imx)/(radius(imn)+radius(imx))
            ar=sqrt(rbar*comp)
            nu = (6.-2.*skshear)/(6.-skshear)
            fnij = -(sqrt(2.)/3./(1.-nu*nu))*ar*comp + gamma*ar*vnij
            slipij = slipij + sqrt(2.)/(2-nu)/(1+nu)*ar*dt*vtij
        endif !contact_model
        
        ! To increase the friction at the top of the pack
        !                  if (i >= int(n/2) .or. j >= int(n/2)) then
        !                     ftmax = 3.*friction*fnij
        !                  else
        ftmax = friction*fnij
        !                  endif

        if (abs(slipij) > abs(ftmax)) then
            slipij=sign(ftmax,slipij)
        endif
        
        ftij=slipij

        contfn(knt)=fnij
        contft(knt)=slipij

        ! specials: Calculate maximum and mean overlap and normal velocity
        if (ipfc == 1) then
            ncontacts=ncontacts+1
            sumcomp=sumcomp+comp
            if (comp > maxcomp) then
                maxcomp=comp
                refi=imn
                refj=imx
            endif
            sumvn=sumvn+abs(vnij)
            maxvn=max(maxvn,abs(vnij))

            contang(knt)=tang*datan(ryij/rxij)
        endif

        !      resolve forces back into x and y and add to sums for atom imn
        fxij = (rxij*fnij-ryij*ftij)/rij
        fyij = (ryij*fnij+rxij*ftij)/rij
        
        fx(imn) = fx(imn) + fxij
        fy(imn) = fy(imn) + fyij
        
        tq(imn) = tq(imn) + ftij

        !   resolve opposite forces back into x and y and add to sums for imx
        fx(imx)   = fx(imx) - fxij
        fy(imx)   = fy(imx) - fyij
        tq(imx)   = tq(imx) + ftij

        histdy=(ytop-ybot)/my 
                       
        iycelln = 1+int((ry(imn)-ybot)/histdy)
        iycellx = 1+int((ry(imx)-ybot)/histdy)
        if(iycelln > my) iycelln = my
        if(iycellx > my) iycellx = my
        if((iycelln < 1) .or. (iycellx < 1)) print*, 'WARNING: iycelln, iycellx',iycelln,iycellx
        
        sxxt(iycelln) = sxxt(iycelln) + fxij*rxij/2d0
        sxyt(iycelln) = sxyt(iycelln) + fxij*ryij/2d0  
        syxt(iycelln) = syxt(iycelln) + fyij*rxij/2d0
        syyt(iycelln) = syyt(iycelln) + fyij*ryij/2d0
        sxxt(iycellx) = sxxt(iycellx) + fxij*rxij/2d0
        sxyt(iycellx) = sxyt(iycellx) + fxij*ryij/2d0  
        syxt(iycellx) = syxt(iycellx) + fyij*rxij/2d0
        syyt(iycellx) = syyt(iycellx) + fyij*ryij/2d0

        if (crush1 == 1) then 
            xm= fxij
            ym= fyij
            xi= -xm
            yi= -ym
            if (ipf == 1) then 
                if (rx(imn) < rx(imx)) then
                    rx1=rx(imn)+xl
                    rx2=rx(imx)
                else
                    rx1=rx(imn)
                    rx2=rx(imx)+xl
                endif
            else
                rx1=rx(imn)
                rx2=rx(imx)
            endif

            call stress(imn,imx,rx1,rx2,xm,ym,xi,yi)
        endif

    else

        ! particles are not interacting
        contfn(knt)=0.
        contft(knt)=0.
        if (ipfc == 1) then
            contang(knt)=0.
        endif
    endif
    
end subroutine

! the fluid pressure force is gradP*m/(rho * (1-phi)) or in a different way: gradP*V/N where V is a cell volume and N in 
! the number of particle in the cell given in the array numpar
! after using non-dimensionalization we get that the fluid force is grad'P'*m'/(1-phi) * (3*pi/32)

!=======================================================================
subroutine fluidf(fbc)
      
    integer ix,ixvxh,iy,iyvyh,k,MyCount
    integer fbc

    real(8) yb1,xb1,dxr,dxrvx,dyr,dyrvy,xl,addx,addy
    real(8) prefac
    
    logical :: stop_yes = .false.

    prefac = pi/6.
    dxc=(xright-xleft)/mx
    dyc=(ytop-ybot)/my
    !vol =dxc*dyc
    ! inner grains
    !do concurrent (k = nbound(4) + 1: n)
    do concurrent (k = 1: n)
      if (gtype(k) /= 0) then
          if (rx(k) < xleft .or. rx(k) > xright &
               .or. ry(k) < ybot .or. ry(k) > ytop+ftiny) then
              print*, k, rx(k),ry(k),xleft,xright,ybot,ytop
              stop_yes = .true.
          endif         
        
          if((intruder==1).and.(k==n))then
              xl=xright-xleft   
              addx=0
              addy=0
              MyCount=0
              
                do ix=1,mx
                  do iy=1,my+1
                      xb1=xleft+(ix-1)*dxc
                      yb1=ybot+(iy-1)*dyc
                      dxr=rx(n)-xb1
                      dyr=ry(n)-yb1
                      if (ib(3) < 0) then
                        if(abs(dxr)>(0.5d0*xl)) dxr=dxr-sign(xl,dxr)
                      endif
                      if((dxr**2+dyr**2)<rintruder**2)then
                          MyCount=MyCount+1
                          addx=addx+gradpx(ix,iy)
                          addy=addy+gradpy(ix,iy)
                      endif
                  enddo
              enddo

              addx=addx/MyCount
              addy=addy/MyCount
              addx = addx*pi*radius(k)**2
              addy = addy*pi*radius(k)**2
          else
              ix = int((rx(k) - xleft)/dxc) + 1
              if (ix > mx) ix = mx
              if (rx(k) < (xleft + 0.5d0*dxc)) then
                ixvxh = mx
              else
                ixvxh = int((rx(k)-(xleft + 0.5d0*dxc))/dxc)+1
              endif
              iy = int((ry(k) - ybot)/dyc) + 1
              if (iy > my) iy = my
              if (ry(k) < (ybot + 0.5d0*dyc)) then
                iyvyh = 1
              else
                iyvyh=int((ry(k)-(ybot + 0.5d0*dyc))/dyc)+2
              endif

              !xb1=xleft+(ix-1)*dxc
              !yb1=ybot+(iy-1)*dyc
              !dxr=(rx(k)-xb1)/dxc
              !dyr=(ry(k)-yb1)/dyc
              dxr=(rx(k)-xleft)/dxc - (ix-1)
              dyr=(ry(k)-ybot)/dyc - (iy-1)
              if (rx(k) < xleft + 0.5d0*dxc) then
                 dxrvx = (rx(k) - xleft)/dxc + 0.5d0
              else
                 dxrvx = (rx(k) - xleft)/dxc - (ixvxh - 0.5d0)
              endif
              dyrvy = (ry(k) - ybot)/dyc - (iyvyh - 1.5d0)

!              addx = (1.-dyr)*(1.-dxr)*gradpx(ix,iy)/solfra(ix,iy) + &
!                  dyr*(1.-dxr)*gradpx(ix,iy+1)/solfra(ix,iy+1) + &
!                  dyr*dxr*gradpx(ix+1,iy+1)/solfra(ix+1,iy+1) + &
!                  (1.-dyr)*dxr*gradpx(ix+1,iy)/solfra(ix+1,iy)
!              addy = (1.-dyr)*(1-dxr)*gradpy(ix,iy)/solfra(ix,iy) + &
!                  dyr*(1.-dxr)*gradpy(ix,iy+1)/solfra(ix,iy+1) + &
!                  dyr*dxr*gradpy(ix+1,iy+1)/solfra(ix+1,iy+1) + &
!                  (1.-dyr)*dxr*gradpy(ix+1,iy)/solfra(ix+1,iy)
              addx=(1.-dyr)*(1.-dxrvx)*gradpx(ixvxh,iy)/solfra(ixvxh+1,iy)+ &
                  dyr*(1.-dxrvx)*gradpx(ixvxh,iy+1)/solfra(ixvxh+1,iy+1)+ &
                  (1.-dyr)*dxrvx*gradpx(ixvxh+1,iy)/solfra(ixvxh+1,iy)+ &
                  dyr*dxrvx*gradpx(ixvxh+1,iy+1)/solfra(ixvxh+1,iy+1)
              addy=(1.-dyrvy)*(1.-dxr)*gradpy(ix,iyvyh)/solfra(ix,iyvyh)+ &
                  dyrvy*(1.-dxr)*gradpy(ix,iyvyh+1)/solfra(ix,iyvyh)+ &
                  (1.-dyrvy)*dxr*gradpy(ix+1,iyvyh)/solfra(ix+1,iyvyh)+ &
                  dyrvy*dxr*gradpy(ix+1,iyvyh+1)/solfra(ix+1,iyvyh)

              addx = addx*prefac*8*radius(k)**3
              addy = addy*prefac*8*radius(k)**3

              ! do not impose P grad force where bc conditions dictate grad P = 0
              if (bdgrain(k) > 0) then
                if (fbc == 2 .or. fbc == 6 .or. ((fbc == 4).and.(k > nbound(1))) .or. ((fbc == 5).and.(k <= nbound(1)))) then
                  addy = 0d0  ! will be resolved below
                endif
              endif

          endif            
        
          fx(k) = fx(k) - addx                        
          fy(k) = fy(k) - addy
      endif
      !print*, 'XX addy,i = ',addy,k,dyrvy,dxr,gradpy(ix,iyvyh),gradpy(ix,iyvyh+1),gradpy(ix+1,iyvyh),gradpy(ix+1,iyvyh+1), &
      ! solfra(ix,iyvyh),solfra(ix+1,iyvyh)
    enddo       
    !print*, 'XX fluidf: fx,fy = ',fx(nbound(1)+2),fy(nbound(1)+2)
    
    if(stop_yes) stop     

!    if (fbc == 1 .or. fbc == 4) then
!    !     bottom grains
!      !if (fbc == 1) then
!          do concurrent (k = 1:nbound(1))
!            if ((rx(k) >= xleft).and.(rx(k) <= xright)) then
!              ix = int((rx(k) - xleft)/dxc) + 1
!              iy = 1
!              xb1 = xleft+(ix-1)*dxc
!              yb1 = ybot
!              dxr = (rx(k)-xb1)/dxc
!              dyr = 0.
!              addx = 1.*(1.-dxr)*gradpx(ix,iy)/solfra(ix,iy) + &
!                 dxr*gradpx(ix+1,iy)/solfra(ix+1,iy)
!    !     mass of boundary grains is half
!              addx = addx*prefac*8*radius(k)**3/2
!              fx(k) = fx(k) - addx
!
!              addy = 1.*(1-dxr)*gradpy(ix,iy)/solfra(ix,iy) + &
!                      dxr*gradpy(ix+1,iy)/solfra(ix+1,iy)
!
!              addy = addy*prefac*8*radius(k)**3/2
!              fy(k) = fy(k) - addy
!            endif
!          enddo
!    endif
!    if (fbc == 1 .or. fbc == 5) then
!    !     top grains
!      do concurrent (k = nbound(1) + 1:nbound(2))
!        if ((rx(k) >= xleft).and.(rx(k) <= xright)) then
!          ix = int((rx(k) - xleft)/dxc) + 1
!          iy = my+1
!
!          xb1=xleft+(ix-1)*dxc
!          yb1=ytop
!
!          dxr=(rx(k)-xb1)/dxc
!          dyr= 0.
!
!          addx = 1.*(1.-dxr)*gradpx(ix,iy)/solfra(ix,iy) + &
!                   dxr*gradpx(ix+1,iy)/solfra(ix+1,iy)
!    !     mass of boundary grains is half
!          addx = addx*prefac*8*radius(k)**3/2
!          fx(k) = fx(k) - addx
!          addy = 1.*(1-dxr)*gradpy(ix,iy)/solfra(ix,iy) + &
!                  dxr*gradpy(ix+1,iy)/solfra(ix+1,iy)
!
!          addy = addy*prefac*8*radius(k)**3/2
!
!          fy(k) = fy(k) - addy
!        endif
!      enddo
!    endif

    if (fbc == 2 .or. fbc == 5 .or. fbc == 6) then
    !  bottom grains
      do concurrent (k = 1:nbound(1))
        if ((rx(k) >= xleft).and.(rx(k) <= xright)) then
          ix = int((rx(k) - xleft)/dxc) + 1          
          xb1=xleft+(ix-1)*dxc            
          dxr=(rx(k)-xb1)/dxc
    !     when calculating the force as pressure* surface area
          addy = ((1-dxr)*press(ix,1) + &
                dxr*press(ix+1,1))*pi*radius(k)**2
          fy(k) = fy(k) - addy
        endif
      enddo
    endif
    if (fbc == 2 .or. fbc == 4 .or. fbc == 6) then
    ! top grains
      do concurrent (k = nbound(1)+1:nbound(2))
        if ((rx(k) >= xleft).and.(rx(k) <= xright)) then
          ix = int((rx(k) - xleft)/dxc) + 1          
          iy = my+1
          xb1=xleft+(ix-1)*dxc            
          dxr=(rx(k)-xb1)/dxc
    !     when calculating the force as pressure* surface area
          addy = ((1-dxr)*press(ix,iy) + &
                 dxr*press(ix+1,iy))*pi*radius(k)**2
    !     here we use '+' sign because positive PP sould act to push the top wall up  
          fy(k) = fy(k) + addy
        endif
      enddo
    endif
       
end subroutine

end module
