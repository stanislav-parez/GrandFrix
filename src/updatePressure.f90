module mod_updatepressure
    use mycommons 
    use mod_inputReader
    use mod_outputWriter
    implicit none

    private

    public :: updatepressure

    contains
!ccccccccccccccc 15/2/09
!   File should go with SmoothNonExtended
!   Modified from updatePressure3, but with the removal of the extra layers 
!   at the top and bottom of the visible domain.  

!cccccccccccccccc 16/9/09
!   Modified from updatePressure4.f but instead of regular mean for the permeability 
!   between grid points it uses harmonic mean.
! ***************************************************************************************
!       FROM NUMERICAL RECIPES P. 43
!=======================================================================
subroutine tridag(a,b,c,r,u,n)

    integer,intent(in) :: n
    real(8),intent(in) :: a(n), b(n), c(n), r(n)
    real(8),intent(out) :: u(n)

    integer j    
    real(8) bet,gam(n)

    if(b(1)==0.) then
        print *, 'pb pivot in tridag 1'
        stop
    endif
    bet=b(1)
    u(1)=r(1)/bet
    
    do j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if(bet==0.)then
            print*, 'stopped at ',j
            stop 'pb pivot in tridag'
        endif
        u(j)=(r(j)-a(j)*u(j-1))/bet
    enddo
    
    do j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
    enddo

end subroutine
!=======================================================================
!       FROM NUMERICAL RECIPES P. 68
subroutine cyclic(a,b,c,alpha,beta,r,x,n)

    integer,intent(in) :: n
    real(8),intent(in) :: alpha, beta, a(n), b(n), c(n), r(n)
    real(8),intent(out) :: x(n)

    integer i
    real(8) fact,gamma,bb(n),u(n),z(n)
        
    if(n <= 2) then
        print *, 'n too small in cyclic'
        stop
    endif
    gamma = -b(1)
    bb(1)=b(1)-gamma
    bb(n)=b(n)-alpha*beta/gamma
    
    do i=2,n-1
        bb(i) = b(i)
    enddo
    
    call tridag(a,bb,c,r,x,n)
    u(1)=gamma
    u(n)=alpha
    
    do i=2,n-1
        u(i)=0
    enddo
    
    call tridag(a,bb,c,u,z,n)
    fact=(x(1)+beta*x(n)/gamma)/ &
    (1.+z(1)+beta*z(n)/gamma)
    
    do i=1,n
        x(i)=x(i)-fact*z(i)
    enddo

end subroutine

!=======================================================================
! ************************************************************************************
! dont forget to initialize the pressure matrix *************************************
subroutine updatepressure(bc,numit)

! When we introduce the fluid and calculate the fluid pressure, we must introduce scaling factors.    
! Note that in the dry granular everything is really non dimensional and only when presenting the result for 
! a physcial system we need to scale everything up.
!
! The equation to solve is: beta*Phi*(dp/dt) + nabla u_s - nabla[kappa/mu nabla P] = 0
! we use: P0 = k/x0
!         V0 = x0*sqrt(k/m)
!         t0 = sqrt(m/k)
!         kappa0 = alpha*x0^2, where alpha = 1/(45*12) = 1/540, 
!                  and the permeability is calculated with a radius and not a diamter
!         x0 = d, grain diameter 
! The non dim fluid equation is then:
!  (dp/dt) + (x0/k*beta)*(1/phi)*nabla u_s - sqrt(m/k)*(alpha/beta*mu)*(1/phi) nabla[kappa nabla P] = 0
!  
!  assign values: beta = 4.5e-10 1/Pa, E (bulk modulus of grains) = 8e10 Pa, x0 = d = 1 mm = 1e-3 m 
!  therefore k = 8e10*1e-3 = 8e7 Pa m
!  m = pi*(4/3)*(d/2)^3*rho_g = 1.3823e-6 kg
!  (rho_g = 2640 kg/m^3)
! Rge resulting coefficients are: (x0/k*beta) = 2.78e-2 and sqrt(m/k)*(alpha/beta*mu) = 540.94

! acoef = 1/E*fluid compessibility
!      parameter (acoef=2.78e-2)
!      parameter (acoef=2.777777778e-2)
!      parameter (acoef=0.169635) !for d=1e-2m,E=1.31e10Pa,rho=2640kg/m^3

! bcoef is Pe^-1 = k/(mu * fluid compessibility *d *velocity)
!      parameter (bcoef= 540.94)
!      parameter (bcoef= 540.9411345)
!parameter (bcoef= 13367.782) !for d=1e-2m,E=1.31e10Pa,rho=2640kg/m^3
!      parameter (bcoef= 7218602.28) !for d=1e-2m,E=1.31e10Pa,rho=2640kg/m^3

    integer i,j,debug
    integer bc,bcl
    integer numit
    integer stepprint
    
    real(8) a(MAXD+1),b(MAXD+1),c(MAXD+1)
    real(8) r(MAXD+1),tempP(MAXD+1)
    real(8) halfP(MAXD+1,MAXD+1)
    real(8) csqx,csqy,cdiv,cadx,cady
    real(8) plhli,plhlj,mnhli,mnhlj,coef2
    real(8) alpha,beta
    real(8)  locdt
!    real(8) maxperm
!    real(8) maxdim


    debug = 0
!    maxperm = 0
!    do i = 1,mx+1
!        do j = 1,my+1
!            if (perm(i,j) > maxperm)maxperm=perm(i,j)
!        enddo
!    enddo
    
!    if (dxc > dyc) then
!        maxdim = dxc
!    else
!        maxdim = dyc
!    endif

    bcl = bc

    if(bcl == 3) bcl=1
    locdt = dt/real(TRATIO)

    stepprint = TRATIO

    ! OUTPUT ONLY IF EXTERNAL TOPRINT IS 1 AND IF STEPPRINT == CURRENTINTERATION.
    ! NORMALLY STEPPRINT = TRATIO SO ONLY ONE FILES WILL BE PRINTED FOR TRATIO ITERATIONS.
    ! FOR STEPPRINT = TRATIO/2 WE WILL HAVE TWO PRINTTED FILES FOR TRATIO ITERATION.

    !       first part of ADI update for time n+1/2, for each of the rows.

    csqx = (locdt*bcoef)/(2.*dxc*dxc)
    csqy = (locdt*bcoef)/(2.*dyc*dyc)
    cdiv = (locdt*acoef)/2.
    cadx = locdt/(2.*dxc)
    cady = locdt/(2.*dyc)

    if (yn_explicit) then
      csqx = 2d0*csqx
      csqy = 2d0*csqy
      cdiv = 2d0*cdiv
      cadx = 2d0*cadx
      cady = 2d0*cady
    endif

    ! %%%%%%%%%%%%%%%%%%%% NEUMANN BOUNDARY CONDITION WITH DERIVATIVE = 0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if (bcl == 2) then
        do i=1,mx
            do j=1,my+1
                halfP(i,j)=0.
            enddo
            tempP(i) = 0.
            a(i) = 0.
            b(i) = 0.
            c(i) = 0.
            r(i) = 0.
        enddo
        
        do j=1,my+1
            do i = 1,mx
            coef2 = phimat(i,j)
            if (j /= 1 .and. j /= my+1) then
                plhlj =2.*perm(i,j+1)*perm(i,j)/ &
                (perm(i,j+1)+perm(i,j))
                mnhlj = 2.*perm(i,j-1)*perm(i,j)/ &
                (perm(i,j-1)+perm(i,j))
                r(i)= press(i,j-1)*(csqy*mnhlj)+ &
                press(i,j)* &
                (coef2-csqy*(plhlj+mnhlj))+  &
                press(i,j+1)*(csqy*plhlj) - &
                cdiv*divvel(i,j)
            else

                !       BOTTOM BOUNDARY 
                if (j == 1) then
                    plhlj = 2.*perm(i,j+1)*perm(i,j)/ &
                    (perm(i,j+1)+perm(i,j))
                    mnhlj = 0
                    r(i) = 0. + &
                    press(i,j)* &
                    (coef2-csqy*plhlj)+  &
                    press(i,j+1)*(csqy*plhlj) - &
                    cdiv*divvel(i,j)
                    !       TOP BOUNDARY
                else
                    plhlj = 0.
                    mnhlj = 2.*perm(i,j-1)*perm(i,j)/ &
                    (perm(i,j-1)+perm(i,j))
                    r(i)= press(i,j-1)*(csqy*mnhlj)+ &
                    press(i,j)* &
                    (coef2-csqy*mnhlj)+  &
                    0. - &
                    cdiv*divvel(i,j)

                endif
            endif

            if (i /= 1 .and. i /= mx) then
                plhli = 2.*perm(i+1,j)*perm(i,j)/ &
                (perm(i+1,j)+perm(i,j))
                mnhli = 2.*perm(i-1,j)*perm(i,j)/ &
                (perm(i-1,j)+perm(i,j))
                a(i) =  - (csqx*mnhli)
                b(i) = coef2 + csqx*(plhli + mnhli)
                c(i) =  - (csqx*plhli)
            else
                !       LEFT BOUNDARY
                if (i == 1)then
                    plhli = 2.*perm(i+1,j)*perm(i,j)/ &
                    (perm(i+1,j)+perm(i,j))
                    mnhli = 2.*perm(mx,j)*perm(i,j)/ &
                    (perm(mx,j)+perm(i,j))
                    a(i) = 0.
                    beta =  - (csqx*mnhli)
                    b(i) = coef2 + csqx*(plhli + mnhli)
                    c(1) =  - (csqx*plhli)
                    !       RIGHT BOUNDARY
                else
                    plhli = 2.*perm(1,j)*perm(i,j)/ &
                    (perm(1,j)+perm(i,j))
                    mnhli = 2.*perm(i-1,j)*perm(i,j)/ &
                    (perm(i-1,j)+perm(i,j))
                    a(i) = - (csqx*mnhli)
                    b(i) = coef2 + csqx*(plhli + mnhli)
                    c(i) = 0.
                    alpha = - (csqx*plhli)
                endif
            endif
        enddo
        
        if (debug == 1 .and. j == 4) then
            do i = 1,mx
                print*, 'i,j = ',i,j
                print*,'a(i,j) =',a(i)
                print*,'b(i,j) =',b(i)
                print*,'c(i,j) =',c(i)
                print*,'r(i,j) =',r(i)
            enddo
            print*, 'alpha',alpha
            print*, 'beta',beta
            print*,'calling cyclic'
        endif

        call cyclic(a,b,c,alpha,beta,r,tempP,mx)
        if (debug == 1 .and. j == 4) print*, 'tempP = ('
        do i=1,mx
            if (debug == 1  .and. j == 4) print*,i,',',j,') = ', tempP(i)
            halfP(i,j)=tempP(i)
        enddo
        halfP(mx+1,j) = halfP(1,j)
    enddo

    !       second part of ADI update for the next time step n+1
    do j=1,my+1
        tempP(j) = 0.
        a(j) = 0.
        b(j) = 0.
        c(j) = 0.
        r(j) = 0.
    enddo


    do i = 1,mx
    do j=1,my+1
    coef2 = phimat(i,j)
    if (i /= 1 .and. i /=mx) then
        plhli = 2.*perm(i+1,j)*perm(i,j)/ &
        (perm(i+1,j)+perm(i,j))
        mnhli = 2.*perm(i-1,j)*perm(i,j)/ &
        (perm(i-1,j)+perm(i,j))
        r(j) = halfP(i-1,j)*(csqx* mnhli) + &
        halfP(i,j)* &
        (coef2-csqx*(plhli+mnhli))+  &
        halfP(i+1,j)*(csqx* plhli) - &
        cdiv*divvel(i,j)
    else
        if (i == 1)then
        plhli = 2.*perm(i+1,j)*perm(i,j)/ &
        (perm(i+1,j)+perm(i,j))
        mnhli = 2.*perm(mx,j)*perm(i,j)/ &
        (perm(mx,j)+perm(i,j))
        r(j) = halfP(mx,j)*( csqx* mnhli) + &
        halfP(i,j)* &
        (coef2-csqx*(plhli+mnhli))+  &
        halfP(i+1,j)*(csqx* plhli)- &
        cdiv*divvel(i,j)
        else
        plhli = 2.*perm(1,j)*perm(i,j)/ &
        (perm(1,j)+perm(i,j))
        mnhli = 2.*perm(i-1,j)*perm(i,j)/ &
        (perm(i-1,j)+perm(i,j))
        r(j) = halfP(i-1,j)*( csqx* mnhli) + &
        halfP(i,j)* &
        (coef2-csqx*(plhli+mnhli))+  &
        halfP(1,j)*( csqx* plhli) - &
        cdiv*divvel(i,j)
        endif
    endif
    
    if(j /= 1 .and. j /= my+1) then
        plhlj = 2.*perm(i,j+1)*perm(i,j)/ &
        (perm(i,j+1)+perm(i,j))
        mnhlj = 2.*perm(i,j-1)*perm(i,j)/ &
        (perm(i,j-1)+perm(i,j))
        a(j) =  - (csqy*mnhlj)
        b(j) = coef2 + csqy*(plhlj + mnhlj)
        c(j) =  - (csqy*plhlj)
    else
    
    if (j == 1)then
        plhlj = 2.*perm(i,j+1)*perm(i,j)/ &
        (perm(i,j+1)+perm(i,j))
        a(j)= 0.
        b(j)= coef2 + csqy*plhlj
        c(j) = - csqy*plhlj
    else
        mnhlj = 2.*perm(i,j-1)*perm(i,j)/ &
        (perm(i,j-1)+perm(i,j))
        a(j) = - csqy*mnhlj
        b(j) = coef2 + csqy*mnhlj
        c(j) = 0.
    endif

    endif
    enddo

    call tridag(a,b,c,r,tempP,my+1)
    do j=1,my+1
        press(i,j)=tempP(j)
    enddo
    enddo

    do j = 1,my+1
    press(mx+1,j)=press(1,j)
    enddo
    endif

    ! %%%%%%%%%%%%%%%% DIRICHLET BOUNDARY CONDITION %%%%%%%%%%%%%%%%%%%%%%%%%%%

    if (bc == 1 .or. bc == 3) then

      if (yn_explicit) then

        do i=1,mx
            do j=2,mywl
               halfP(i,j)=0.
            enddo
            press(i,1) = press_bot(i)
            press(i,mywl+1) = ptop
            halfP(i,1)=press(i,1)
            halfP(i,mywl+1)=press(i,mywl+1)
        enddo

        do j=2,mywl
            do i = 1,mx
             coef2 = phimat(i,j)
             plhlj = 2.*perm(i,j+1)*perm(i,j)/ &
                      (perm(i,j+1)+perm(i,j))
             mnhlj = 2.*perm(i,j-1)*perm(i,j)/ &
                      (perm(i,j-1)+perm(i,j))
             halfP(i,j) =  press(i,j-1)*csqy*mnhlj + &
                       press(i,j)*(coef2 - csqy*(plhlj+mnhlj)) + &
                       press(i,j+1)*csqy*plhlj - &
                       cdiv*divvel(i,j)
             if (yn_advect) then
               halfP(i,j) = halfP(i,j) - &
                press(i,j-1)*0.5d0*cady*(1-coef2)*sitevy(i,j) + &
                press(i,j)* &
                 0.5d0*cady*(1-coef2)*(sitevy(i,j)-sitevy(i,j+1)) + &
                press(i,j+1)*0.5d0*cady*(1-coef2)*sitevy(i,j+1)
             endif

             if (i /= 1 .and. i /= mx) then
               plhli = 2.*perm(i+1,j)*perm(i,j)/ &
                        (perm(i+1,j)+perm(i,j))
               mnhli = 2.*perm(i-1,j)*perm(i,j)/ &
                        (perm(i-1,j)+perm(i,j))
               halfP(i,j) = halfP(i,j) + &
                         press(i-1,j)*csqx*mnhli - &
                         press(i,j)*csqx*(plhli+mnhli) + &
                         press(i+1,j)*csqx*plhli
               if (yn_advect) then
                 halfP(i,j) = halfP(i,j) - &
                  press(i-1,j)*0.5d0*cadx*(1-coef2)*sitevx(i-1,j) + &
                  press(i,j)* &
                    0.5d0*cadx*(1-coef2)*(sitevx(i-1,j)-sitevx(i,j)) + &
                  press(i+1,j)*0.5d0*cadx*(1-coef2)*sitevx(i,j)
               endif
             else
               if (i == 1) then
                 plhli = 2.*perm(i+1,j)*perm(i,j)/ &
                          (perm(i+1,j)+perm(i,j))
                 mnhli = 2.*perm(mx,j)*perm(i,j)/ &
                          (perm(mx,j)+perm(i,j))
                 halfP(i,j) = halfP(i,j) + &
                         press(mx,j)*csqx*mnhli - &
                         press(i,j)*csqx*(plhli+mnhli) + &
                         press(i+1,j)*csqx*plhli
                 if (yn_advect) then
                   halfP(i,j) = halfP(i,j) - &
                    press(mx,j)*0.5d0*cadx*(1-coef2)*sitevx(mx,j) + &
                    press(i,j)* &
                      0.5d0*cadx*(1-coef2)*(sitevx(mx,j)-sitevx(i,j)) + &
                    press(i+1,j)*0.5d0*cadx*(1-coef2)*sitevx(i,j)
                 endif
               else
                 plhli = 2.*perm(1,j)*perm(i,j)/ &
                          (perm(1,j)+perm(i,j))
                 mnhli = 2.*perm(i-1,j)*perm(i,j)/ &
                          (perm(i-1,j)+perm(i,j))
                 halfP(i,j) = halfP(i,j) + &
                         press(i-1,j)*csqx*mnhli - &
                         press(i,j)*csqx*(plhli+mnhli) + &
                         press(1,j)*csqx*plhli
                 if (yn_advect) then
                   halfP(i,j) = halfP(i,j) - &
                    press(i-1,j)*0.5d0*cadx*(1-coef2)*sitevx(i-1,j) + &
                    press(i,j)* &
                     0.5d0*cadx*(1-coef2)*(sitevx(i-1,j)-sitevx(i,j)) + &
                    press(1,j)*0.5d0*cadx*(1-coef2)*sitevx(i,j)
                 endif
               endif
             endif

           enddo
         enddo

         do j=2,mywl
            do i = 1,mx
                press(i,j) = halfP(i,j)/phimat(i,j)
            enddo
         enddo
         do j=1,mywl+1
            press(mx+1,j) = press(1,j)
         enddo

      else  ! not yn_explicit

        do i=1,mx
            do j=2,mywl
                halfP(i,j)=0.
            enddo
            tempP(i) = 0.
            a(i) = 0.
            b(i) = 0.
            c(i) = 0.
            r(i) = 0.
            press(i,1) = press_bot(i)
            press(i,mywl+1) = ptop
            halfP(i,1)=press(i,1)
            halfP(i,mywl+1)=press(i,mywl+1)
        enddo
        !print*, 'XX u1: halfP = ',halfP(1:mx+1,my+1),press(1:mx+1,my+1)

        do j=2,mywl
            do i = 1,mx
            coef2 = phimat(i,j)
            plhlj = 2.*perm(i,j+1)*perm(i,j)/ &
            (perm(i,j+1)+perm(i,j))
            mnhlj = 2.*perm(i,j-1)*perm(i,j)/ &
            (perm(i,j-1)+perm(i,j))
            r(i) =  press(i,j-1)*(csqy*mnhlj) + &
            press(i,j)* &
            (coef2 - csqy*(plhlj+mnhlj)) + &
            press(i,j+1)*(csqy*plhlj) - &
            cdiv*divvel(i,j)
!            if (j==mywl) then
!                print*, 'XX r-:',i,r(i),press(i,j-1),press(i,j),press(i,j+1)
!                print*, 'XX r-:',i,csqy,mnhlj,plhlj,coef2
!            endif
            if (yn_advect) then
                r(i)=r(i) - &
                press(i,j-1)*0.5d0*cady*(1-coef2)*sitevy(i,j) + &
                press(i,j)* &
                 0.5d0*cady*(1-coef2)*(sitevy(i,j)-sitevy(i,j+1)) + &
                press(i,j+1)*0.5d0*cady*(1-coef2)*sitevy(i,j+1)
!                if (j==mywl) then
!                    print*, 'XX r+:',i,r(i),press(i,j-1),press(i,j),press(i,j+1)
!                    print*, 'XX r+:',i,cady,coef2,sitevy(i,j),sitevy(i,j+1)
!                endif
            endif

            if (i /= 1 .and. i /= mx) then
                plhli = 2.*perm(i+1,j)*perm(i,j)/ &
                (perm(i+1,j)+perm(i,j))
                mnhli = 2.*perm(i-1,j)*perm(i,j)/ &
                (perm(i-1,j)+perm(i,j))
                a(i) =  - (csqx*mnhli)
                b(i) = coef2 + csqx*(plhli + mnhli)
                c(i) =  - (csqx*plhli)
                if (yn_advect) then
                    a(i) = a(i) + 0.5d0*cadx*(1-coef2)*sitevx(i-1,j)
                    b(i) = b(i) - &
                     0.5d0*cadx*(1-coef2)*(sitevx(i-1,j)-sitevx(i,j))
                    c(i) = c(i) - 0.5d0*cadx*(1-coef2)*sitevx(i,j)
                endif
            else
                !       LEFT BOUNDARY
                if (i == 1)then
                    plhli = 2.*perm(i+1,j)*perm(i,j)/ &
                    (perm(i+1,j)+perm(i,j))
                    mnhli = 2.*perm(mx,j)*perm(i,j)/ &
                    (perm(mx,j)+perm(i,j))
                    a(i) = 0.
                    beta =  - (csqx*mnhli)
                    b(i) = coef2 + csqx*(plhli + mnhli)
                    c(1) =  - (csqx*plhli)
                   if (yn_advect) then
                     beta = beta + 0.5d0*cadx*(1-coef2)*sitevx(mx,j)
                     b(i) = b(i) - &
                      0.5d0*cadx*(1-coef2)*(sitevx(mx,j)-sitevx(i,j))
                     c(i) = c(i) - 0.5d0*cadx*(1-coef2)*sitevx(i,j)
                   endif
                    !       RIGHT BOUNDARY
                else
                    plhli = 2.*perm(1,j)*perm(i,j)/ &
                    (perm(1,j)+perm(i,j))
                    mnhli = 2.*perm(i-1,j)*perm(i,j)/ &
                    (perm(i-1,j)+perm(i,j))
                    a(i) = - (csqx*mnhli)
                    b(i) = coef2 + csqx*(plhli + mnhli)
                    c(i) = 0.
                    alpha = - (csqx*plhli)
                   if (yn_advect) then
                     a(i) = a(i) + 0.5d0*cadx*(1-coef2)*sitevx(i-1,j)
                     b(i) = b(i) - &
                      0.5d0*cadx*(1-coef2)*(sitevx(i-1,j)-sitevx(i,j))
                     alpha = alpha - 0.5d0*cadx*(1-coef2)*sitevx(i,j)
                   endif
                endif
            endif
            enddo

!            if (j==mywl) then
!                print*, 'XX a:',a
!                print*, 'XX b:',b
!                print*, 'XX c:',c
!                print*, 'XX alpha:',alpha
!                print*, 'XX beta:',beta
!                print*, 'XX r:',r
!                print*, 'XX tempP:',tempP
!                print*, 'XX mx:',mx
!            endif
            call cyclic(a,b,c,alpha,beta,r,tempP,mx)
        
            do i=1,mx
                halfP(i,j)=tempP(i)
            enddo
            halfP(mx+1,j)=halfP(1,j)
            !print*, 'XX u2',halfP(1:mx+1,my)
        enddo

        do j=1,mywl+1
            tempP(j) = 0.
            a(j) = 0.
            b(j) = 0.
            c(j) = 0.
            r(j) = 0.
        enddo
        
        do i = 1,mx
            do j=2,mywl
                coef2 = phimat(i,j)
                if (i /= 1 .and. i /=mx) then
                    plhli = 2.*perm(i+1,j)*perm(i,j)/ &
                    (perm(i+1,j)+perm(i,j))
                    mnhli = 2.*perm(i-1,j)*perm(i,j)/ &
                    (perm(i-1,j)+perm(i,j))
                    r(j) = halfP(i-1,j)*(csqx* mnhli) + &
                    halfP(i,j)*(coef2-csqx*(plhli+mnhli))+ &
                    halfP(i+1,j)*( csqx* plhli) - &
                    cdiv*divvel(i,j)
                    if (yn_advect) then
                      r(j) = r(j) - &
                      halfP(i-1,j)*0.5d0*cadx*(1-coef2)*sitevx(i-1,j) + &
                      halfP(i,j)* &
                       0.5d0*cadx*(1-coef2)*(sitevx(i-1,j)-sitevx(i,j)) + &
                      halfP(i+1,j)*0.5d0*cadx*(1-coef2)*sitevx(i,j)
                    endif
                else
                    if (i == 1)then
                        plhli = 2.*perm(i+1,j)*perm(i,j)/ &
                        (perm(i+1,j)+perm(i,j))
                        mnhli = 2.*perm(mx,j)*perm(i,j)/ &
                        (perm(mx,j)+perm(i,j))
                        r(j) = halfP(mx,j)*( csqx* mnhli) + &
                        halfP(i,j)*(coef2 - csqx*(plhli+mnhli)) + &
                        halfP(i+1,j)*( csqx* plhli) - &
                        cdiv*divvel(i,j)
                       if (yn_advect) then
                         r(j) = r(j) - &
                         halfP(mx,j)*0.5d0*cadx*(1-coef2)*sitevx(mx,j) + &
                         halfP(i,j)* &
                          0.5d0*cadx*(1-coef2)*(sitevx(mx,j)-sitevx(i,j)) + &
                         halfP(i+1,j)*0.5d0*cadx*(1-coef2)*sitevx(i,j)
                       endif
                    else
                        plhli = 2.*perm(1,j)*perm(i,j)/ &
                        (perm(1,j)+perm(i,j))
                        mnhli = 2.*perm(i-1,j)*perm(i,j)/ &
                        (perm(i-1,j)+perm(i,j))
                        r(j) = halfP(i-1,j)*( csqx* mnhli) + &
                        halfP(i,j)*(coef2 - csqx*(plhli+mnhli)) + &
                        halfP(1,j)*( csqx* plhli) - &
                        cdiv*divvel(i,j)
                       if (yn_advect) then
                         r(j) = r(j) - &
                         halfP(i-1,j)*0.5d0*cadx*(1-coef2)*sitevx(i-1,j) + &
                         halfP(i,j)* &
                          0.5d0*cadx*(1-coef2)*(sitevx(i-1,j)-sitevx(i,j)) + &
                         halfP(1,j)*0.5d0*cadx*(1-coef2)*sitevx(i,j)
                       endif
                    endif
                endif

                plhlj = 2.*perm(i,j+1)*perm(i,j)/ &
                (perm(i,j+1)+perm(i,j))
                mnhlj = 2.*perm(i,j-1)*perm(i,j)/ &
                (perm(i,j-1)+perm(i,j))
                a(j) =  - (csqy*mnhlj)
                b(j) = coef2 + csqy*(plhlj + mnhlj)
                c(j) =  - (csqy*plhlj)
                if (yn_advect) then
                    a(j) = a(j) + 0.5d0*cady*(1-coef2)*sitevy(i,j)
                    b(j) = b(j) - &
                     0.5d0*cady*(1-coef2)*(sitevy(i,j)-sitevy(i,j+1))
                    c(j) = c(j) - 0.5d0*cady*(1-coef2)*sitevy(i,j+1)
                endif
            enddo
            
            r(1) = halfP(i,1)
            a(1) = 0
            b(1) = 1
            c(1) = 0
            r(mywl+1) = halfP(i,mywl+1)
            a(mywl+1) = 0
            b(mywl+1) = 1
            c(mywl+1) = 0
            !print*, 'XX u3-:',press(1:mx+1,my)
            call tridag(a,b,c,r,tempP,mywl+1)
            do j=1,mywl+1
                press(i,j)=tempP(j)
            enddo
            !print*, 'XX u3+:',press(1:mx+1,my)
        enddo

        do j = 1,mywl+1
            press(mx+1,j)=press(1,j)
        enddo

      endif ! yn_explicit
    endif

    ! %%%%%%%%%%%%%%%% MIXED B.C. DIRICHLET TOP AND NEUMANN BOTTOM %%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! %%%%  at the bottom 0 pressure and at the top 0 pressure gradient %%%%%%%%%%%%%%%%%%%%
    if (bc == 4) then
        do i=1,mx
            do j=2,my+1
                halfP(i,j)=0.
            enddo
            tempP(i) = 0.
            a(i) = 0.
            b(i) = 0.
            c(i) = 0.
            r(i) = 0.
            press(i,1) = ptop   ! called ptop but is set at the bottom
            halfP(i,1)=press(i,1)
        enddo

        do j=2,my+1
          do i = 1,mx
            coef2 = phimat(i,j)
            if (j /= my+1) then
                plhlj =2.*perm(i,j+1)*perm(i,j)/ &
                (perm(i,j+1)+perm(i,j))
                mnhlj = 2.*perm(i,j-1)*perm(i,j)/ &
                (perm(i,j-1)+perm(i,j))
                r(i)= press(i,j-1)*(csqy*mnhlj)+ &
                press(i,j)* &
                (coef2-csqy*(plhlj+mnhlj))+  &
                press(i,j+1)*(csqy*plhlj) - &
                cdiv*divvel(i,j)
            else
                plhlj = 0.
                mnhlj = 2.*perm(i,j-1)*perm(i,j)/ &
                (perm(i,j-1)+perm(i,j))
                r(i)= press(i,j-1)*(csqy*mnhlj)+ &
                press(i,j)* &
                (coef2-csqy*mnhlj)+  &
                0. - &
                cdiv*divvel(i,j)
            endif

            if (i /= 1 .and. i /= mx) then
                plhli = 2.*perm(i+1,j)*perm(i,j)/ &
                (perm(i+1,j)+perm(i,j))
                mnhli = 2.*perm(i-1,j)*perm(i,j)/ &
                (perm(i-1,j)+perm(i,j))
                a(i) =  - (csqx*mnhli)
                b(i) = coef2 + csqx*(plhli + mnhli)
                c(i) =  - (csqx*plhli)
            else
                !       LEFT BOUNDARY
                if (i == 1)then
                    plhli = 2.*perm(i+1,j)*perm(i,j)/ &
                    (perm(i+1,j)+perm(i,j))
                    mnhli = 2.*perm(mx,j)*perm(i,j)/ &
                    (perm(mx,j)+perm(i,j))
                    a(i) = 0.
                    beta =  - (csqx*mnhli)
                    b(i) = coef2 + csqx*(plhli + mnhli)
                    c(1) =  - (csqx*plhli)
                    !       RIGHT BOUNDARY
                else
                    plhli = 2.*perm(1,j)*perm(i,j)/ &
                    (perm(1,j)+perm(i,j))
                    mnhli = 2.*perm(i-1,j)*perm(i,j)/ &
                    (perm(i-1,j)+perm(i,j))
                    a(i) = - (csqx*mnhli)
                    b(i) = coef2 + csqx*(plhli + mnhli)
                    c(i) = 0.
                    alpha = - (csqx*plhli)
                endif
            endif
          enddo

          call cyclic(a,b,c,alpha,beta,r,tempP,mx)

          do i=1,mx
            halfP(i,j)=tempP(i)
          enddo
          halfP(mx+1,j) = halfP(1,j)
        enddo

        !       second part of ADI update for the next time step n+1
        do j=1,my+1
            tempP(j) = 0.
            a(j) = 0.
            b(j) = 0.
            c(j) = 0.
            r(j) = 0.
        enddo

        do i = 1,mx
          do j=2,my+1
            coef2 = phimat(i,j)
            if (i /= 1 .and. i /=mx) then
                plhli = 2.*perm(i+1,j)*perm(i,j)/ &
                (perm(i+1,j)+perm(i,j))
                mnhli = 2.*perm(i-1,j)*perm(i,j)/ &
                (perm(i-1,j)+perm(i,j))
                r(j) = halfP(i-1,j)*(csqx* mnhli) + &
                halfP(i,j)* &
                (coef2-csqx*(plhli+mnhli))+  &
                halfP(i+1,j)*(csqx* plhli) - &
                cdiv*divvel(i,j)
            else
              if (i == 1)then
                plhli = 2.*perm(i+1,j)*perm(i,j)/ &
                (perm(i+1,j)+perm(i,j))
                mnhli = 2.*perm(mx,j)*perm(i,j)/ &
                (perm(mx,j)+perm(i,j))
                r(j) = halfP(mx,j)*( csqx* mnhli) + &
                halfP(i,j)* &
                (coef2-csqx*(plhli+mnhli))+  &
                halfP(i+1,j)*(csqx* plhli)- &
                cdiv*divvel(i,j)
              else
                plhli = 2.*perm(1,j)*perm(i,j)/ &
                (perm(1,j)+perm(i,j))
                mnhli = 2.*perm(i-1,j)*perm(i,j)/ &
                (perm(i-1,j)+perm(i,j))
                r(j) = halfP(i-1,j)*( csqx* mnhli) + &
                halfP(i,j)* &
                (coef2-csqx*(plhli+mnhli))+  &
                halfP(1,j)*( csqx* plhli) - &
                cdiv*divvel(i,j)
              endif
            endif

            if(j /= my+1) then
                plhlj = 2.*perm(i,j+1)*perm(i,j)/ &
                (perm(i,j+1)+perm(i,j))
                mnhlj = 2.*perm(i,j-1)*perm(i,j)/ &
                (perm(i,j-1)+perm(i,j))
                a(j) =  - (csqy*mnhlj)
                b(j) = coef2 + csqy*(plhlj + mnhlj)
                c(j) =  - (csqy*plhlj)
            else
                mnhlj = 2.*perm(i,j-1)*perm(i,j)/ &
                (perm(i,j-1)+perm(i,j))
                a(j) = - csqy*mnhlj
                b(j) = coef2 + csqy*mnhlj
                c(j) = 0.
            endif
          enddo

          r(1) = halfP(i,1)
          a(1) = 0
          b(1) = 1
          c(1) = 0

          call tridag(a,b,c,r,tempP,my+1)
          do j=1,my+1
            press(i,j)=tempP(j)
          enddo
        enddo

        do j = 1,my+1
          press(mx+1,j)=press(1,j)
        enddo
    endif

    ! %%%%  at the bottom 0 pressure gradient and at the top 0 pressure %%%%%%%%%%%%%%%%%%%%

    if (bc == 5) then
        do i=1,mx
            do j=1,mywl
                halfP(i,j)=0.
            enddo
            tempP(i) = 0.
            a(i) = 0.
            b(i) = 0.
            c(i) = 0.
            r(i) = 0.
            press(i,mywl+1) = ptop
            halfP(i,mywl+1)=press(i,mywl+1)

        enddo
        do j=1,mywl
            do i = 1,mx
                coef2 = phimat(i,j)
                if (j /= 1) then
                    plhlj = 2.*perm(i,j+1)*perm(i,j)/ &
                    (perm(i,j+1)+perm(i,j))
                    mnhlj = 2.*perm(i,j-1)*perm(i,j)/ &
                    (perm(i,j-1)+perm(i,j))

                    r(i)= press(i,j-1)*(csqy*mnhlj)+ &
                    press(i,j)* &
                    (coef2-csqy*(plhlj+mnhlj))+  &
                    press(i,j+1)*(csqy*plhlj) - &
                    cdiv*divvel(i,j)
                else
                    !       BOTTOM BOUNDARY
                    plhlj = 2.*perm(i,j+1)*perm(i,j)/ &
                    (perm(i,j+1)+perm(i,j))
                    mnhlj = 0
                    r(i) = 0. + &
                    press(i,j)* &
                    (coef2-csqy*plhlj)+ &
                    press(i,j+1)*(csqy*plhlj) - &
                    cdiv*divvel(i,j)
                endif
                if (i /= 1 .and. i /= mx) then
                    plhli = 2.*perm(i+1,j)*perm(i,j)/ &
                    (perm(i+1,j)+perm(i,j))
                    mnhli = 2.*perm(i-1,j)*perm(i,j)/ &
                    (perm(i-1,j)+perm(i,j))

                    a(i) =  - (csqx*mnhli)
                    b(i) = coef2 + csqx*(plhli + mnhli)
                    c(i) =  - (csqx*plhli)
                else
                    !       LEFT BOUNDARY
                    if (i == 1)then
                        plhli = 2.*perm(i+1,j)*perm(i,j)/ &
                        (perm(i+1,j)+perm(i,j))
                        mnhli = 2.*perm(mx,j)*perm(i,j)/ &
                        (perm(mx,j)+perm(i,j))
                        a(i) = 0.
                        beta =  - (csqx*mnhli)
                        b(i) = coef2 + csqx*(plhli + mnhli)
                        c(1) =  - (csqx*plhli)
                        !       RIGHT BOUNDARY
                    else
                        plhli = 2.*perm(1,j)*perm(i,j)/ &
                        (perm(1,j)+perm(i,j))
                        mnhli = 2.*perm(i-1,j)*perm(i,j)/ &
                        (perm(i-1,j)+perm(i,j))
                        a(i) = - (csqx*mnhli)
                        b(i) = coef2 + csqx*(plhli + mnhli)
                        c(i) = 0.
                        alpha = - (csqx*plhli)
                    endif
                endif
            enddo


            call cyclic(a,b,c,alpha,beta,r,tempP,mx)
            if (debug == 1 .and. j == 4) print*, 'tempP = ('
            
            do i=1,mx
                if (debug == 1  .and. j == 4) print*,i,',',j,') = ', tempP(i)
                halfP(i,j)=tempP(i)
            enddo
            
            halfP(mx+1,j) = halfP(1,j)
        enddo

        !       second part of AFI update for the next time step n+1
        do j=1,mywl+1
            tempP(j) = 0.
            a(j) = 0.
            b(j) = 0.
            c(j) = 0.
            r(j) = 0.
        enddo


        do i = 1,mx
            do j=1,mywl
                coef2 = phimat(i,j)
                if (i /= 1 .and. i /=mx) then
                    plhli = 2.*perm(i+1,j)*perm(i,j)/ &
                    (perm(i+1,j)+perm(i,j))
                    mnhli = 2.*perm(i-1,j)*perm(i,j)/ &
                    (perm(i-1,j)+perm(i,j))
                    r(j) = halfP(i-1,j)*(csqx* mnhli) + &
                    halfP(i,j)* &
                    (coef2-csqx*(plhli+mnhli))+  &
                    halfP(i+1,j)*(csqx* plhli) - &
                    cdiv*divvel(i,j)
                else
                    if (i == 1)then
                        plhli = 2.*perm(i+1,j)*perm(i,j)/ &
                        (perm(i+1,j)+perm(i,j))
                        mnhli = 2.*perm(mx,j)*perm(i,j)/ &
                        (perm(mx,j)+perm(i,j))
                        r(j) = halfP(mx,j)*( csqx* mnhli) + &
                        halfP(i,j)* &
                        (coef2-csqx*(plhli+mnhli))+  &
                        halfP(i+1,j)*(csqx* plhli)- &
                        cdiv*divvel(i,j)
                    else
                        plhli = 2.*perm(1,j)*perm(i,j)/ &
                        (perm(1,j)+perm(i,j))
                        mnhli = 2.*perm(i-1,j)*perm(i,j)/ &
                        (perm(i-1,j)+perm(i,j))
                        r(j) = halfP(i-1,j)*( csqx* mnhli) + &
                        halfP(i,j)* &
                        (coef2-csqx*(plhli+mnhli))+  &
                        halfP(1,j)*( csqx* plhli) - &
                        cdiv*divvel(i,j)
                    endif
                endif
                
                if(j /= 1) then
                    plhlj = 2.*perm(i,j+1)*perm(i,j)/ &
                    (perm(i,j+1)+perm(i,j))
                    mnhlj = 2.*perm(i,j-1)*perm(i,j)/ &
                    (perm(i,j-1)+perm(i,j))
                    a(j) =  - (csqy*mnhlj)
                    b(j) = coef2 + csqy*(plhlj + mnhlj)
                    c(j) =  - (csqy*plhlj)
                else
                    ! bottom boundary
                    plhlj = 2.*perm(i,j+1)*perm(i,j)/ &
                    (perm(i,j+1)+perm(i,j))
                    a(j)= 0.
                    b(j)= coef2 + csqy*plhlj
                    c(j) = - csqy*plhlj
                endif
            enddo
            r(mywl+1) = halfP(i,mywl+1)
            a(mywl+1) = 0
            b(mywl+1) = 1
            c(mywl+1) = 0
            call tridag(a,b,c,r,tempP,mywl+1)
            do j=1,mywl+1
                press(i,j)=tempP(j)
            enddo
        enddo

        do j = 1,mywl+1
            press(mx+1,j)=press(1,j)
        enddo
    endif

    if (bc == 6) then
        do i=1,mx+1
            do j=1,mywl+1
                press(i,j)=pbot
            enddo
        enddo
    endif

!    print*, 'XX u4:',press(1:mx+1,my)
!    print*, 'XX u5:',press(1:mx+1,my+1)
!    print*, 'XX u6:',press(1:mx+1,my+2)
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!    do j = 1,my+1
!        do i = 2,mx
!            gradpx(i,j) = (press(i+1,j) - press(i-1,j))/ &
!            (2.*dxc)
!        enddo
!        gradpx(1,j) = (press(2,j) - press(mx,j))/ &
!        (2.*dxc)
!        gradpx(mx+1,j) = gradpx(1,j)
!    enddo
    do j = 1,my+1
        do i = 1,mx-1
            gradpx(i,j) = (press(i+1,j) - press(i,j))/dxc
            plhli = 2.*perm(i+1,j)*perm(i,j)/(perm(i+1,j)+perm(i,j))
            vol_flux_x(i,j) = sitevx(i,j) - plhli/visco*gradpx(i,j)
        enddo
        gradpx(mx,j) = (press(1,j) - press(mx,j))/dxc
        plhli = 2.*perm(1,j)*perm(mx,j)/(perm(1,j)+perm(mx,j))
        vol_flux_x(mx,j) = sitevx(mx,j) - plhli/visco*gradpx(mx,j)
        gradpx(mx+1,j) = gradpx(1,j)   ! for easier implementation of forces in force module
    enddo

    do i = 1,mx + 1
!        do j = 2,my
!            gradpy(i,j) = (press(i,j+1)-press(i,j-1))/ &
!            (2.*dyc)
!        enddo
        do j = 2,my+1
            gradpy(i,j) = (press(i,j)-press(i,j-1))/dyc
            plhli = 2.*perm(i,j)*perm(i,j-1)/(perm(i,j)+perm(i,j-1))
            vol_flux_y(i,j) = sitevy(i,j) - plhli/visco*gradpy(i,j)
        enddo
        gradpy(i,1) = 0d0
        vol_flux_y(i,1) = 0d0
        gradpy(i,my+2) = 0d0
        vol_flux_y(i,my+2) = 0d0

!        !       neumann = zero flux b.c.
!        if (bcl == 2) then
!            gradpy(i,1) = (press(i,2) - press(i,1))/(dyc)
!            gradpy(i,my+1)=(press(i,my+1)-press(i,my))/(dyc)
!        else
!            if (bcl == 1) then
!                !       dirichlet
!                gradpy(i,1) =  (press(i,2) - press(i,1))/(dyc)
!                gradpy(i,my+1) = (press(i,my+1) - press(i,my))/(dyc)
!            else
!                !       mixed
!                gradpy(i,1) = (press(i,2) - press(i,1))/(dyc)
!                gradpy(i,my+1) = (press(i,my+1) - press(i,my))/(dyc)
!            endif
!        endif
    enddo

    ! CALCULATE VOLUME FLUX FROM EACH CELL
    max_flux_x = 0d0
    do j = 1,my+1
        diff_vol_flux_x(1,j)=(vol_flux_x(1,j)-vol_flux_x(mx,j))/dxc
        max_flux_x=max(max_flux_x,abs(diff_vol_flux_x(1,j)))
        do i = 2,mx
            diff_vol_flux_x(i,j)=(vol_flux_x(i,j)-vol_flux_x(i-1,j))/dxc
            max_flux_x=max(max_flux_x,abs(diff_vol_flux_x(i,j)))
        enddo
    enddo

    max_flux_y = 0d0
    do i = 1,mx
        do j = 1,my+1
            diff_vol_flux_y(i,j)=(vol_flux_y(i,j+1)-vol_flux_y(i,j))/dyc
            max_flux_y=max(max_flux_y,abs(diff_vol_flux_y(i,j)))
        enddo
    enddo

    do j = 1,my+1
        do i = 1,mx
            pressave(i,j)=pressave(i,j)+press(i,j)
        enddo
    enddo
    
end subroutine

end module
