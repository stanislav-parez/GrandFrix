module mod_rw_conf
    use mycommons 
    use mod_pssubs
    use mod_inputReader
    use mod_outputWriter
    implicit none

    integer,public :: shearcolor
    integer :: restart_id = 0
    
    real(8) :: rd(256,2) = 0.
    real(8) :: bl(256,2) = 0.
    
    real(8),public :: afrac1,ftsum,c1,c2,prop1,prop2
    real(8),public :: diffsum1,diffsum2,diffrac1,diffrac2

    private

    ! public procedures
    public :: write_restart, pscriptplot, pscriptplot2, &
              readrestart

    contains
!      SUBROUTINE READCN READS BINARY UNFOMATTED DATA OF CNFILE
!
!       SUBROUTINE WRITCN WRITES BINARY UNFOMATTED DATA IN CNFILE,
!     IF FORM=1. ASCII IF FORM=0.
!
!     SUBROUTINE SAVE  SAVES CURRENT CONFIGURATION FOR FUTURE CHECKING. 
!=======================================================================
SUBROUTINE readrestart(title,ifluid) 

    integer i,k
    integer ifluid
    INTEGER :: CNUNIT,stat
    
    real(8) topg
    
    character(len=50) title
    character(len=7) ident

!    *******************************************************************
!    ** SUBROUTINE TO READ IN INITIAL CONFIGURATION FROM UNIT 10      **
!    *******************************************************************
    
    open(newunit=CNUNIT,file=trim(title),form='unformatted')
    rewind(cnunit)
    read(cnunit) ident
    read(cnunit) n,nbound(1:4)
    read(cnunit) ybot_out,ytop_out,xleft_out,xright_out
    read(cnunit) tau,tauadd
    
    if (ident == 'flforce') then
        read(cnunit) (rx(i),ry(i),radius(i),rnt(i),vx(i),vy(i),w(i), &
        fx(i),fy(i),tq(i),bdgrain(i),i=1,n)
    else
        read(cnunit) (rx(i),ry(i),radius(i),rnt(i),vx(i),vy(i),w(i), &
        fx(i),fy(i),tq(i),bdgrain(i),gtype(i),i=1,n)
    endif
    read(cnunit) (color(i),i=1,n)

    !  if reading old style n*n slip array, then convert to
    !   new style

    read(cnunit) contactknt
    do k=1,contactknt
        read(cnunit) contacti(k),contactj(k),contfn(k),contft(k)
    enddo

    read(cnunit,iostat=stat) dspring,dwall
    !if (stat /= 0) exit

    !if((ifluid == 0).and.(irstart/=2))then
        !do k = 1,MAXD+1
        do k = 1,my+2
            read(cnunit) (press(i,k),i=1,mx+1)
            read(cnunit) (gradpx(i,k),i=1,mx+1)
            read(cnunit) (gradpy(i,k),i=1,mx+1)
        enddo
    !endif
    
    !if (ib(2) == 5) read(cnunit) dspring,dwall,fspring
    
    CLOSE ( UNIT = CNUNIT )
    
    if((intruder==1).and.(radius(n)<(rintruder-1d-6))) then
        n=n+1
        radius(n)=rintruder
        
        if((rintruder+1)>(0.5d0*(xright-xleft)))then
            print*, 'intruder too big'
            stop
        endif
        
        rx(n)=0.5*(xleft+xright)+1d-4
        topg=ry(1)+radius(1)
        
        do i=nbound(4)+1,n-1
            topg=max(topg,ry(i)+radius(i))
        enddo         
        
        ry(n)=topg+radius(n)+1d0
        rnt(n)=0d0
        vx(n)=0d0
        vy(n)=0d0
        w(n)=0d0
        fx(n)=0d0
        fy(n)=0d0
        tq(n)=0d0
        bdgrain(n)=0
        gtype(n)=1
        color(n)=color(n-1)
    endif
    
end subroutine
!=======================================================================
SUBROUTINE write_restart(sufix)

    integer i,k
    INTEGER :: CNUNIT
    
    character(len=50) res_name
    character(len=*),intent(in),optional :: sufix
    character(len=7) ident

    if(present(sufix)) then
        write(res_name,'(2g0)') "restart_",sufix
    else
        restart_id = restart_id + 1
        write(res_name,'(2g0)') "restart_",restart_id
    endif
    
    OPEN ( newUNIT = CNUNIT, FILE = trim(res_name) , FORM = 'UNFORMATTED' )

    ident="gtype"
    write(cnunit) ident
    write(cnunit) n,(nbound(i),i=1,4)
    write(cnunit) ybot_out,ytop_out,xleft_out,xright_out
    write(cnunit) tau,tauadd

    write(cnunit) (rx(i),ry(i),radius(i),rnt(i),vx(i),vy(i),w(i), &
    fx(i),fy(i),tq(i),bdgrain(i),gtype(i),i=1,n)
    write(cnunit) (color(i),i=1,n)

    write(cnunit) contactknt
    
    do k=1,contactknt
        write(cnunit) contacti(k),contactj(k),contfn(k),contft(k)
    enddo

    write(cnunit) dspring,dwall

    do k = 1,my+2
        write(cnunit) (press(i,k),i=1,mx+1)
        write(cnunit) (gradpx(i,k),i=1,mx+1)
        write(cnunit) (gradpy(i,k),i=1,mx+1)
    enddo
    
    !if (ib(2) == 5) write(cnunit) dspring,dwall,fspring

    CLOSE ( UNIT = CNUNIT )

end subroutine
!=======================================================================
subroutine pscriptplot(next,igif,fnn,inext)
    !*******************************************************************
    !*******************************************************************
    !   DRAW A POSTSCRIPT PAGE !!!!
    !*******************************************************************
    !*******************************************************************

    integer logging
    integer psunit, i,j,k,index
    integer icnt,ikolor,next,n1,n2,n3,n4
    integer inext
    
    real(8) pmin,p1
    real(8) olap(MAXN*10),maxlap,minlap,wscale,width
    real(8) maxpr,pr(maxn)
    real(8) xl,sqrt,tforce
    real(8) amax,fmax,pmax
    real(8) gs,red,blue,green
    
    logical igif
    
    character(len=3) ann
    character(len=4) fnn
    character(len=30) pspage
    character(len=40) com1
    character(len=12) gfile
    character a1,a2,a3,a4
      
    22      format(e12.4,1x,e12.4,a)
    28      format(6e12.4,a)
    29      format(e12.4,a)


    ! set up red-blue color map 
    if (bl(1,1) == 0.) then
        do i=1,128
            rd(i,1) = .3+.7*dble(i-1)/127.
            bl(i,1) = .3+.7*(127.-dble(i-1))/127.
            rd(i,2) = .3+.7*dble(i-1)/127.
            bl(i,2) = .3+.7*(127.-dble(i-1))/127.
        enddo
    endif

    n1=(next/1000)
    n2=(next-n1*1000)/100
    n3=(next-n1*1000-n2*100)/10
    n4=(next-n1*1000-n2*100-n3*10)
    write(a1,fmt='(i1)') n1
    write(a2,fmt='(i1)') n2
    write(a3,fmt='(i1)') n3
    write(a4,fmt='(i1)') n4
    fnn = a1//a2//a3//a4
    
    if (inext == 0) then
        pspage='bpage'//fnn//".ps"
    else
        ann=char(inext+64)
        pspage='cpage'//fnn//ann
    endif
    
    ! start new postscript page and plot box  
    call newpsdoc(psunit,xleft_out,xright_out,ybot_out,ytop_out,pspage,0,tau)
    print*, 'height = ',ytop-ybot


    ! scale  forces and velocities by max. value at this time
    !  (not ideal, but what else to do?)
    amax=0.
    fmax=0.
    pmax=0.
    pmin=1.e+06

    xl = xright_out-xleft_out

    !   plot pressures and stress directions

    ! set the maximum pressure for scaling the color of particles
    ! i.e. all above pmax will be pink
    !            pmax=2.*(2.*presy)
    !            pmax=3.*(2.*presy)

    !  loop through all the contacts, find the overlaps

!    dxi = real(mx)/(xright-xleft)
!    dyi = real(my)/(ytop-ybot)
    maxlap=0.
    minlap=0.
    maxpr=0.
    
    do i=1,n
        pr(i)=0.
    enddo
    
    do k=1,contactknt
        i=contacti(k)
        j=contactj(k)
        pr(i)=pr(i)-contfn(k)
        pr(j)=pr(j)-contfn(k)
        tforce=contfn(k)*contfn(k)+contft(k)*contft(k)
        olap(k)=sqrt(tforce)

        !  FOR LINES REPRESENTING FRICTION

        maxlap=max(maxlap,olap(k))
    enddo
    
    do i=1,n
        maxpr=max(maxpr,pr(i))
    enddo

    do i=1,n
        if (gtype(i) /= 0) then
          if (bdgrain(i) == 0) then
            if (maxpr > 0.) then
                p1=pr(i)/maxpr
            else
                p1=0.
            endif
            if (p1 > 1.) then
                red=1.
                blue=.6
                green=.6
            else
                index =1+ max(0,nint(127*p1))
                red = rd(index,gtype(i))
                blue = bl(index,gtype(i))
                gs=blue-1.8*(1.-blue)
                green = max(0.0d0,gs)
            endif
          else
            red=.3
            blue=.3
            green=.3
          endif

            ikolor=1

            write(psunit,28)  radius(i), &
            rx(i),ry(i), &
            red,green,blue,' tbgrain'
            !if (ib(3) == -1) then
            !      duplicate atoms wrapping around x direction
            if (rx(i) > xright_out-radius(i)) then
                write(psunit,28)  radius(i), &
                rx(i)-xl,ry(i), &
                red,green,blue,' tbgrain'
            else if (rx(i) < xleft_out+radius(i)) then
                write(psunit,28)  radius(i), &
                rx(i)+xl,ry(i), &
                red,green,blue,' tbgrain'
            endif
            !endif

        endif
    enddo
    
    icnt=0
    ! draw the network of contacts
    write(psunit,*) '0 0 0 setrgbcolor' 

    wscale=0.25

    ! set logging=1 for run with gravity, takes the log of
    ! overlaps when drawing force lines
    logging=0
    if (logging == 1) then
        minlap=dlog10(minlap)
        maxlap=dlog10(maxlap)
    endif

    width=1.
    do k=1,contactknt
        if (olap(k) > 0.) then

        icnt=icnt+1
        i=contacti(k)
        j=contactj(k)

        if (maxlap /= minlap) then
            if (logging == 1) then
                width=dlog10(olap(k))/maxlap
            else
                width=olap(k)/maxlap
            endif
        endif

        if (width < 0.) width=0.
        
        write(psunit,29) wscale*width,' setlinewidth'

        if (abs(rx(j)-rx(i)) > .5*xl) then
            if (rx(i) < 0.) then
                write(psunit,22,advance='no') rx(i),ry(i),' moveto'
                write(psunit,22) rx(j)-xl,ry(j),' lineto stroke'
                write(psunit,22,advance='no') rx(i)+xl,ry(i),' moveto'
                write(psunit,22) rx(j),ry(j),' lineto stroke'
            else if (rx(i) > 0.) then
                write(psunit,22,advance='no') rx(i)-xl,ry(i),' moveto'
                write(psunit,22) rx(j),ry(j),' lineto stroke'
                write(psunit,22,advance='no') rx(i),ry(i),' moveto'
                write(psunit,22) rx(j)+xl,ry(j),' lineto stroke'
            endif
        else
            write(psunit,22,advance='no') rx(i),ry(i),' moveto'
            write(psunit,22) rx(j),ry(j),' lineto stroke'
            endif
        endif
    enddo

    print*,' number of contacts = ',icnt

    call endpspage(psunit)
    
    close(psunit)

    if (igif) then
        gfile = ""!" >& jpegjunk"
        com1='sh mkjpeg.sh '//trim(pspage)//trim(gfile)
        print*, com1
        call execute_command_line(com1)
    endif


end subroutine
!=======================================================================
subroutine pscriptplot2(next,igif,fnn)
    !*******************************************************************
    !*******************************************************************
    !  POSTSCRIPT PAGE: OPTION 2!!!!
    !*******************************************************************
    !*******************************************************************

    integer psunit, i,k,index
    integer icnt,ikolor,next,n1,n2,n3,n4
    integer nsum1,nsum2,g1,g2
    integer ksmall,kbig,nonzeroc
    integer np1,np2,nps1,nps2,ibnet

    real(8) amax,fmax,pmax
    real(8) gs,red,blue,green
    real(8) pmin,p1
    real(8) maxlap
    real(8) maxpr,pr(maxn)
    real(8) xl,sqrt
    real(8) gmn1,gmn2,gsum1,gsum2
    real(8) sigd(maxn),sigt(maxn),sign(maxn)
    real(8) ft,ftmax(maxn)
    real(8) bnet
    
    logical igif
    
    character(len=4) fnn
    character(len=30) pspage
    character(len=35) com1
    character(len=50) gfile
    character a1,a2,a3,a4

    28      format(3e12.4,3f7.4/2e11.3,a)
    
    
    ibnet = 0
    
    ! set up color map
    ! set up color map
    if (bl(1,1) == 0.) then
        do i=1,128
            rd(i,1) = .3+.7*dble(i-1)/127.
            bl(i,1) = .3+.7*(127.-dble(i-1))/127.
        enddo
    endif

    n1=(next/1000)
    n2=(next-n1*1000)/100
    n3=(next-n1*1000-n2*100)/10
    n4=(next-n1*1000-n2*100-n3*10)
    write(a1,fmt='(i1)') n1
    write(a2,fmt='(i1)') n2
    write(a3,fmt='(i1)') n3
    write(a4,fmt='(i1)') n4
    fnn = a1//a2//a3//a4
    pspage='espage'//fnn

    ! start new postscript page and plot box  
    call newpsdoc(psunit,xleft_out,xright_out,ybot_out,ytop_out,pspage,0,tau)
    print*, 'height = ',ytop-ybot

    ! scale  forces and velocities by max. value at this time
    !  (not ideal, but what else to do?)
    amax=0.
    fmax=0.
    pmax=0.
    g1=0.
    g2=0.
    pmin=1.e+06

    xl = xright_out-xleft_out

    if (higr < higr2) then
        do i=1,n
            if (2.*radius(i) > higr) then
                gtype(i)=2
            else
                gtype(i)=1
            endif
        enddo
    else
        do i=1,n
            if (2.*radius(i) > higr2) then
                gtype(i)=1
            else
                gtype(i)=2
            endif
        enddo
    endif

    !   plot differential stresses and stress directions

!    dxi = real(mx)/(xright-xleft)
!    dyi = real(my)/(ytop-ybot)
    maxlap=0.
    maxpr=0.
    
    do i=1,n
        sigd(i)=0.5*(sigmax(i)-sigmin(i))
        sigt(i)=-0.5*(sigmax(i)+sigmin(i))
        if (sigt(i) /= 0.) then
            sign(i)=sigd(i)/sigt(i)
        else
            sign(i)=1.
        endif
    enddo

    !  what do you want to plot? 
    do i=1,n
        pr(i)=sign(i)
        maxpr=max(pr(i),maxpr)
    enddo

    do i=1,n
        if (gtype(i) == 1) then
            call write_to('diffstress1.out', 2*radius(i),sigd(i),sigt(i),sign(i),sigmin(i), sigmax(i))
        endif
    enddo         
    call close_file('diffstress1.out') 
    
    do i=1,n
        if (gtype(i) == 2) then
            call write_to('diffstress2.out', 2*radius(i),sigd(i),sigt(i),sign(i),sigmin(i), sigmax(i))
        endif
    enddo          
    call close_file('diffstress2.out')

    gsum1=0.
    gsum2=0.
    nsum1=0
    nsum2=0
    do i=1,n
        if (gtype(i) == 1) then 
            gsum1=gsum1+pr(i)
            nsum1=nsum1+1
        else
            gsum2=gsum2+pr(i)
            nsum2=nsum2+1
        endif
    enddo
    gmn1=gsum1/nsum1
    gmn2=gsum2/nsum2

    gsum1=0.
    gsum2=0.
    do i=1,n
        if (gtype(i) == 1) then 
            gsum1=gsum1+(pr(i)-gmn1)**2
        else
            gsum2=gsum2+(pr(i)-gmn2)**2
        endif
    enddo
    gsum1=gsum1/nsum1
    gsum2=gsum2/nsum2
    call write_to('diffstress.summary', 'population 1:')
    call write_to('diffstress.summary', '   mean:', gmn1, '   std.dev.:', sqrt(gsum1))
    call write_to('diffstress.summary', 'population 2:')
    call write_to('diffstress.summary', '   mean:', gmn2, '   std.dev.:', sqrt(gsum2))
    call write_to('diffstress.summary', 'ratio: ', gmn1/gmn2)
    call close_file('diffstress.summary')

    do i=1,n
        ftmax(i)=0.
    enddo

    nonzeroc=0
    do k=1,contactknt
        if (contfn(k) /= 0.) then
            ft=sqrt(contfn(k)*contfn(k)+contft(k)*contft(k))
            write(42,*) nonzeroc+1,ft
            ftmax(contacti(k))=max(ft,ftmax(contacti(k)))
            ftmax(contactj(k))=max(ft,ftmax(contactj(k)))
            ftsum=ftsum+ft
            nonzeroc=nonzeroc+1
        endif
    enddo
    do i=1,n
        write(41,*) 2*radius(i),ftmax(i)
    enddo
    ftsum=ftsum/dble(nonzeroc)
    do i=1,n
        if (ftmax(i) >= 2.*ftsum) then
            if (gtype(i) == 1) then
                g1=g1+1
            else
                g2=g2+1
            endif
        endif
    enddo

    do k=1,contactknt
        ksmall=0
        kbig=0
        if (contfn(k) /= 0.) then
            if (ft < 0.5*ftsum) ksmall=ksmall+1
            if (ft > 2.*ftsum) kbig=kbig+1
            if (gtype(contacti(k))*gtype(contactj(k))==4) ibnet=ibnet+1
        endif
    enddo

    !open(unit=48,file='forcechain.summary')
    do i=1,n
        if (gtype(i) == 1) then
            c1=c1+coordn(i)
        else
            c2=c2+coordn(i)
        endif
    enddo
    c1=c1/dble(nsum1)
    c2=c2/dble(nsum2)
     
    call write_to("forcechain.summary", "coordination of  pop. 1 and pop.2:", c1, c2)      
    call write_to("forcechain.summary", "average contact force:", ftsum) 
    prop1=dble(g1)/dble(nsum1)
    prop2=dble(g2)/dble(nsum2)
    call write_to("forcechain.summary", "proportion in chains from pop. 1, pop. 2:", prop1, prop2)

    call write_to("forcechain.summary", "small, big and total contact forces:", ksmall, kbig, nonzeroc)
    call write_to("forcechain.summary", "small/total  big/total:", dble(ksmall)/dble(nonzeroc), dble(kbig)/dble(nonzeroc))

    call close_file("forcechain.summary")

    np1=0
    np2=0
    nps1=0
    nps2=0
    diffsum1=0.
    diffsum2=0.

    do i=1,n
        if (gtype(i) == 1) then
            np1=np1+1
            diffsum1=diffsum1+sign(i)
            if (sign(i) > 1.) nps1=nps1+1
        else
            np2=np2+1
            diffsum2=diffsum2+sign(i)
            if (sign(i) > 1.) nps2=nps2+1
        endif
    enddo
    diffsum1=diffsum1/dble(np1)
    diffsum2=diffsum2/dble(np2)
    diffrac1=dble(nps1)/dble(np1)
    diffrac2=dble(nps2)/dble(np2)
    bnet=dble(ibnet)/dble(np2)

    close(48)
    
    call write_to("force.summary", afrac1,ftsum,c1,c2,prop1,prop2,bnet)
    call close_file("force.summary")
    
    call write_to("stress.summary", afrac1, diffsum1, diffsum2, diffrac1, diffrac2)
    call close_file("stress.summary")

    do i=1,n
        if (maxpr > 0.) then
            p1=pr(i)/maxpr
        else
            p1=0.
        endif
        if (p1 > 1.) then
            red=1.
            blue=.6
            green=.6
        else
            index = max(0,nint(127*p1))
            red = rd(index+1,1)
            blue = bl(index+1,1)
            gs=blue-1.8*(1.-blue)
            green = max(0.0d0,gs)
        endif


        ikolor=1

        theta(i)=cron(i)*180./pi
        write(psunit,28)  radius(i), &
        rx(i),ry(i), &
        red,green,blue,p1,theta(i),' tagrain'

        !if (ib(3) == -1) then
        !      duplicate atoms wrapping around x direction
        if (rx(i) > xright_out-radius(i)) then
            write(psunit,28)  radius(i), &
            rx(i)-xl,ry(i), &
            red,green,blue,p1,theta(i),' tagrain'
        else if (rx(i) < xleft_out+radius(i)) then
            write(psunit,28)  radius(i), &
            rx(i)+xl,ry(i), &
            red,green,blue,p1,theta(i),' tagrain'
        endif
        !endif
    enddo
    
    icnt=0

    !           print*,' number of contacts = ',icnt

    call endpspage(psunit)
    close(psunit)

    if (igif) then
        gfile = ""
        gfile = ' trigrainxxx.gif'
        gfile(10:12) = pspage(7:9)
        com1='sh mkjpeg.sh '//trim(pspage)//trim(gfile)
        print*, com1
        call execute_command_line(com1)
    endif


end subroutine

end module
