module mod_pssubs
    implicit none

    private
    
    public :: newpsdoc, endpspage

    contains

!=======================================================================
subroutine newpsdoc(psunit,xleft,xright,ybot,ytop,pspage,ibox,time)

    integer thispage
    integer psunit,ibox
    
    real(8) xleft,xright,ybot,ytop,xr,xl,sc,sch,scv,xshift,xmid
    real(8) yshift,ymid,yb,yt,time
    real(8) x1,x2,y1,y2,xd
    
    character(len=30) pspage
    
    10   format(a)
    20   format(e12.5,2x,e12.5,a)
    25   format(a,i4)
    28   format(a,4(f5.0,1x))
    40   format(f8.3,f8.3,a)
    45   format(i5,i5,a)
    50   format(a,f8.0,a)

    sch=550./(xright-xleft)/3
    scv=550./(ytop-ybot)
    sc=min(sch,scv)
    
    xmid=.5*(xleft+xright)
    ymid=.5*(ytop+ybot)

    yshift=275.-ymid

    if (ibox == 0) then
        thispage=0

        sc=sc*1.3
        xshift=100.-xleft*sc
        yshift=100.-ybot*sc
        open(newunit=psunit,file=trim(pspage))
        write(psunit,10) '%!PS'
        write(psunit,25) '%%Pages: 1'
        x1=xleft*sc+xshift
        y1=ybot*sc+yshift
        x2=xright*sc+xshift
        y2=ytop*sc+yshift
        xd=.15*(x2-x1)
        write(psunit,28) '%%BoundingBox: ',x1-xd,y1-xd,x2+xd,y2+xd

        write(psunit,10) '/grain {'
        write(psunit,10) ' /lw exch def'
        write(psunit,10) ' /bl exch def'
        write(psunit,10) ' /gr exch def'
        write(psunit,10) ' /rd exch def'
        write(psunit,10) ' /rad exch def'
        write(psunit,10) ' /y exch def'
        write(psunit,10) ' /x exch def'
        write(psunit,10) ' rd gr bl setrgbcolor'
        write(psunit,10) ' newpath x y rad 0 360 arc closepath fill'
        ! uncomment the next three lines for drawing a black line around grains
        write(psunit,10) ' lw setlinewidth'
        write(psunit,10) ' 0 0 0 setrgbcolor'
        write(psunit,10) ' newpath x y rad 0 360 arc closepath stroke'
        write(psunit,10) '} def'

        write(psunit,10) '/arrow {'
        write(psunit,10) ' /bl exch def'
        write(psunit,10) ' /gr exch def'
        write(psunit,10) ' /rd exch def'
        write(psunit,10) ' /r exch def'
        write(psunit,10) ' /a exch def'
        write(psunit,10) ' /y exch def'
        write(psunit,10) ' /x exch def'
        write(psunit,10) ' gsave'
        write(psunit,10) ' x y translate a a scale r rotate'
        write(psunit,10) ' rd gr bl vector'
        write(psunit,10) ' grestore'
        write(psunit,10) '} def'

        write(psunit,10) '/vector {'
        write(psunit,10) ' /blue exch def'
        write(psunit,10) ' /green exch def'
        write(psunit,10) ' /red exch def'
        write(psunit,10) ' red green blue setrgbcolor'  
        write(psunit,10) ' 0 0 1 0 360 arc fill' 
        write(psunit,10) ' 0 setlinewidth '
        write(psunit,10) ' 0 0 moveto 10 0 lineto stroke'
        write(psunit,10) '} def'

        write(psunit,10) '/vector2 {'
        write(psunit,10) ' 0 0 moveto 1 0 lineto stroke'
        write(psunit,10) '} def'

        write(psunit,10) '/cross'
        write(psunit,10) '{'
        write(psunit,10) ' /t1 exch def'
        write(psunit,10) ' /s1 exch def'
        write(psunit,10) ' /r exch def'
        write(psunit,10) ' /y1 exch def'
        write(psunit,10) ' /x1 exch def'

        write(psunit,10) ' /t180 {t1 180 add} def'
        write(psunit,10) ' 0 0 0 setrgbcolor'
        write(psunit,10) ' s1 setlinewidth'     
        write(psunit,10) ' gsave x1 y1 translate r r scale t1 rotate'
        write(psunit,10) ' vector2 grestore'
        write(psunit,10) ' gsave x1 y1 translate r r scale t180 rotate'
        write(psunit,10) ' vector2 grestore'
        write(psunit,10) '} def'


        write(psunit,10) '/tagrain'
        write(psunit,10) '{'
        write(psunit,10) ' /th exch def'
        write(psunit,10) ' /s1 exch def'
        write(psunit,10) ' /bl exch def'
        write(psunit,10) ' /gr exch def'
        write(psunit,10) ' /rd exch def'
        write(psunit,10) ' /y1 exch def'
        write(psunit,10) ' /x1 exch def'
        write(psunit,10) ' /rad exch def'
        write(psunit,10) ' x1 y1 rad rd gr bl grain'

        write(psunit,10) '} def'
        write(psunit,10) ' '

        write(psunit,10) '/tbgrain'
        write(psunit,10) '{'
        write(psunit,10) ' /bl exch def'
        write(psunit,10) ' /gr exch def'
        write(psunit,10) ' /rd exch def'
        write(psunit,10) ' /y1 exch def'
        write(psunit,10) ' /x1 exch def'
        write(psunit,10) ' /rad exch def'
        write(psunit,10) ' x1 y1 rad rd gr bl 0 grain'

        write(psunit,10) '} def'
        write(psunit,10) ' '
        write(psunit,10) '/tbgrain2'
        write(psunit,10) '{'
        write(psunit,10) ' /bl exch def'
        write(psunit,10) ' /gr exch def'
        write(psunit,10) ' /rd exch def'
        write(psunit,10) ' /y1 exch def'
        write(psunit,10) ' /x1 exch def'
        write(psunit,10) ' /rad exch def'
        write(psunit,10) ' x1 y1 rad rd gr bl 0.1 grain'

        write(psunit,10) '} def'
        write(psunit,10) ' '
        thispage=thispage+1
        write(psunit,25) '%%Page:',thispage
        write(psunit,10) '2 setlinewidth'
        write(psunit,45) nint(x1), nint(y1-.5*xd),' moveto'
        write(psunit,10) '/Helvetica findfont 15 scalefont setfont'
        write(psunit,50) '(T = ',time,') show'


        write(psunit,10) 'gsave'
        write(psunit,20) xshift,yshift,' translate '
        write(psunit,40) sc,sc,' scale'
        write(psunit,10) ' 0 setlinewidth'
        write(psunit,10) 'newpath'
        yb=ybot
        yt=ytop
        xl=xleft
        xr=xright
        write(psunit,20) xl,yb, ' moveto'
        write(psunit,20) xl,yt, ' lineto'
        write(psunit,20) xr,yt, ' lineto'
        write(psunit,20) xr,yb, ' lineto'
        write(psunit,10) 'closepath clip stroke'
        write(psunit,10) '1 1 1 setrgbcolor'
        write(psunit,20) xl,yb, ' moveto'
        write(psunit,20) xl,yt, ' lineto'
        write(psunit,20) xr,yt, ' lineto'
        write(psunit,20) xr,yb, ' lineto'
        write(psunit,10) 'closepath fill'


    else if (ibox == 1) then

        write(psunit,10) 'grestore'
        write(psunit,10) 'gsave'
        xshift=350.
        write(psunit,20) xshift+xmid,yshift,' translate '
        write(psunit,40) sc,sc,' scale'
        write(psunit,10) ' 0 setlinewidth'
        write(psunit,10) 'newpath'
        yb=ybot
        yt=ytop
        xl=xleft
        xr=xright
        write(psunit,20) xl,yb, ' moveto'
        write(psunit,20) xl,yt, ' lineto'
        write(psunit,20) xr,yt, ' lineto'
        write(psunit,20) xr,yb, ' lineto'
        write(psunit,10) 'closepath clip stroke'

    else if (ibox == 2) then

        write(psunit,10) 'grestore'
        write(psunit,10) 'gsave'
        xshift=250.+1.1*(xright-xleft)*sc
        write(psunit,20) xshift+1.2*xmid,yshift,' translate '
        write(psunit,40) sc,sc,' scale'
        write(psunit,10) ' 0 setlinewidth'
        write(psunit,10) 'newpath'
        yb=ybot
        yt=ytop
        xl=xleft
        xr=xright
        write(psunit,20) xl,yb, ' moveto'
        write(psunit,20) xl,yt, ' lineto'
        write(psunit,20) xr,yt, ' lineto'
        write(psunit,20) xr,yb, ' lineto'
        write(psunit,10) 'closepath clip stroke'

    endif

end subroutine
!=======================================================================
    !subroutine endpspage(psunit,thispage)
subroutine endpspage(psunit)

    !integer psunit,thispage
    integer psunit

    write(psunit,'(a)') 'showpage'
    write(psunit,'(a)') 'grestore'

end subroutine

end module
