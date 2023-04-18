module mod_postScriptTool
!***********************************************************************
!
! MODULE FOR CSH SCRIPT WRITING
!
! Martin Svoboda
! svobod.martin@gmail.com .or. svobodam@icpf.cas.cz 
! 2019
!
!***********************************************************************
    implicit none
    
    private
    
    public :: postScriptTool_init, postScriptTool_clean
    
    contains
!=======================================================================
subroutine postScriptTool_init(resolution)

   integer :: resolution
   
   call create_mkJpeg(resolution)
   call create_mkJpeg_all(resolution)
   call create_pstojpeg_dws_sh()
   
   call execute_command_line("")
    
end subroutine 
!=======================================================================
subroutine postScriptTool_clean()
   
   call execute_command_line("rm mkjpeg.sh mkjpeg_all.sh pstojpeg.dws.sh")
    
end subroutine          
!=======================================================================
subroutine create_mkJpeg(resolution)

    integer :: resolution
    integer :: nu
    
    open(newunit=nu, file="mkjpeg.sh")
    write(nu,'(a)') "#!/bin/bash"
    write(nu,'(a,i0,a)') "csh -f ./pstojpeg.dws.sh -portrait -xsize ", resolution, " ${1%.*}"
    close(nu)
    
end subroutine    
!=======================================================================
subroutine create_mkJpeg_all(resolution)

    integer :: resolution
    integer :: nu
    
    open(newunit=nu, file="mkjpeg_all.sh")
    write(nu,'(a)') "#!/bin/bash"
    write(nu,'(a)') "for f in bpage*.ps; do"
    write(nu,'(a,i0,a)') "    csh -f ./pstojpeg.dws.sh -portrait -xsize ", resolution, " ${f%.*}"
    write(nu,'(a)') "done"
    close(nu)

    
end subroutine    
!=======================================================================
subroutine create_pstojpeg_dws_sh()

    integer :: nu
    
    open(newunit=nu, file="pstojpeg.dws.sh")
    write(nu,'(a)') '#!/bin/csh -f'
    write(nu,'(a)') '#'
    write(nu,'(a)') '#    Uses ghostscript to translate an Encapsulated PostScript file to'
    write(nu,'(a)') '#    Portable Anymap format file(s).'
    write(nu,'(a)') '#    pstopnm will create as many files as the number of pages in'
    write(nu,'(a)') '#    the Postscript document.  The name of the files will be'
    write(nu,'(a)') '#    psfile001.ppm, psfile002.ppm, etc.'
    write(nu,'(a)') '#    The ouput files will contain the area inside the BoundingBox.'
    write(nu,'(a)') '#    If BoundingBox parameters are not found in the PostScript'
    write(nu,'(a)') '#    document, default values are used.'
    write(nu,'(a)') '#'
    write(nu,'(a)') '#'
    write(nu,'(a)') '#    Usage: pstopnm [-forceplain] [-help] [-llx s] [-lly s]'
    write(nu,'(a)') '#               [-urx s] [-ury s] [-nocrop] [-pbm|-pgm|-ppm]'
    write(nu,'(a)') '#               [-verbose] [-xborder n] [-xmax n] [-xsize n]'
    write(nu,'(a)') '#               [-yborder n] [-ymax n] [-ysize n]'
    write(nu,'(a)') '#               [-portrait] [-landscape] psfile[.ps]'
    write(nu,'(a)') '#'
    write(nu,'(a)') '#    Copyright (C) 1992 by Alberto Accomazzi, Smithsonian Astrophysical'
    write(nu,'(a)') '#    Observatory (alberto@cfa.harvard.edu).'
    write(nu,'(a)') '#'
    write(nu,'(a)') '#    Permission to use, copy, modify, and distribute this software and its'
    write(nu,'(a)') '#    documentation for any purpose and without fee is hereby granted,'
    write(nu,'(a)') '#    provided that the above copyright notice appear in all copies and'
    write(nu,'(a)') '#    that both that copyright notice and this permission notice appear'
    write(nu,'(a)') '#    in supporting documentation.  This software is provided "as is"'
    write(nu,'(a)') '#    without express or implied warranty.'
    write(nu,'(a)') '#'
    write(nu,'(a)') 'set noglob'
    write(nu,*)
    write(nu,'(a)') 'set progname = $0'
    write(nu,'(a)') 'set progname = $progname:t'
    write(nu,'(a)') 'set filtertail = ""'
    write(nu,'(a)') 'set filterhead = "jpeg"'
    write(nu,'(a)') 'set xsize = 0'
    write(nu,'(a)') 'set ysize = 0'
    write(nu,'(a)') 'set xres = ""'
    write(nu,'(a)') 'set yres = ""'
    write(nu,*)
    write(nu,'(a)') '# default values: max image x and y sizes'
    write(nu,'(a)') 'set xmax = 612'
    write(nu,'(a)') 'set ymax = 792'
    write(nu,'(a)') '# default values: image area fits in a 8.5x11 sheet with 1 inch border'
    write(nu,'(a)') 'set llx = 72'
    write(nu,'(a)') 'set lly = 72'
    write(nu,'(a)') 'set urx = 540'
    write(nu,'(a)') 'set ury = 720'
    write(nu,'(a)') '# default values: x and y borders are 10% of x and y size'
    write(nu,'(a)') 'set xborder = "0.1"'
    write(nu,'(a)') 'set yborder = "0.1"'
    write(nu,'(a)') '# default values: orientation is unknown'
    write(nu,'(a)') 'set orient = 0'
    write(nu,*)
    write(nu,'(a)') 'set psfile = ""'
    write(nu,'(a)') 'set USAGE = "Usage: $progname [-forceplain] [-help] [-llx s] [-lly s]\'
    write(nu,'(a)') '[-urx s] [-ury s] [-landscape] [-portrait]\'
    write(nu,'(a)') '[-nocrop] [-pbm|-pgm|-ppm] [-verbose] [-xborder s] [-xmax s]\'
    write(nu,'(a)') '[-xsize s] [-yborder s] [-ymax s] [-ysize s] psfile[.ps]"'
    write(nu,'(a)') "alias usage 'echo $USAGE; exit 1'"
    write(nu,*)
    write(nu,'(a)') 'while ($#argv > 0)'
    write(nu,'(a)') '    switch ($argv[1])'
    write(nu,'(a)') '    case -h*:    # -help'
    write(nu,'(a)') '    usage'
    write(nu,'(a)') '    breaksw'
    write(nu,'(a)') '    case -pbm:'
    write(nu,'(a)') '    case -pgm:'
    write(nu,'(a)') '    case -ppm:'
    write(nu,'(a)') '    set filterhead = `echo "$argv[1]" | sed "s/-//1"`'
    write(nu,'(a)') '    breaksw'
    write(nu,'(a)') '    case -llx:'
    write(nu,'(a)') '    shift argv'
    write(nu,'(a)') '    if ($#argv == 0) eval usage'
    write(nu,'(a)') '    set llx = `(echo "scale=4";echo "$argv[1] * 72")|bc -l`'
    write(nu,'(a)') '    set nobb'
    write(nu,'(a)') '    breaksw'
    write(nu,'(a)') '    case -lly:'
    write(nu,'(a)') '    shift argv'
    write(nu,'(a)') '    if ($#argv == 0) eval usage'
    write(nu,'(a)') '    set lly = `(echo "scale=4";echo "$argv[1] * 72")|bc -l`'
    write(nu,'(a)') '    set nobb'
    write(nu,'(a)') '    breaksw'
    write(nu,'(a)') '    case -urx:'
    write(nu,'(a)') '    shift argv'
    write(nu,'(a)') '    if ($#argv == 0) eval usage'
    write(nu,'(a)') '    set urx = `(echo "scale=4";echo "$argv[1] * 72")|bc -l`'
    write(nu,'(a)') '    set nobb'
    write(nu,'(a)') '    breaksw'
    write(nu,'(a)') '    case -ury:'
    write(nu,'(a)') '    shift argv'
    write(nu,'(a)') '    if ($#argv == 0) eval usage'
    write(nu,'(a)') '    set ury = `(echo "scale=4";echo "$argv[1] * 72")|bc -l`'
    write(nu,'(a)') '    set nobb'
    write(nu,'(a)') '    breaksw'
    write(nu,'(a)') '    case -no*:    # -nocrop'
    write(nu,'(a)') '    set nocrop'
    write(nu,'(a)') '    breaksw'
    write(nu,'(a)') '    case -xs*:    # -xsize'
    write(nu,'(a)') '    shift argv'
    write(nu,'(a)') '    if ($#argv == 0) eval usage'
    write(nu,'(a)') '    @ xsize = $argv[1]'
    write(nu,'(a)') '    breaksw'
    write(nu,'(a)') '    case -ys*:    # -ysize'
    write(nu,'(a)') '    shift argv'
    write(nu,'(a)') '    if ($#argv == 0) eval usage'
    write(nu,'(a)') '    @ ysize = $argv[1]'
    write(nu,'(a)') '    breaksw'
    write(nu,'(a)') '    case -xm*:    # -xmax'
    write(nu,'(a)') '    shift argv'
    write(nu,'(a)') '    if ($#argv == 0) eval usage'
    write(nu,'(a)') '    @ xmax = $argv[1]'
    write(nu,'(a)') '    breaksw'
    write(nu,'(a)') '    case -ym*:    # -ymax'
    write(nu,'(a)') '    shift argv'
    write(nu,'(a)') '    if ($#argv == 0) eval usage'
    write(nu,'(a)') '    @ ymax = $argv[1]'
    write(nu,'(a)') '    breaksw'
    write(nu,'(a)') '    case -xb*:    # -xborder'
    write(nu,'(a)') '    shift argv'
    write(nu,'(a)') '    if ($#argv == 0) eval usage'
    write(nu,'(a)') '    set xborder = $argv[1]'
    write(nu,'(a)') '    breaksw'
    write(nu,'(a)') '    case -yb*:    # -yborder'
    write(nu,'(a)') '    shift argv'
    write(nu,'(a)') '    if ($#argv == 0) eval usage'
    write(nu,'(a)') '    set yborder = $argv[1]'
    write(nu,'(a)') '    breaksw'
    write(nu,'(a)') '    case -f*:    # -forceplain'
    write(nu,'(a)') '    set filtertail = ""'
    write(nu,'(a)') '    breaksw'
    write(nu,'(a)') '    case -v*:    # -verbose'
    write(nu,'(a)') '    set verb'
    write(nu,'(a)') '    breaksw'
    write(nu,'(a)') '    case -po*:    # -portrait'
    write(nu,'(a)') '    set orient = 1'
    write(nu,'(a)') '    breaksw'
    write(nu,'(a)') '    case -la*:    # -landscape'
    write(nu,'(a)') '    set orient = 2'
    write(nu,'(a)') '    breaksw'
    write(nu,'(a)') '    case -*:'
    write(nu,'(a)') '    echo "${progname}: Unknown option $argv[1]"'
    write(nu,'(a)') '    usage'
    write(nu,'(a)') '    breaksw'
    write(nu,'(a)') '    default:    # input file'
    write(nu,'(a)') '    set psfile = $argv[1]'
    write(nu,'(a)') '    set psfilePS = ${argv[1]}.ps'
    write(nu,'(a)') '    set ppmfile = `basename $argv[1] .ps`'
    write(nu,'(a)') '    breaksw'
    write(nu,'(a)') '    endsw'
    write(nu,'(a)') '    shift argv'
    write(nu,'(a)') 'end'
    write(nu,*)
    write(nu,'(a)') 'if ($psfilePS =~ "") eval usage'
    write(nu,'(a)') 'if (! -f $psfilePS) then'
    write(nu,'(a)') '    echo "${progname}: file $psfilePS not found"'
    write(nu,'(a)') '    usage'
    write(nu,'(a)') 'endif'
    write(nu,*)
    write(nu,'(a)') 'set bb = `grep "%%BoundingBox" $psfilePS`'
    write(nu,'(a)') 'if ($?nobb == 0 && $#bb == 5) then'
    write(nu,'(a)') '    set llx = $bb[2]'
    write(nu,'(a)') '    set lly = $bb[3]'
    write(nu,'(a)') '    set urx = $bb[4]'
    write(nu,'(a)') '    set ury = $bb[5]'
    write(nu,'(a)') 'else'
    write(nu,'(a)') '    if ($?nobb == 0) \'
    write(nu,'(a)') '    echo "${progname}: warning: BoundingBox not found in input file"'
    write(nu,'(a)') 'endif'
    write(nu,*)
    write(nu,'(a)') 'set tmpsx = `(echo "scale=4";echo "$urx - $llx")|bc -l`'
    write(nu,'(a)') 'set tmpsy = `(echo "scale=4";echo "$ury - $lly")|bc -l`'
    write(nu,*)
    write(nu,'(a)') '# see if orientation was specified'
    write(nu,'(a)') 'if ($orient == 0) then'
    write(nu,'(a)') '    # no orientation was specified; compute default orientation'
    write(nu,'(a)') '    set tmpx = 0'
    write(nu,'(a)') '    set tmpy = 0'
    write(nu,'(a)') '    set tmpsx1 = $tmpsx:r'
    write(nu,'(a)') '    set tmpsy1 = $tmpsy:r'
    write(nu,'(a)') '    # default is landscape mode'
    write(nu,'(a)') '    set orient = 2'
    write(nu,'(a)') '    if ($xsize == 0 && $ysize == 0) then'
    write(nu,'(a)') '    set tmpx = $xmax'
    write(nu,'(a)') '    set tmpy = $ymax'
    write(nu,'(a)') '    else'
    write(nu,'(a)') '    if ($xsize != 0) set tmpx = $xsize'
    write(nu,'(a)') '    if ($ysize != 0) set tmpy = $ysize'
    write(nu,'(a)') '    endif'
    write(nu,'(a)') '    if ($tmpx == 0 || $tmpy == 0) then'
    write(nu,'(a)') '    # only one size was specified'
    write(nu,'(a)') '    if ($tmpsy1 > $tmpsx1) set orient = 1'
    write(nu,'(a)') '    else'
    write(nu,'(a)') '    # no size or both sizes were specified'
    write(nu,'(a)') '    if ($tmpsy1 > $tmpsx1 && $tmpy > $tmpx) set orient = 1'
    write(nu,'(a)') '    if ($tmpsx1 > $tmpsy1 && $tmpx > $tmpy) set orient = 1'
    write(nu,'(a)') '    endif'
    write(nu,'(a)') 'endif'
    write(nu,*)
    write(nu,'(a)') '# now reset BoundingBox llc and total size to take into account margin'
    write(nu,'(a)') 'set llx = `(echo "scale=4";echo "$llx - $tmpsx * $xborder")|bc -l`'
    write(nu,'(a)') 'set lly = `(echo "scale=4";echo "$lly - $tmpsy * $yborder")|bc -l`'
    write(nu,'(a)') 'set urx = `(echo "scale=4";echo "$urx + $tmpsx * $xborder")|bc -l`'
    write(nu,'(a)') 'set ury = `(echo "scale=4";echo "$ury + $tmpsy * $yborder")|bc -l`'
    write(nu,'(a)') '# compute image area size'
    write(nu,'(a)') 'set sx = `(echo "scale=4";echo "$tmpsx + 2 * $xborder * $tmpsx")|bc -l`'
    write(nu,'(a)') 'set sy = `(echo "scale=4";echo "$tmpsy + 2 * $yborder * $tmpsy")|bc -l`'
    write(nu,*)
    write(nu,'(a)') 'if ($orient != 1) then'
    write(nu,'(a)') '    # render image in landscape mode'
    write(nu,'(a)') '    set tmpsx = $sx'
    write(nu,'(a)') '    set sx = $sy'
    write(nu,'(a)') '    set sy = $tmpsx'
    write(nu,'(a)') 'endif'
    write(nu,*)
    write(nu,'(a)') '# if xsize or ysize was specified, compute resolution from them'
    write(nu,'(a)') 'if ($xsize != 0) set xres = `(echo "scale=4";echo "$xsize *72 / $sx")|bc -l`'
    write(nu,'(a)') 'if ($ysize != 0) set yres = `(echo "scale=4";echo "$ysize *72 / $sy")|bc -l`'
    write(nu,*)
    write(nu,'(a)') 'if ($xres =~ "" && $yres !~ "") then'
    write(nu,'(a)') '    # ysize was specified, xsize was not; compute xsize based on ysize'
    write(nu,'(a)') '    set xres = $yres'
    write(nu,'(a)') '    set xsize = `(echo "scale=4";echo "$sx * $xres /72 + 0.5")|bc -l`'
    write(nu,'(a)') '    set xsize = $xsize:r'
    write(nu,'(a)') 'else'
    write(nu,'(a)') '    if ($yres =~ "" && $xres !~ "") then'
    write(nu,'(a)') '    # xsize was specified, ysize was not; compute ysize based on xsize'
    write(nu,'(a)') '    set yres = $xres'
    write(nu,'(a)') '    set ysize = `(echo "scale=4";echo "$sy * $yres /72 + 0.5")|bc -l`'
    write(nu,'(a)') '    set ysize = $ysize:r'
    write(nu,'(a)') '    else'
    write(nu,'(a)') '    if ($xres =~ "" && $yres =~ "") then'
    write(nu,'(a)') '        # neither xsize nor ysize was specified; compute them from'
    write(nu,'(a)') '        # xmax and ymax'
    write(nu,'(a)') '        set xres = `(echo "scale=4";echo "$xmax *72/$sx")|bc -l`'
    write(nu,'(a)') '        set yres = `(echo "scale=4";echo "$ymax *72/$sy")|bc -l`'
    write(nu,'(a)') '        set xres = `(echo "scale=4";echo "if($xres>$yres)$yres";echo "if($yres>$xres)$xres")|bc -l`'
    write(nu,'(a)') '        set yres = $xres'
    write(nu,'(a)') '        if ($?nocrop) then'
    write(nu,'(a)') '        # keep output file dimensions equal to xmax and ymax'
    write(nu,'(a)') '        set xsize = $xmax'
    write(nu,'(a)') '        set ysize = $ymax'
    write(nu,'(a)') '        else'
    write(nu,'(a)') '        set xsize = `(echo "scale=4";echo "$sx * $xres /72+0.5")|bc -l`'
    write(nu,'(a)') '        set ysize = `(echo "scale=4";echo "$sy * $yres /72+0.5")|bc -l`'
    write(nu,'(a)') '        endif'
    write(nu,'(a)') '        set xsize = $xsize:r'
    write(nu,'(a)') '        set ysize = $ysize:r'
    write(nu,'(a)') '    endif'
    write(nu,'(a)') '    endif'
    write(nu,'(a)') 'endif'
    write(nu,*)
    write(nu,'(a)') '# translate + rotate image, if necessary'
    write(nu,'(a)') 'if ($orient == 1) then'
    write(nu,'(a)') '    # portrait mode'
    write(nu,'(a)') '    # adjust offsets'
    write(nu,'(a)') '    set llx = `(echo "scale=4";echo "$llx - ($xsize *72/$xres - $sx)/2")|bc -l`'
    write(nu,'(a)') '    set lly = `(echo "scale=4";echo "$lly - ($ysize *72/$yres - $sy)/2")|bc -l`'
    write(nu,'(a)') '    set pstrans = "$llx neg $lly neg translate"'
    write(nu,'(a)') 'else'
    write(nu,'(a)') '    # landscape mode'
    write(nu,'(a)') '    # adjust offsets'
    write(nu,'(a)') '    set llx = `(echo "scale=4";echo "$llx - ($ysize *72/$yres - $sy)/2")|bc -l`'
    write(nu,'(a)') '    set ury = `(echo "scale=4";echo "$ury + ($xsize *72/$xres - $sx)/2")|bc -l`'
    write(nu,'(a)') '    set pstrans = "90 rotate $llx neg $ury neg translate"'
    write(nu,'(a)') 'endif'
    write(nu,*)
    write(nu,*)
    write(nu,'(a)') 'if ($?verb) then'
    write(nu,'(a)') '    echo "sx = $sx"'
    write(nu,'(a)') '    echo "sy = $sy"'
    write(nu,'(a)') '    echo "xres    = $xres"'
    write(nu,'(a)') '    echo "yres    = $yres"'
    write(nu,'(a)') '    echo "xsize = $xsize"'
    write(nu,'(a)') '    echo "ysize = $ysize"'
    write(nu,'(a)') '    echo -n "orientation "'
    write(nu,'(a)') '    if ($orient == 1) then'
    write(nu,'(a)') '    echo "portrait"'
    write(nu,'(a)') '    else'
    write(nu,'(a)') '    echo "landscape"'
    write(nu,'(a)') '    endif'
    write(nu,'(a)') '    echo "PS header: $pstrans"'
    write(nu,'(a)') 'endif'
    write(nu,*)
    write(nu,'(a)') 'echo "${progname}: writing $filterhead file(s)"'
    write(nu,'(a)') 'echo "filterhead: ${filterhead} filtertail: ${filtertail}"'
    write(nu,*)
    write(nu,'(a)') 'echo $pstrans | \'
    write(nu,'(a)') '    gs -sDEVICE=${filterhead}${filtertail} \'
    write(nu,'(a)') '       -sOutputFile=$psfile.jpeg \'
    write(nu,'(a)') '       -g${xsize}x${ysize} \'
    write(nu,'(a)') '       -r${xres}x${yres} \'
    write(nu,'(a)') '       -q - $psfilePS \'
    write(nu,*)
    write(nu,*)
    write(nu,*)
    write(nu,*)

    close(nu)

end subroutine    

end module
