module mod_outputWriter
use iso_fortran_env,only:i4 => int32,&
                         r4 => real32, r8 => real64
!***********************************************************************
!
! MODULE FOR OUTPUT WRITING
!
! Martin Svoboda
! svobod.martin@gmail.com .or. svobodam@icpf.cas.cz 
! 2019
!
!***********************************************************************

    implicit none
   
    private
    
    interface write_to
        procedure :: output_real, output_2real, output_3real, output_4real, &
                     output_char, output_char_int, output_char_real, output_char_2real, &
                     output_char_3real, output_char_real_char_real, output_5real, &
                     output_6real, output_7real, output_char_3int, output_real_int, &
                     output_5real_3int, output_5real_3int_2real, &
                     output_2int_3real, output_int_real_int, output_10real
    end interface
    
    public :: write_to, close_file
    
    

contains
!=======================================================================
integer function open_file(fname,styleIn) result(num)
    character(*), intent(in) :: fname
    character(*), optional, intent(in) :: styleIn
    integer(kind=i4) :: ios
    character(len=50) :: style=""

    if(present(styleIn)) then
    style = styleIn
    else
    style = "unknown"
    endif

    inquire(file=fname,number=num)
    if(num/=-1) return
    open(newunit=num, file=fname, iostat=ios, status=style)
  
end function
!=======================================================================
subroutine close_file(fname)
    character(*), intent(in) :: fname
    integer(kind=i4) :: ios, num

    inquire(file=fname,number=num)
    
    if(num/=-1) then
        close(num, iostat=ios)
    else
        print *, "close_file error: there is not this file"
        stop
    endif
  
end subroutine
!=======================================================================
subroutine output_char(fname,var)

    character(*),intent(in) :: fname
    character(*),intent(in) :: var
    integer(kind=i4) :: num

    num = open_file(trim(fname))

    write(num,'(*(g0,1x))') var

end subroutine
!=======================================================================
subroutine output_real_int(fname,var1,var2)

    character(*),intent(in) :: fname
    real(kind=r8),intent(in) :: var1
    integer(kind=i4),intent(in) :: var2
    integer(kind=i4) :: num

    num = open_file(trim(fname))

    write(num,'(*(g0,1x))') var1, var2

end subroutine
!=======================================================================
subroutine output_char_int(fname,var1,var2)

    character(*),intent(in) :: fname
    character(*),intent(in) :: var1
    integer(kind=i4),intent(in) :: var2
    integer(kind=i4) :: num

    num = open_file(trim(fname))

    write(num,'(*(g0,1x))') var1, var2

end subroutine
!=======================================================================
subroutine output_char_3int(fname,var1,var2,var3,var4)

    character(*),intent(in) :: fname
    character(*),intent(in) :: var1
    integer(kind=i4),intent(in) :: var2,var3,var4
    integer(kind=i4) :: num

    num = open_file(trim(fname))

    write(num,'(*(g0,1x))') var1, var2, var3, var4

end subroutine
!=======================================================================
subroutine output_char_real_char_real(fname,var1,var2,var3,var4)

    character(*),intent(in) :: fname
    character(*),intent(in) :: var1,var3
    real(kind=r8),intent(in) :: var2,var4
    integer(kind=i4) :: num

    num = open_file(trim(fname))

    write(num,'(*(g0,1x))') var1, var2, var3, var4

end subroutine
!=======================================================================
subroutine output_char_2real(fname,var1,var2,var3)

    character(*),intent(in) :: fname
    character(*),intent(in) :: var1
    real(kind=r8),intent(in) :: var2,var3
    integer(kind=i4) :: num

    num = open_file(trim(fname))

    write(num,'(*(g0,1x))') var1, var2, var3

end subroutine
!=======================================================================
subroutine output_char_3real(fname,var1,var2,var3,var4)

    character(*),intent(in) :: fname
    character(*),intent(in) :: var1
    real(kind=r8),intent(in) :: var2,var3,var4
    integer(kind=i4) :: num

    num = open_file(trim(fname))

    write(num,'(*(g0,1x))') var1, var2, var3, var4

end subroutine
!=======================================================================
subroutine output_char_real(fname,var1,var2)

    character(*),intent(in) :: fname
    character(*),intent(in) :: var1
    real(kind=r8),intent(in) :: var2
    integer(kind=i4) :: num

    num = open_file(trim(fname))

    write(num,'(*(g0,1x))') var1, var2

end subroutine
!=======================================================================
subroutine output_real(fname,var)

    character(*),intent(in) :: fname
    real(kind=r8),intent(in) :: var
    integer(kind=i4) :: num

    num = open_file(trim(fname))

    write(num,'(*(g0,1x))') var

end subroutine
!=======================================================================
subroutine output_2real(fname,var1,var2)

    character(*),intent(in) :: fname
    real(kind=r8),intent(in) :: var1, var2
    integer(kind=i4) :: num

    num = open_file(trim(fname))

    write(num,'(*(g0,1x))') var1, var2

end subroutine
!=======================================================================
subroutine output_int_real_int(fname,var1,var2,var3)

    character(*),intent(in) :: fname
    integer(kind=i4),intent(in) :: var1, var3
    real(kind=r8),intent(in) :: var2
    integer(kind=i4) :: num

    num = open_file(trim(fname))

    write(num,'(*(g0,1x))') var1, var2, var3

end subroutine
!=======================================================================
subroutine output_3real(fname,var1,var2,var3)

    character(*),intent(in) :: fname
    real(kind=r8),intent(in) :: var1, var2, var3
    integer(kind=i4) :: num

    num = open_file(trim(fname))

    write(num,'(*(g0,1x))') var1, var2, var3

end subroutine
!=======================================================================
subroutine output_4real(fname,var1,var2,var3,var4)

    character(*),intent(in) :: fname
    real(kind=r8),intent(in) :: var1, var2, var3, var4
    integer(kind=i4) :: num

    num = open_file(trim(fname))

    write(num,'(*(g0,1x))') var1, var2, var3, var4

end subroutine
!=======================================================================
subroutine output_5real(fname,var1,var2,var3,var4,var5)

    character(*),intent(in) :: fname
    real(kind=r8),intent(in) :: var1, var2, var3, var4, var5
    integer(kind=i4) :: num

    num = open_file(trim(fname))

    write(num,'(*(g0,1x))') var1, var2, var3, var4, var5

end subroutine
!=======================================================================
subroutine output_2int_3real(fname,var1,var2,var3,var4,var5)

    character(*),intent(in) :: fname
    integer(kind=i4),intent(in) :: var1, var2
    real(kind=r8),intent(in) :: var3, var4, var5
    integer(kind=i4) :: num

    num = open_file(trim(fname))

    write(num,'(*(g0,1x))') var1, var2, var3, var4, var5

end subroutine
!=======================================================================
subroutine output_6real(fname,var1,var2,var3,var4,var5,var6)

    character(*),intent(in) :: fname
    real(kind=r8),intent(in) :: var1, var2, var3, var4, var5, var6
    integer(kind=i4) :: num

    num = open_file(trim(fname))

    write(num,'(*(g0,1x))') var1, var2, var3, var4, var5, var6

end subroutine
!=======================================================================
subroutine output_7real(fname,var1,var2,var3,var4,var5,var6,var7)

    character(*),intent(in) :: fname
    real(kind=r8),intent(in) :: var1, var2, var3, var4, var5, var6, var7
    integer(kind=i4) :: num

    num = open_file(trim(fname))

    write(num,'(*(g0,1x))') var1, var2, var3, var4, var5, var6, var7

end subroutine
!=======================================================================
subroutine output_5real_3int(fname,var1,var2,var3,var4,var5,var6,var7,var8)

    character(*),intent(in) :: fname
    real(kind=r8),intent(in) :: var1, var2, var3, var4, var5
    integer(kind=i4),intent(in) :: var6, var7, var8
    integer(kind=i4) :: num

    num = open_file(trim(fname))

    write(num,'(*(g0,1x))') var1, var2, var3, var4, var5, var6, var7, var8

end subroutine
!=======================================================================
subroutine output_5real_3int_2real(fname,var1,var2,var3,var4,var5,var6,var7,var8,var9,var10)

    character(*),intent(in) :: fname
    real(kind=r8),intent(in) :: var1, var2, var3, var4, var5, var9, var10
    integer(kind=i4),intent(in) :: var6, var7, var8
    integer(kind=i4) :: num

    num = open_file(trim(fname))

    write(num,'(*(g0,1x))') var1, var2, var3, var4, var5, var6, var7, var8, var9, var10

end subroutine
!=======================================================================
subroutine output_10real(fname,var1,var2,var3,var4,var5,var6,var7,var8,var9,var10)

    character(*),intent(in) :: fname
    real(kind=r8),intent(in) :: var1, var2, var3, var4, var5, var6, var7, var8, var9, var10
    integer(kind=i4) :: num

    num = open_file(trim(fname))

    write(num,'(*(g0,1x))') var1, var2, var3, var4, var5, var6, var7, var8, var9, var10

end subroutine

end module
