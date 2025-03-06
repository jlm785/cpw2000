!>  Converts a string to upper case
!>
!>  \author       https://rosettacode.org/wiki/String_case#Fortran
!>  \version      5.11
!>  \date         27 February 2025.
!>  \copyright    GNU Free Document License 1.3


subroutine string_to_upper(str)

  character(*), intent(in out)      :: str                               !<  string to be converted
  integer :: i

  do i = 1, len(str)
    select case(str(i:i))
      case("a":"z")
        str(i:i) = achar(iachar(str(i:i))-32)
    end select
  end do

  return

end subroutine string_to_upper
