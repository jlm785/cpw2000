!>  Module with the handle of the Intel MKL library DFTI descriptor.
!>
!>  \author       José Luís Martins
!>  \version      5.13
!>  \date         16 March 2010.
!>  \copyright    GNU Public License v2

module mkl_handle

  use mkl_dfti

  type(DFTI_DESCRIPTOR), POINTER :: hand

end module mkl_handle
