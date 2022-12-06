!------------------------------------------------------------!
! This file is distributed as part of the cpw2000 code and   !
! under the terms of the GNU General Public License. See the !
! file `LICENSE' in the root directory of the cpw2000        !
! distribution, or http://www.gnu.org/copyleft/gpl.txt       !
!                                                            !
! The webpage of the cpw2000 code is not yet written         !
!                                                            !
! The cpw2000 code is hosted on GitHub:                      !
!                                                            !
! https://github.com/jlm785/cpw2000                          !
!------------------------------------------------------------!

!>  Reads the file "filename" (default PW_RHO_V.DAT) and returns the
!>  array dimensions for the crystal geometry
!>
!>  \author       Jose Luis Martins
!>  \version      5.06
!>  \date         6 December 2022.
!>  \copyright    GNU Public License v2

subroutine pw_rho_v_in_size_xtal(filename, io, mxdtyp, mxdatm)

! Extracted from pw_rho_v_in_size

  implicit none

!  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  character(len=*), intent(in)       ::  filename                        !<  name of file
  integer, intent(in)                ::  io                              !<  number of tape associated with filename.

! output
  integer, intent(out)               ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(out)               ::  mxdatm                          !<  array dimension of types of atoms

! local allocatable array

  integer,allocatable                ::  natom(:)                        !  number of atoms of type i

! local variables

  integer      ::  ioerr

! counters

  integer    ::  i


  open(unit=io,file=trim(filename),status='old',form='UNFORMATTED')

  read(io,IOSTAT=ioerr) mxdtyp


  allocate(natom(mxdtyp))

  read(io) (natom(i),i=1,mxdtyp)

  mxdatm = natom(1)
  do i = 1,mxdtyp
    if(mxdatm < natom(i)) mxdatm = natom(i)
  enddo

  deallocate(natom)

  close(unit=io)

  return

end subroutine pw_rho_v_in_size_xtal


