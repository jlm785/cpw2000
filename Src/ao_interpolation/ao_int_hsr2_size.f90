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

!>  Reads the version and array dimensions at the beginning of
!>  hsr2.dat.  Must be followed by ao_int_hsr2_read on the same unit.
!>
!>  \author       claude
!>  \version      5.13
!>  \date         23 October 2026.
!>  \copyright    GNU Public License v2

subroutine ao_int_hsr2_size(io, lok, ntype, mxdatm, mxdlao, mxdlqp)

! Written 23 October 2026. claude

  implicit none

! input

  integer, intent(in)                ::  io                              !<  tape number of hsr2.dat (open, unformatted, at the beginning)

! output

  logical, intent(out)               ::  lok                             !<  the file has the expected version
  integer, intent(out)               ::  ntype                           !<  number of types of atoms
  integer, intent(out)               ::  mxdatm                          !<  array dimension of number of atoms of a given type
  integer, intent(out)               ::  mxdlao                          !<  array dimension of orbital per atom type
  integer, intent(out)               ::  mxdlqp                          !<  array dimension of the radial functions

! local variables

  character(len=24)    ::  version
  integer              ::  ioerr

! parameters

  character(len=24), parameter  ::  HSR2VER = 'cpw2000 hsr2 version 1  '


  lok = .FALSE.
  ntype = 0
  mxdatm = 0
  mxdlao = 0
  mxdlqp = 0

  read(io,iostat=ioerr) version
  if(ioerr /= 0) return
  if(version /= HSR2VER) return

  read(io,iostat=ioerr) ntype, mxdatm, mxdlao, mxdlqp
  if(ioerr /= 0) return

  lok = .TRUE.

  return

end subroutine ao_int_hsr2_size
