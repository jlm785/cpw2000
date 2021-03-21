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

!>  Reads the data in the band structure file  
!>  and calculates required array sizes

  subroutine dos_read_size(filename, io, nrk, mxdbnd, nx)

! Written November 16, 2013. jlm
! Modified 12 June 2014
! Modified, documentation, 19 September 2020. JLM
! Modified, new file format, 13 december 2020. JLM
! copyright  J.L.Martins, INESC-MN.

! version 4.99 of cpw

  implicit none

! input

  integer, intent(in)                ::  io                              !<  tape numbers
  character(len=*), intent(in)       ::  filename                        !<  input file

! output

  integer, intent(out)               ::  nrk                             !<  size of k-points in irreducible wedge of Z
  integer, intent(out)               ::  mxdbnd                          !<  size of number of bands
  integer, intent(out)               ::  nx(3)                           !<  number of k-points in each direction in regular grid

! local allocatable variable

  integer, allocatable  ::  nband(:)

! local variables

  integer       ::  ios

! counters

  integer                            ::  irk

  open(unit=io, file=trim(filename), status='OLD', form='UNFORMATTED', iostat=ios)

  if(ios /= 0) then
    write(6,*) '  Stopped in dos_read_size: ', trim(filename), ' not found'

    stop

  endif

  read(io) nrk

  allocate(nband(nrk))

  read(io) nband(1:nrk)
  read(io) nx(1),nx(2),nx(3)

  mxdbnd = 1
  do irk = 1,nrk
    if(nband(irk) > mxdbnd) mxdbnd = nband(irk)
  enddo

  close(unit=io)

  deallocate(nband)

  return
  end subroutine dos_read_size
