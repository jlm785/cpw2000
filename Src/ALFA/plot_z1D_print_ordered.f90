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

!>  prints the list of atoms ordered by position along the third axis.
!>
!>  \author       Jose Luis Martins
!>  \version      5.06
!>  \date         2 March 2023.
!>  \copyright    GNU Public License v2

subroutine plot_z1D_print_ordered(ipr, io, ntype, natom, nameat, rat,    &
       ntot, iptype, ipnatom,                                            &
       mxdtyp, mxdatm)

! adapted from print_crystal, 2 March 2023. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of number of atoms of a given type

  integer, intent(in)                ::  ipr                             !<  should be equal to one if information is to be printed.
  integer, intent(in)                ::  io                              !<  default tape number
  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i
  character(len=2), intent(in)       ::  nameat(mxdtyp)                  !<  chemical symbol for the type i
  real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)            !<  lattice coordinates of atom j of type i

! output

  integer, intent(out)                ::  ntot                           !<  total number of atoms
  integer, intent(out)                ::  iptype(mxdatm*mxdtyp)          !<  points to atom type
  integer, intent(out)                ::  ipnatom(mxdatm*mxdtyp)         !<  points to atom of a given type

!  allocatable arrays

  real(REAL64), allocatable           ::  zcoord(:)                      !  z (third) lattice coordinate of aton in the primitive cell
  integer, allocatable                ::  indx(:)                        !  index by increasing zcoord

  integer, allocatable                ::  indnt(:), indj(:)              !  temporary indices

!  local variables


! counters

  integer       ::  nt, j, icount


  ntot = 0
  do nt = 1,ntype
    ntot = ntot + natom(nt)
  enddo

  allocate(zcoord(ntot))
  allocate(indx(ntot), indnt(ntot), indj(ntot))

  icount = 0
  do nt = 1,ntype
    do j = 1,natom(nt)
      icount = icount + 1
      zcoord(icount) = rat(3,j,nt) - floor(rat(3,j,nt)+0.000001)
      indnt(icount) = nt
      indj(icount) = j
    enddo
  enddo

  call sort(ntot, zcoord, indx)

  do icount = 1,ntot
    iptype(icount) = indnt(indx(icount))
    ipnatom(icount) = indj(indx(icount))
  enddo

  if(ipr == 1) then
    write(io,*)
    write(io,'(5x,"index",3x,"atom type",5x,"lattice coordinates")')
    write(io,*)
    do icount = 1,ntot
      write(6,'(5x,i5,3x,a2,1x,i5,3x,3f10.4)')                           &
               icount, nameat(iptype(icount)), ipnatom(icount),          &
               (rat(j,ipnatom(icount),iptype(icount)),j=1,2),            &
               zcoord(indx(icount))
    enddo
    write(io,*)
  endif

  deallocate(zcoord)
  deallocate(indx, indnt, indj)

  return
end subroutine plot_z1D_print_ordered
