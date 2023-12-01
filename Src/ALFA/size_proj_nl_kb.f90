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

!>  calculates the number of Kleinman-Bylander projectors
!>
!>  \author       Jose Luis Martins
!>  \version      5.09
!>  \date         December 20, 2013. 30 November 2023.
!>  \copyright    GNU Public License v2

subroutine size_proj_nl_kb(ntype, natom, nkb, nanl, nanlso, nanlspin,    &
              mxdtyp)

! written December 20, 2013. jlm
! adapted from proj_nl_kb.
! modified February 10, 2014. jlm
! Modified, documentation, January 2020. JLM
! Added nanlspin, 30 November 2023.   NEW API   NEW API

  implicit none

!  integer, parameter          :: REAL64 = selected_real_kind(12)
 integer, parameter          ::  LMAX = 3

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms

  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i

  integer, intent(in)                ::  nkb(0:LMAX,-1:1,mxdtyp)         !<  KB pseudo.  normalization for atom k, ang. mom. l

! output

  integer, intent(out)               ::  nanl                            !<  number of Kleinman-Bylander projectors without spin-orbit.
  integer, intent(out)               ::  nanlso                          !<  number of Kleinman-Bylander projectors with spin-orbit.
  integer, intent(out)               ::  nanlspin                        !<  number of Kleinman-Bylander projectors with or without spin-orbit depending on pseudopotential

! counters

  integer    ::  k, l, ic

  nanl = 0
  do k = 1,ntype
    do l = 0,LMAX
      if(nkb(l,0,k) /= 0) nanl = nanl + (2*l+1)*natom(k)
    enddo
  enddo

  nanlso = 0
  do k=1,ntype
    if(nkb(0, 1,k) /= 0) nanlso = nanlso + 2*natom(k)
    do l = 1,LMAX
      if(nkb(l,-1,k) /= 0) nanlso = nanlso + (2*l)*natom(k)
      if(nkb(l, 1,k) /= 0) nanlso = nanlso + (2*l+2)*natom(k)
    enddo
  enddo

  nanlspin = 0
  do k=1,ntype

    ic = 0
    if(nkb(0, 1,k) /= 0) ic = ic+1
    do l = 1,LMAX
      if(nkb(l,-1,k) /= 0) ic = ic+1
      if(nkb(l, 1,k) /= 0) ic = ic+1
    enddo

    if(ic == 0) then
      do l = 0,LMAX
        if(nkb(l,0,k) /= 0) nanlspin = nanlspin + (2*l+1)*natom(k)
      enddo
    else
      if(nkb(0, 1,k) /= 0) nanlspin = nanlspin + 2*natom(k)
      do l = 1,LMAX
        if(nkb(l,-1,k) /= 0) nanlspin = nanlspin + (2*l)*natom(k)
        if(nkb(l, 1,k) /= 0) nanlspin = nanlspin + (2*l+2)*natom(k)
      enddo
    endif
  enddo

  return

end subroutine size_proj_nl_kb
