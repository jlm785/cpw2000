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

!>  prints information about the nearest neighbors of each atom
!>  obtained with a fuzzy Voronoi algorithm
!>
!>  \author       Jose Luis Martins
!>  \version      5.01
!>  \date         22 April 2021
!>  \copyright    GNU Public License v2

subroutine voronoi_neighb_print(adot, ntype, natom, nameat, rat,         &
      nneighb, neighbtype, rneighb, wneighb,                             &
      mxdtyp, mxdatm, mxdnb)

! copyright INESC-MN/Jose Luis Martins

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of number of atoms of a given type
  integer, intent(in)                ::  mxdnb                           !<  array dimension of maximum number neighbors

  integer, intent(in)                ::  ntype                           !< number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !< number of atoms of type i
  character(len=2),intent(in)        ::  nameat(mxdtyp)                  !<  chemical symbol of type of atom
  real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)            !< lattice coordinates of atom j of type i

  real(REAL64), intent(in)           ::  adot(3,3)                       !< metric in real space

  integer, intent(in)                ::  nneighb(mxdatm,mxdtyp)          !<  number of neighbors for each atom
  integer, intent(in)                ::  neighbtype(mxdnb,mxdatm,mxdtyp) !<  type of nth neighbors for each atom  return
  real(REAL64), intent(in)           ::  rneighb(3,mxdnb,mxdatm,mxdtyp)  !<  relative position of the n-th neighbor
  real(REAL64), intent(in)           ::  wneighb(mxdnb,mxdatm,mxdtyp)    !<  strength of the n-th neighbor bond

! local variables

  integer           ::  io                   !  default output
  real(REAL64)      ::  avec(3,3), bvec(3,3)
  real(REAL64)      ::  rcar(3)
  real(REAL64)      ::  dist
  integer           ::  nclose

! counters

  integer        ::  nt, ja
  integer        ::  n, j


  io = 6

  call adot_to_avec_sym(adot,avec,bvec)

  write(io,*)
  write(io,*)
  write(io,*) '    neighbour analysis of the structure'
  write(io,*) '    uses a fuzzy Voronoi polyhedra algorithm (Becke)'
  write(io,*)

  do nt = 1,ntype
    write(io,*)
    write(io,'(4x,"There are ",i5," atoms of ",a2)') natom(nt), nameat(nt)
    write(io,*)
    do ja = 1,natom(nt)
      do j = 1,3
        rcar(j) = avec(j,1)*rat(1,ja,nt) + avec(j,2)*rat(2,ja,nt) +      &
                  avec(j,3)*rat(3,ja,nt)
      enddo
      nclose = 0
      do n = 1,nneighb(ja,nt)
        if(wneighb(n,ja,nt) > 0.5) nclose = nclose + 1
      enddo
      write(io,*)
      write(io,'(8x,a2," atom at ",3f10.3," has ",i2," + ",i2,           &
             & " neighbors")') nameat(nt),(rcar(j),j=1,3),               &
                         nclose, nneighb(ja,nt)-nclose
      write(io,*)
      write(io,'(8x,"Element",4x,"Distance",4x,"Bond strength",20x,"Direction")')
      write(io,*)
      do n = 1,nneighb(ja,nt)
        dist = rneighb(1,n,ja,nt)*rneighb(1,n,ja,nt) +                   &
               rneighb(2,n,ja,nt)*rneighb(2,n,ja,nt) +                   &
               rneighb(3,n,ja,nt)*rneighb(3,n,ja,nt)
        dist = sqrt(dist)
        if(n == nclose +1) write(io,*)
        write(io,'(12x,a2,2x,f10.3,4x,f10.3,11x,3f10.3)')                &
             nameat(neighbtype(n,ja,nt)), dist, wneighb(n,ja,nt),        &
                     (rneighb(j,n,ja,nt),j=1,3)
      enddo
      write(io,*)
    enddo
    write(io,*)
  enddo


  return
end subroutine voronoi_neighb_print
