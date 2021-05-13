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

!>  Identifies the nearest neighbors of each atom
!>  using a fuzzy Voronoi algorithm
!>
!>  \author       Jose Luis Martins
!>  \version      5.01
!>  \date         22 April 2021
!>  \copyright    GNU Public License v2

subroutine voronoi_neighb(adot, ntype, natom, nameat, rat,               &
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

! output

  integer, intent(out)               ::  nneighb(mxdatm,mxdtyp)          !<  number of neighbors for each atom
  integer, intent(out)               ::  neighbtype(mxdnb,mxdatm,mxdtyp) !<  type of nth neighbors for each atom  return
  real(REAL64), intent(out)          ::  rneighb(3,mxdnb,mxdatm,mxdtyp)  !<  relative position of the n-th neighbor
  real(REAL64), intent(out)          ::  wneighb(mxdnb,mxdatm,mxdtyp)    !<  strength of the n-th neighbor bond

! allocatable arrays

  real(REAL64), allocatable          ::  rtry(:,:,:)     !  cartesian coordinates of centered atoms
  real(REAL64), allocatable          ::  dist(:)         !  distance between atoms
  real(REAL64), allocatable          ::  rvec(:,:)       !  vector connecting atoms
  integer, allocatable               ::  indx(:)         !  ordering of closest atoms
  real(REAL64), allocatable          ::  radius(:)       !  conventional atomic radius
  integer, allocatable               ::  itype(:)        !  type of atom

  real(REAL64), allocatable          ::  rcar(:,:)
  real(REAL64), allocatable          ::  rbs(:)

  integer, allocatable               ::  idneib(:)
  real(REAL64), allocatable          ::  rneib(:)
  real(REAL64), allocatable          ::  wneib(:)

! local variables

  integer           ::  natot                   !  total number of atoms
  real(REAL64)      ::  avec(3,3), bvec(3,3)
  integer           ::  nsize                   !  size of temporary array

  integer           ::  nneib

! parameters

  real(REAL64), parameter    ::  UM = 1.0_REAL64
  real(REAL64), parameter    ::  ZERO = 0.0_REAL64
  integer, parameter         ::  NSHELL = 2               !  number of explored cells in each direction (1 should be enough...)

! counters

  integer        ::  nt, ja
  integer        ::  nt2, ja2
  integer        ::  n, j
  integer        ::  icount
  integer        ::  i1, i2, i3

  call adot_to_avec_sym(adot,avec,bvec)

  allocate(rtry(3,mxdatm,mxdtyp))
  allocate(radius(mxdtyp))

  natot = 0
  do nt = 1,ntype
    natot = natot + natom(nt)
    call p_tbl_radii(nameat(nt),radius(nt))
    do ja = 1,natom(nt)
      do n = 1,3
        rtry(n,ja,nt) = ZERO
        do j = 1,3
          rtry(n,ja,nt) = rtry(n,ja,nt) + avec(n,j)*                     &
                             (rat(j,ja,nt)-nint(rat(j,ja,nt)))
        enddo
      enddo
    enddo
  enddo

  nsize = natot*(2*NSHELL+1)*(2*NSHELL+1)*(2*NSHELL+1)
  allocate(dist(nsize))
  allocate(rvec(3,nsize))
  allocate(indx(nsize))
  allocate(itype(nsize))

  allocate(rcar(3,MXDNB))
  allocate(rbs(MXDNB))

  allocate(idneib(MXDNB))
  allocate(rneib(MXDNB))
  allocate(wneib(MXDNB))

! loop over reference atom

  do nt = 1,ntype
  do ja = 1,natom(nt)
    icount = 0

!   loop over other atoms and cells

    do nt2 = 1,ntype
    do ja2 = 1,natom(nt2)
      do i1 = -NSHELL,NSHELL
      do i2 = -NSHELL,NSHELL
      do i3 = -NSHELL,NSHELL

        icount = icount + 1
        do j = 1,3
          rvec(j,icount) = rtry(j,ja2,nt2) - rtry(j,ja,nt)               &
                    - i1*avec(j,1) - i2*avec(j,2) - i3*avec(j,3)
        enddo
        dist(icount) = rvec(1,icount)*rvec(1,icount) +                   &
                       rvec(2,icount)*rvec(2,icount) +                   &
                       rvec(3,icount)*rvec(3,icount)
        itype(icount) = nt2
      enddo
      enddo
      enddo

    enddo
    enddo

!   sorts the neighbors and keeps the first MXDNB

    call sort(icount,dist,indx)

    do j = 1,MXDNB
      rcar(1,j) = rvec(1,indx(j))
      rcar(2,j) = rvec(2,indx(j))
      rcar(3,j) = rvec(3,indx(j))
      rbs(j) = radius(itype(indx(j)))
    enddo

!   applies the Voronoi fuzzy algorithm

    call voronoi_neighb_one(MXDNB, MXDNB, rcar, rbs,                     &
       nneib, idneib, rneib, wneib)

!   recovers information

    nneighb(ja,nt) = nneib
    do j = 1,nneib
      neighbtype(j,ja,nt) = itype(indx(idneib(j)))
      rneighb(1,j,ja,nt) = rcar(1,idneib(j))
      rneighb(2,j,ja,nt) = rcar(2,idneib(j))
      rneighb(3,j,ja,nt) = rcar(3,idneib(j))
      wneighb(j,ja,nt) = wneib(j)*(rneib(1)/rneib(j))**2.5
    enddo

  enddo
  enddo

  deallocate(rtry)
  deallocate(radius)

  deallocate(dist)
  deallocate(rvec)
  deallocate(indx)
  deallocate(itype)

  deallocate(rcar)
  deallocate(rbs)

  deallocate(idneib)
  deallocate(rneib)
  deallocate(wneib)

  return
end subroutine voronoi_neighb
