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

!>  identifies the neighbours of an atom using a fuzzy Voronoy
!>  polyhedra algorithm (JChemPhys 88, 2547, 1988).
!>  Version for atom 1!!!!
!>
!>  \author       Jose Luis Martins
!>  \version      5.01
!>  \date         12 october 2006, modernized 20 April 2021
!>  \copyright    GNU Public License v2

subroutine voronoi_neighb_one(nneigh, maxnb, rcar, rbs,                  &
    nneib, idneib, rneib, wneib)

  implicit none

  integer, parameter  :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)           ::  nneigh                                !<  number of atoms
  integer, intent(in)           ::  maxnb                                 !<  maximum number of neighbors

  real(REAL64), intent(in)      ::  rcar(3,nneigh)                        !<  cartesian coordinate of atom i
  real(REAL64), intent(in)      ::  rbs(nneigh)                           !<  size of atom i (Bragg-Slater table for example)

! output

  integer, intent(out)          ::  nneib                                 !<  number of neighbours of atom 1
  integer, intent(out)          ::  idneib(maxnb)                         !<  j-th neighbour of atom 1
  real(REAL64), intent(out)     ::  rneib(maxnb)                          !<  distance between atom 1 and idneib(j)
  real(REAL64), intent(out)     ::  wneib(maxnb)                          !<  strength of bond between atom 1 and idneib(j)

! allocatable arrays

  real(REAL64), allocatable     ::  vpw(:)

! local variables

  integer            ::  kmax
  real(REAL64)       ::  d, xk, rk(3), dxk
  real(REAL64)       ::  vmax, vm1, v0, vp1

! parameters

  real(REAL64), parameter    ::  UM = 1.0_REAL64
  real(REAL64), parameter    ::  ZERO = 0.0_REAL64
  real(REAL64), parameter    ::  THRESH = 0.1_REAL64
  integer, parameter         ::  NK = 50

! counters

  integer     ::  j, k


  if(nneigh <= 1) then
    nneib = 0
  else

    allocate(vpw(nneigh))

    nneib = 0

    do j = 2,nneigh

      d = (rcar(1,1)-rcar(1,j))*(rcar(1,1)-rcar(1,j)) +              &
          (rcar(2,1)-rcar(2,j))*(rcar(2,1)-rcar(2,j)) +              &
          (rcar(3,1)-rcar(3,j))*(rcar(3,1)-rcar(3,j))

      d = sqrt(d)

      vmax = zero
      do k = 0,NK

        xk = 0.2 + 0.6*real(k)/real(NK)
        rk(1) = xk*rcar(1,1) + (UM-xk)*rcar(1,j)
        rk(2) = xk*rcar(2,1) + (UM-xk)*rcar(2,j)
        rk(3) = xk*rcar(3,1) + (UM-xk)*rcar(3,j)
        call voronoi_fuzzy(rk,nneigh,rcar,rbs,vpw)

        if(vpw(1)*vpw(j) > vmax) then
          vmax = vpw(1)*vpw(j)
          kmax = k
        endif

      enddo

      if(vmax > THRESH/4) then

        xk = 0.2 + 0.6*real(kmax-1)/real(NK)
        rk(1) = xk*rcar(1,1) + (UM-xk)*rcar(1,j)
        rk(2) = xk*rcar(2,1) + (UM-xk)*rcar(2,j)
        rk(3) = xk*rcar(3,1) + (UM-xk)*rcar(3,j)
        call voronoi_fuzzy(rk,nneigh,rcar,rbs,vpw)
        vm1 = vpw(1)*vpw(j)

        xk = 0.2 + 0.6*real(kmax)/real(NK)
        rk(1) = xk*rcar(1,1) + (UM-xk)*rcar(1,j)
        rk(2) = xk*rcar(2,1) + (UM-xk)*rcar(2,j)
        rk(3) = xk*rcar(3,1) + (UM-xk)*rcar(3,j)
        call voronoi_fuzzy(rk,nneigh,rcar,rbs,vpw)
        v0 = vpw(1)*vpw(j)

        xk = 0.2 + 0.6*real(kmax+1)/real(NK)
        rk(1) = xk*rcar(1,1) + (UM-xk)*rcar(1,j)
        rk(2) = xk*rcar(2,1) + (UM-xk)*rcar(2,j)
        rk(3) = xk*rcar(3,1) + (UM-xk)*rcar(3,j)
        call voronoi_fuzzy(rk,nneigh,rcar,rbs,vpw)
        vp1 = vpw(1)*vpw(j)

        dxk = (vp1 - vm1) / (2*(2*v0 - vp1 -vm1))
        xk = 0.2 + 0.6*(real(kmax)+dxk)/real(NK)
        rk(1) = xk*rcar(1,1) + (UM-xk)*rcar(1,j)
        rk(2) = xk*rcar(2,1) + (UM-xk)*rcar(2,j)
        rk(3) = xk*rcar(3,1) + (UM-xk)*rcar(3,j)
        call voronoi_fuzzy(rk,nneigh,rcar,rbs,vpw)
        vmax = 4*vpw(1)*vpw(j)

      endif

      if(vmax > THRESH) then
        nneib = nneib+1

        if(nneib > maxnb) then
          write(6,*) '  STOPPED in cluster_neighbor'
          write(6,*) '  increase the value of maxnb'

          stop

        endif

        idneib(nneib) = j
        rneib(nneib) = d
        wneib(nneib) = vmax
      endif
    enddo

    call heapsort_par(nneib, rneib, wneib, idneib)

    deallocate(vpw)

  endif

  return
end subroutine voronoi_neighb_one
