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

!>  Calculates derivatives with respect to wave-vector k.
!>  Hamiltonian acting on wave-function, d H / d k, energy levels d E_i / d k,
!>  wave-functions d | psi_i > / d k.
!>  The last two may not be accurate if degeneracies are present,
!>  or input wave-function is not accurate.
!>
!>  \author       Jose Luis Martins
!>  \version      5.06
!>  \date         18 January 2023.
!>  \copyright    GNU Public License v2

subroutine berry_derivative(rkpt, mtxd, neig, isort, ekpg, lpsi,         &
    psi, ei, dhdkpsi, dpsidk, deidk, bcurv, bmag,                        &
    ng, kgv,                                                             &
    vscr, kmscr,                                                         &
    nqnl, delqnl, vkb, nkb,                                              &
    ntype, natom, rat, adot,                                             &
    mxdtyp, mxdatm, mxdlqp, mxddim, mxdbnd, mxdgve, mxdscr)


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)
  integer, parameter          :: LMAX = 3                                !   hard coded max ang. mom.

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of number of atoms of a given type
  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands
  integer, intent(in)                ::  mxdlqp                          !<  array dimension for local potential
  integer, intent(in)                ::  mxdgve                          !<  array dimension of G-space vectors
  integer, intent(in)                ::  mxdscr                          !<  array dimension for screening potential

  real(REAL64), intent(in)           ::  rkpt(3)                         !<  k-point reciprocal lattice coordinates
  integer, intent(in)                ::  mtxd                            !<  wavefunction dimension
  integer, intent(in)                ::  neig                            !<  wavefunction dimension

  logical, intent(in)                ::  lpsi                            !<  also calculates the wave-function derivative and Berry curvature

  integer, intent(in)                ::  isort(mxddim)                   !<  G-vector corresponding to coefficient i of wavefunction
  real(REAL64), intent(in)           ::  ekpg(mxddim)                    !<  kinetic energy (hartree) of k+g-vector of row/column i

  complex(REAL64), intent(in)        ::  psi(mxddim, mxdbnd)             !<  |psi> (in principle eigen-functions)
  real(REAL64), intent(in)           ::  ei(mxddim)                      !<  eigenvalue no. i. (hartree)

  integer, intent(in)                ::  ng                              !<  size of g-space
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  G-vectors in reciprocal lattice coordinates

  integer, intent(in)                ::  kmscr(7)                        !<  max value of kgv(i,n) used for
  real(REAL64), intent(in)           ::  vscr(mxdscr)                    !<  screened

  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i
  real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)            !<  lattice coordinates of atom j of type i
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in real space
  integer, intent(in)                ::  nqnl(mxdtyp)                    !<  number of points for pseudo interpolation for atom k
  real(REAL64), intent(in)           ::  delqnl(mxdtyp)                  !<  step used in the pseudo interpolation for atom k
  real(REAL64), intent(in)        ::  vkb(-2:mxdlqp,0:LMAX,-1:1,mxdtyp)  !<  (1/q**l) * KB nonlocal pseudo. for atom k, ang. mom. l. NOT normalized to vcell, hartree
  integer, intent(in)                ::  nkb(0:LMAX,-1:1,mxdtyp)         !<  KB pseudo.  normalization for atom k, ang. mom. l

! output

  complex(REAL64), intent(out)       ::  dhdkpsi(mxddim, mxdbnd, 3)      !<  d H d k | psi>  (lattice coordinates)
  complex(REAL64), intent(out)       ::  dpsidk(mxddim, mxdbnd, 3)       !<  d |psi> / d k  (lattice coordinates)
  real(REAL64), intent(out)          ::  deidk(mxdbnd, 3)                !<  d  E / d k  (lattice coordinates)

  real(REAL64), intent(out)          ::  bcurv(3,mxdbnd)                 !<  Berry curvature  (lattice coordinates)
  real(REAL64), intent(out)          ::  bmag(3,mxdbnd)                  !<  Orbital magnetization  (lattice coordinates)

! local allocatable variables

  real(REAL64), allocatable          ::  xnlkb(:)
  complex(REAL64), allocatable       ::  anlga(:,:)

  complex(REAL64), allocatable       ::  danlgadrk(:,:,:)                !  d anlga / d rkpt
  complex(REAL64), allocatable       ::  d2anlgadrk2(:,:,:,:)            !  d 2 anlga / d 2 rkp (unused)

  complex(REAL64), allocatable       ::  hdpsidk(:,:,:)                  !  (H - E_n) ( d | psi_n > / d k )

! local variables

  integer           ::  nder                                             !  order of derivative

  integer           ::  nanl, nanlso
  integer           ::  mxdanl

! parameters

  real(REAL64), parameter       ::  ZERO = 0.0_REAL64
  real(REAL64), parameter       ::  TOL = 1.0E-7_REAL64
  real(REAL64), parameter       ::  PI = 3.14159265358979323846_REAL64

! counters

  integer    ::  j, n

! external

  complex(REAL64) ,external   :: zdotc



! projector and derivatives

  call size_proj_nl_kb(ntype, natom, nkb, nanl, nanlso, mxdtyp)

  mxdanl = nanl

  allocate(anlga(mxddim,mxdanl))
  allocate(xnlkb(mxdanl))

  allocate(danlgadrk(mxddim,mxdanl,3))
  allocate(d2anlgadrk2(1,1,3,3))

  nder = 1

  call proj_nl_kb_der_c16(rkpt, mtxd, isort, nanl, nder,                 &
      ng, kgv,                                                           &
      nqnl, delqnl, vkb, nkb,                                            &
      ntype, natom, rat, adot,                                           &
      anlga, xnlkb, danlgadrk, d2anlgadrk2,                              &
      mxdtyp, mxdatm, mxdlqp, mxddim, mxdanl, mxdgve)

  call berry_dhdk_psi(rkpt, adot, mtxd, neig, psi, dhdkpsi       ,       &
      kgv, isort,                                                        &
      nanl, anlga, xnlkb, danlgadrk,                                     &
      mxddim, mxdbnd, mxdgve, mxdanl)

  do n = 1,neig
    do j = 1,3
      deidk(n,j) = real( zdotc(mtxd,psi(:,n),1,dhdkpsi(:,n,j),1), REAL64)
    enddo
  enddo

  if(lpsi) then

    call berry_stern_solve(mtxd, neig, psi, ei, dhdkpsi, dpsidk, TOL,    &
       isort, ekpg,                                                      &
       vscr, kmscr,                                                      &
       ng, kgv,                                                          &
       nanl, anlga, xnlkb,                                               &
       mxddim, mxdbnd, mxdgve, mxdscr, mxdanl)

    do n = 1,neig

      bcurv(1,n) = -2*dimag(zdotc(mtxd, dpsidk(:,n,2), 1, dpsidk(:,n,3), 1))
      bcurv(2,n) = -2*dimag(zdotc(mtxd, dpsidk(:,n,3), 1, dpsidk(:,n,1), 1))
      bcurv(3,n) = -2*dimag(zdotc(mtxd, dpsidk(:,n,1), 1, dpsidk(:,n,2), 1))

    enddo

    allocate(hdpsidk(mxddim,mxdbnd,3))

    do j = 1,3

      call hk_psi_c16(mtxd, neig, dpsidk(:,:,j), hdpsidk(:,:,j), .TRUE., &
             ng, kgv,                                                    &
             ekpg, isort, vscr, kmscr,                                   &
             anlga, xnlkb, nanl,                                         &
             mxddim, mxdbnd, mxdanl, mxdgve, mxdscr)

      do n = 1,neig

        call zaxpy(mtxd, cmplx(-ei(n),ZERO,REAL64), dpsidk(:,n,j), 1, hdpsidk(:,n,j), 1)

      enddo

    enddo

    do n = 1,neig

      bmag(1,n) = dimag(zdotc(mtxd, dpsidk(:,n,2), 1, hdpsidk(:,n,3), 1)) / 2
      bmag(2,n) = dimag(zdotc(mtxd, dpsidk(:,n,3), 1, hdpsidk(:,n,1), 1)) / 2
      bmag(3,n) = dimag(zdotc(mtxd, dpsidk(:,n,1), 1, hdpsidk(:,n,2), 1)) / 2

    enddo

    deallocate(hdpsidk)

  endif

  deallocate(anlga)
  deallocate(xnlkb)

  deallocate(danlgadrk)
  deallocate(d2anlgadrk2)

  return

end subroutine berry_derivative
