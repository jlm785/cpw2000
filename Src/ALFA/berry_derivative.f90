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
!>  It also calculates the tensors needed for Bery curvature,
!>  orbital magnetization, inverse effective mass, and the quantum metric.
!>
!>  The results may not be accurate if
!>  the input wave-function is not accurate.
!>
!>  \author       Jose Luis Martins
!>  \version      5.11
!>  \date         18 January 2023. 11 April 2024.
!>  \copyright    GNU Public License v2

subroutine berry_derivative(rkpt, mtxd, neig, isort, ekpg, lpsi,         &
    nlevel, levdeg, leveigs,                                             &
    psi, ei,                                                             &
    dhdkpsi, dpsidk, psidhdkpsi, tfqg, tgammamf, td2hdk2,                &
    ng, kgv,                                                             &
    vscr, kmscr,                                                         &
    nqnl, delqnl, vkb, nkb,                                              &
    ntype, natom, rat, adot,                                             &
    mxdtyp, mxdatm, mxdlqp, mxddim, mxdbnd, mxdgve, mxdscr,              &
    mxdlev, mxddeg)

! Written by Jose Luis Martins, 18 january 2023.
! Final debugged version 3 November 2023. JLM
! Modified, nanlspin, 30 November 2023. JLM
! Corrected psidhdkpsi_sp. 5 March 2024. JLM
! Modified, dimensions (t)bcurv, (t)qmetric. 4 April 2024. JLM
! Modified, complex tensors, quantities calculated elsewhere. 11 April 2024. JLM



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
  integer, intent(in)                ::  mxdlev                          !<  array dimension for number of levels
  integer, intent(in)                ::  mxddeg                          !<  array dimension for number of levels

  real(REAL64), intent(in)           ::  rkpt(3)                         !<  k-point reciprocal lattice coordinates
  integer, intent(in)                ::  mtxd                            !<  wavefunction dimension
  integer, intent(in)                ::  neig                            !<  number of bands

  logical, intent(in)                ::  lpsi                            !<  also calculates the wave-function derivative and Berry curvature

  integer, intent(in)                ::  isort(mxddim)                   !<  G-vector corresponding to coefficient i of wavefunction
  real(REAL64), intent(in)           ::  ekpg(mxddim)                    !<  kinetic energy (hartree) of k+g-vector of row/column i

  integer, intent(in)                ::  nlevel                          !<  number of energy levels
  integer, intent(in)                ::  levdeg(mxdlev)                  !<  degeneragy of level
  integer, intent(in)                ::  leveigs(mxdlev,mxddeg)          !<  points to degenerate level

  complex(REAL64), intent(in)        ::  psi(mxddim, mxdbnd)             !<  |psi> (in principle eigen-functions)
  real(REAL64), intent(in)           ::  ei(mxdbnd)                      !<  eigenvalue no. i. (hartree)

  integer, intent(in)                ::  ng                              !<  size of g-space
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  G-vectors in reciprocal lattice coordinates

  integer, intent(in)                ::  kmscr(7)                        !<  max value of kgv(i,n) used for
  real(REAL64), intent(in)           ::  vscr(mxdscr)                    !<  screened local potential

  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i
  real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)            !<  lattice coordinates of atom j of type i
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in real space

  integer, intent(in)                ::  nqnl(mxdtyp)                    !<  number of points for pseudo interpolation for atom k
  real(REAL64), intent(in)           ::  delqnl(mxdtyp)                  !<  step used in the pseudo interpolation for atom k
  real(REAL64), intent(in)        ::  vkb(-2:mxdlqp,0:LMAX,-1:1,mxdtyp)  !<  (1/q**l) * KB nonlocal pseudo. for atom k, ang. mom. l. NOT normalized to vcell, hartree
  integer, intent(in)                ::  nkb(0:LMAX,-1:1,mxdtyp)         !<  KB pseudo.  normalization for atom k, ang. mom. l

! output

  complex(REAL64), intent(out)       ::  dhdkpsi(mxddim, mxdbnd, 3)      !<  d H / d k | psi>  (lattice coordinates)
  complex(REAL64), intent(out)       ::  dpsidk(mxddim, mxdbnd, 3)       !<  d |psi> / d k  (lattice coordinates)

  complex(REAL64), intent(out)       ::  psidhdkpsi(mxddeg,mxddeg,3,mxdlev)  !<  <psi_n| d H / d k |psi_m>  for each energy level (lattice coordinates)

  complex(REAL64), intent(out)       ::  tfqg(3,3,mxddeg,mxddeg,mxdlev)      !<  Tensor F quantum geomtric, curvature and metric, (lattice coordinates)
  complex(REAL64), intent(out)       ::  tgammamf(3,3,mxddeg,mxddeg,mxdlev)  !<  Tensor Gamma-F, orbital magnetization and contribution to effective mass (lattice coordinates)
  complex(REAL64), intent(out)       ::  td2hdk2(3,3,mxddeg,mxddeg,mxdlev)   !<  <psi| d^2 H / d k^2 |psi> (lattice coordinates)

! local allocatable arrays

  real(REAL64), allocatable          ::  xnlkb(:)                        !  KB normalization without spin-orbit
  complex(REAL64), allocatable       ::  anlga(:,:)                      !  KB projectors without spin-orbit

  complex(REAL64), allocatable       ::  danlgadrk(:,:,:)                !  d anlga / d rkpt
  complex(REAL64), allocatable       ::  d2anlgadrk2(:,:,:,:)            !  d 2 anlga / d 2 rkp (unused)

  complex(REAL64), allocatable       ::  hmedpsidk(:,:,:)                !  (H - E_n) ( d | psi_n > / d k )

  complex(REAL64), allocatable       ::  vnl0(:,:)                       !  <Psi|V_NL|Psi>
  complex(REAL64), allocatable       ::  dvnl0drk(:,:,:)                 !  d <Psi|V_NL|Psi> d k
  complex(REAL64), allocatable       ::  d2vnl0drk2(:,:,:,:)             !  d^2 <Psi|V_NL|Psi> d k^2


! local variables

  integer           ::  nder                                             !  order of derivative

  integer           ::  nanl, nanlso, nanlspin                           !  number of projectors
  integer           ::  mxdanl                                           !  maximum number of projectors
  real(REAL64)      ::  vcell, bdot(3,3)                                 !  cell volume, metric reciprocal space

! parameters

  real(REAL64), parameter       ::  ZERO = 0.0_REAL64
  real(REAL64), parameter       ::  TOL = 1.0E-9_REAL64

! counters

  integer    ::  i, j, n, m
  integer    ::  nl, nk, mk

! external

  complex(REAL64) ,external   :: zdotc


!  allocate(dmat(mxddeg,mxddeg,3,3))

! projector and derivatives

  call size_proj_nl_kb(ntype, natom, nkb, nanl, nanlso, nanlspin,        &
      mxdtyp)

  mxdanl = nanl

  allocate(anlga(mxddim,mxdanl))
  allocate(xnlkb(mxdanl))

  allocate(danlgadrk(mxddim,mxdanl,3))
  allocate(d2anlgadrk2(mxddim,mxdanl,3,3))

  nder = 1
  if(lpsi) then
    nder = 2
  else
    nder = 1
  endif

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

  do nl = 1,nlevel
    do nk = 1,levdeg(nl)
      n = leveigs(nl,nk)
      do mk = 1,levdeg(nl)
        m = leveigs(nl,mk)
        do j = 1,3
          psidhdkpsi(nk,mk,j,nl) = zdotc(mtxd,psi(:,n),1,dhdkpsi(:,m,j),1)
        enddo
      enddo
    enddo
  enddo

  if(lpsi) then

!   solves Sternheimer equation

    call berry_stern_solve(mtxd, neig, psi, ei, dhdkpsi, dpsidk, TOL,    &
       nlevel, levdeg, leveigs,                                          &
       isort, ekpg,                                                      &
       vscr, kmscr,                                                      &
       ng, kgv,                                                          &
       nanl, anlga, xnlkb,                                               &
       mxddim, mxdbnd, mxdgve, mxdscr, mxdanl, mxdlev, mxddeg)

!   geometric quantities

    do nl = 1,nlevel
      do nk = 1,levdeg(nl)
        n = leveigs(nl,nk)
        do mk = 1,levdeg(nl)
          m = leveigs(nl,mk)

          do i = 1,3
          do j = 1,3
            tfqg(i,j,nk,mk,nl) = zdotc(mtxd, dpsidk(:,n,i), 1, dpsidk(:,m,j), 1)
          enddo
          enddo

        enddo
      enddo

    enddo

!   geometric*H quantities

    allocate(hmedpsidk(mxddim,mxdbnd,3))

    do j = 1,3

      call hk_psi_c16(mtxd, neig, dpsidk(:,:,j), hmedpsidk(:,:,j), .TRUE.,  &
             ng, kgv,                                                       &
             ekpg, isort, vscr, kmscr,                                      &
             anlga, xnlkb, nanl,                                            &
             mxddim, mxdbnd, mxdanl, mxdgve, mxdscr)

      do n = 1,neig

        call zaxpy(mtxd, cmplx(-ei(n),ZERO,REAL64), dpsidk(:,n,j), 1, hmedpsidk(:,n,j), 1)

      enddo

    enddo

    allocate(vnl0(mxdbnd,mxdbnd))
    allocate(dvnl0drk(mxdbnd,mxdbnd,3))
    allocate(d2vnl0drk2(mxdbnd,mxdbnd,3,3))

!   extra calculations may be needed later

    call psi_vnl_psi_der(mtxd, neig, nanl, psi, nder,                    &
        vnl0, dvnl0drk, d2vnl0drk2,                                      &
        anlga, xnlkb, danlgadrk, d2anlgadrk2,                            &
        mxddim, mxdbnd, mxdanl)

    call adot_to_bdot(adot, vcell, bdot)

    do nl = 1,nlevel
      do nk = 1,levdeg(nl)
        n = leveigs(nl,nk)
        do mk = 1,levdeg(nl)
          m = leveigs(nl,mk)

          do i = 1,3
          do j = 1,3
            tgammamf(i,j,nk,mk,nl) = zdotc(mtxd, dpsidk(:,n,i), 1, hmedpsidk(:,m,j), 1)
          enddo
          enddo

          do i = 1,3
          do j = 1,3
            td2hdk2(i,j,nk,mk,nl) = d2vnl0drk2(n,m,i,j)
            if(nk == mk) then
              td2hdk2(i,j,nk,mk,nl) = td2hdk2(i,j,nk,mk,nl) + bdot(i,j)
            endif
          enddo
          enddo

        enddo
      enddo

    enddo

    deallocate(vnl0)
    deallocate(dvnl0drk)
    deallocate(d2vnl0drk2)

    deallocate(hmedpsidk)

  endif

  deallocate(anlga)
  deallocate(xnlkb)

  deallocate(danlgadrk)
  deallocate(d2anlgadrk2)

  return

end subroutine berry_derivative
