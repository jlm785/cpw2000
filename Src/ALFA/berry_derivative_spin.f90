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
!>  wave-functions d | psi_sp_i > / d k.
!>  It also calculates the tensors needed for Bery curvature,
!>  orbital magnetization, inverse effective mass, and the quantum metric.
!>
!>  The results may not be accurate if
!>  input wave-function is not accurate.
!>
!>  spin-wave-function version
!>
!>  \author       Jose Luis Martins
!>  \version      5.11
!>  \date         15 December 2023. 11 April 2024.
!>  \copyright    GNU Public License v2

subroutine berry_derivative_spin(rkpt, mtxd, neig, isort, ekpg, lpsi,    &
    nlevel, levdeg, leveigs,                                             &
    psi_sp, ei,                                                          &
    dhdkpsi_sp, dpsidk_sp, psidhdkpsi_sp, tfqg, tgammamf, td2hdk2,       &
    ng, kgv,                                                             &
    vscr_sp, kmscr, nsp,                                                 &
    nqnl, delqnl, vkb, nkb,                                              &
    ntype, natom, rat, adot,                                             &
    mxdtyp, mxdatm, mxdlqp, mxddim, mxdbnd, mxdgve, mxdscr,              &
    mxdlev, mxddeg, mxdnsp)

! Adapted from the non-spin version 15 December 2023. JLM
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
  integer, intent(in)                ::  mxdnsp                          !<  array dimension for number of spin components (1,2,4)

  real(REAL64), intent(in)           ::  rkpt(3)                         !<  k-point reciprocal lattice coordinates
  integer, intent(in)                ::  mtxd                            !<  wavefunction dimension (not including spin)
  integer, intent(in)                ::  neig                            !<  number of bands (including spin)

  logical, intent(in)                ::  lpsi                            !<  also calculates the wave-function derivative and Berry curvature

  integer, intent(in)                ::  isort(mxddim)                   !<  G-vector corresponding to coefficient i of wavefunction
  real(REAL64), intent(in)           ::  ekpg(mxddim)                    !<  kinetic energy (hartree) of k+g-vector of row/column i

  integer, intent(in)                ::  nlevel                          !<  number of energy levels
  integer, intent(in)                ::  levdeg(mxdlev)                  !<  degeneragy of level
  integer, intent(in)                ::  leveigs(mxdlev,mxddeg)          !<  points to degenerate level

  complex(REAL64), intent(in)        ::  psi_sp(2*mxddim, mxdbnd)        !<  |psi_sp> (in principle eigen-functions)
  real(REAL64), intent(in)           ::  ei(mxdbnd)                      !<  eigenvalue no. i. (hartree)

  integer, intent(in)                ::  ng                              !<  size of g-space
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  G-vectors in reciprocal lattice coordinates

  integer, intent(in)                ::  kmscr(7)                        !<  max value of kgv(i,n) used for
  integer, intent(in)                ::  nsp                             !<  number of spin components ox xc-potential (1,2,4)
  real(REAL64), intent(in)           ::  vscr_sp(mxdscr,mxdnsp)          !<  screened local potential (including spin)

  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i
  real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)            !<  lattice coordinates of atom j of type i
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in real space

  integer, intent(in)                ::  nqnl(mxdtyp)                    !<  number of points for pseudo interpolation for atom k
  real(REAL64), intent(in)           ::  delqnl(mxdtyp)                  !<  step used in the pseudo interpolation for atom k
  real(REAL64), intent(in)        ::  vkb(-2:mxdlqp,0:LMAX,-1:1,mxdtyp)  !<  (1/q**l) * KB nonlocal pseudo. for atom k, ang. mom. l. NOT normalized to vcell, hartree
  integer, intent(in)                ::  nkb(0:LMAX,-1:1,mxdtyp)         !<  KB pseudo.  normalization for atom k, ang. mom. l

! output

  complex(REAL64), intent(out)       ::  dhdkpsi_sp(2*mxddim, mxdbnd, 3) !<  d H / d k | psi_sp>  (lattice coordinates)
  complex(REAL64), intent(out)       ::  dpsidk_sp(2*mxddim, mxdbnd, 3)  !<  d |psi_sp> / d k  (lattice coordinates)

  complex(REAL64), intent(out)       ::  psidhdkpsi_sp(mxddeg,mxddeg,3,mxdlev)  !<  <psi_sp_n| d H / d k |psi_sp_m>  for each energy level (lattice coordinates)

  complex(REAL64), intent(out)       ::  tfqg(3,3,mxddeg,mxddeg,mxdlev)      !<  Tensor F quantum geomtric, curvature and metric, (lattice coordinates)
  complex(REAL64), intent(out)       ::  tgammamf(3,3,mxddeg,mxddeg,mxdlev)  !<  Tensor Gamma-F, orbital magnetization and contribution to effective mass (lattice coordinates)
  complex(REAL64), intent(out)       ::  td2hdk2(3,3,mxddeg,mxddeg,mxdlev)   !<  <psi| d^2 H / d k^2 |psi> (lattice coordinates)

! local allocatable arrays

  real(REAL64), allocatable          ::  xnlkbsp(:)                      !  KB normalization with spin-orbit
  complex(REAL64), allocatable       ::  anlsp(:,:)                      !  KB projectors with spin-orbit

  complex(REAL64), allocatable       ::  danlspdrk(:,:,:)                !  d anlsp / d rkpt
  complex(REAL64), allocatable       ::  d2anlspdrk2(:,:,:,:)            !  d 2 anlsp / d 2 rkp (unused)

  complex(REAL64), allocatable       ::  hmedpsidk_sp(:,:,:)             !  (H - E_n) ( d | psi_sp_n > / d k )

  complex(REAL64), allocatable       ::  vnl0(:,:)                       !  <Psi|V_NL|Psi>
  complex(REAL64), allocatable       ::  dvnl0drk(:,:,:)                 !  d <Psi|V_NL|Psi> d k
  complex(REAL64), allocatable       ::  d2vnl0drk2(:,:,:,:)             !  d^2 <Psi|V_NL|Psi> d k^2

! local variables

  integer           ::  nder                                             !  order of derivative

  integer           ::  nanl, nanlso, nanlsp                             !  number of projectors
  integer           ::  mxdasp                                           !  maximum number of projectors
  real(REAL64)      ::  vcell, bdot(3,3)                                 !  cell volume, metric reciprocal space

! parameters

  real(REAL64), parameter       ::  ZERO = 0.0_REAL64
  real(REAL64), parameter       ::  TOL = 1.0E-12_REAL64

! counters

  integer    ::  i, j, n, m
  integer    ::  nl, nk, mk

! external

  complex(REAL64) ,external   :: zdotc


!  allocate(dmat(mxddeg,mxddeg,3,3))

! projector and derivatives

  call size_proj_nl_kb(ntype, natom, nkb, nanl, nanlso, nanlsp,          &
      mxdtyp)

  mxdasp = nanlsp

  allocate(anlsp(2*mxddim,mxdasp))
  allocate(xnlkbsp(mxdasp))

  allocate(danlspdrk(2*mxddim,mxdasp,3))
  allocate(d2anlspdrk2(2*mxddim,mxdasp,3,3))

  nder = 1
  if(lpsi) then
    nder = 2
  else
    nder = 1
  endif

  call proj_nl_kb_so_der_c16(rkpt, mtxd, isort,                          &
      nanlsp, nder,                                                      &
      ng, kgv,                                                           &
      nqnl, delqnl, vkb, nkb,                                            &
      ntype, natom, rat, adot,                                           &
      anlsp, xnlkbsp, danlspdrk, d2anlspdrk2,                            &
      mxdtyp, mxdatm, mxdlqp, mxddim, mxdasp, mxdgve)

  call berry_dhdk_psi_spin(rkpt, adot, mtxd, neig, psi_sp, dhdkpsi_sp,   &
      kgv, isort,                                                        &
      nanlsp, anlsp, xnlkbsp, danlspdrk,                                 &
      mxddim, mxdbnd, mxdgve, mxdasp)

  do nl = 1,nlevel
    do nk = 1,levdeg(nl)
      n = leveigs(nl,nk)
      do mk = 1,levdeg(nl)
        m = leveigs(nl,mk)
        do j = 1,3
          psidhdkpsi_sp(nk,mk,j,nl) = zdotc(2*mtxd, psi_sp(:,n), 1, dhdkpsi_sp(:,m,j), 1)
        enddo
      enddo
    enddo
  enddo

  if(lpsi) then

!   solves Sternheimer equation

    call berry_stern_solve_spin(mtxd, neig, psi_sp, ei, dhdkpsi_sp,      &
       dpsidk_sp, TOL,                                                   &
       nlevel, levdeg, leveigs,                                          &
       isort, ekpg,                                                      &
       vscr_sp, kmscr, nsp,                                              &
       ng, kgv,                                                          &
       nanlsp, anlsp, xnlkbsp,                                           &
       mxddim, mxdbnd, mxdgve, mxdscr, mxdasp, mxdlev, mxddeg, mxdnsp)

!   geometric quantities

    do nl = 1,nlevel
      do nk = 1,levdeg(nl)
        n = leveigs(nl,nk)
        do mk = 1,levdeg(nl)
          m = leveigs(nl,mk)

          do i = 1,3
          do j = 1,3
            tfqg(i,j,nk,mk,nl) = zdotc(2*mtxd, dpsidk_sp(:,n,i), 1, dpsidk_sp(:,m,j), 1)
          enddo
          enddo

        enddo
      enddo

    enddo

!   geometric*H quantities

    allocate(hmedpsidk_sp(2*mxddim,mxdbnd,3))

    do j = 1,3

      call hk_psi_spin_c16(mtxd, neig, dpsidk_sp(:,:,j), hmedpsidk_sp(:,:,j), .TRUE.,  &
             ng, kgv,                                                                &
             ekpg, isort, vscr_sp, kmscr, nsp,                                       &
             anlsp, xnlkbsp, nanlsp,                                                 &
             mxddim, mxdbnd, mxdasp, mxdgve, mxdscr, mxdnsp)

      do n = 1,neig

        call zaxpy(2*mtxd, cmplx(-ei(n),ZERO,REAL64), dpsidk_sp(:,n,j), 1,           &
                    hmedpsidk_sp(:,n,j), 1)

      enddo

    enddo

    allocate(vnl0(mxdbnd,mxdbnd))
    allocate(dvnl0drk(mxdbnd,mxdbnd,3))
    allocate(d2vnl0drk2(mxdbnd,mxdbnd,3,3))

!   extra calculations may be needed later

    call psi_vnl_psi_der(2*mtxd, neig, nanlsp, psi_sp, nder,             &
        vnl0, dvnl0drk, d2vnl0drk2,                                      &
        anlsp, xnlkbsp, danlspdrk, d2anlspdrk2,                          &
        2*mxddim, mxdbnd, mxdasp)

    call adot_to_bdot(adot, vcell, bdot)

    do nl = 1,nlevel
      do nk = 1,levdeg(nl)
        n = leveigs(nl,nk)
        do mk = 1,levdeg(nl)
          m = leveigs(nl,mk)

          do i = 1,3
          do j = 1,3
            tgammamf(i,j,nk,mk,nl) = zdotc(2*mtxd, dpsidk_sp(:,n,i), 1, hmedpsidk_sp(:,m,j), 1)
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

    deallocate(hmedpsidk_sp)

  endif

  deallocate(anlsp)
  deallocate(xnlkbsp)

  deallocate(danlspdrk)
  deallocate(d2anlspdrk2)

  return

end subroutine berry_derivative_spin
