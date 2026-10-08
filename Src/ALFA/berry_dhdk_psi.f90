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

!>  Calculates the derivatives of V_NL|Psi>  + K |psi> for a separable non-local
!>  pseudopotential with respect to the k-vector.
!>
!>  \author       Jose Luis Martins, Carlos Loia Reis
!>  \version      5.13
!>  \date         9 January 2023. 8 October 2026.
!>  \copyright    GNU Public License v2


subroutine berry_dhdk_psi(rkpt, adot, mtxd, neig, psi, dhdkpsi,          &
    kgv, isort,                                                          &
    nanl, anlga, xnlkb, danlgadrk,                                       &
    lgks, vtaumsh, kmscr,                                                &
    mxddim, mxdbnd, mxdgve, mxdanl, mxdscr)

! adapted from psi_vnl_psi_der, psi_p_psi and CLR phonon hk_psi_nl_lr_c16
! Non-local part in a separate subroutine. 15 December 2023. JLM
! Generalized Kohn-Sham meta-GGA, vtau term. 8 October 2026. JLM+claude


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands
  integer, intent(in)                ::  mxdgve                          !<  array dimension of G-space vectors
  integer, intent(in)                ::  mxdanl                          !<  array dimension of number of projectors
  integer, intent(in)                ::  mxdscr                          !<  array dimension of vtaumsh

  real(REAL64), intent(in)           ::  rkpt(3)                         !<  k-point reciprocal lattice coordinates
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in real space

  integer, intent(in)                ::  mtxd                            !<  wavefunction dimension
  integer, intent(in)                ::  neig                            !<  wavefunction dimension

  complex(REAL64), intent(in)        ::  psi(mxddim, mxdbnd)             !<  |psi> (in principle eigen-functions)

  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  G-vectors in reciprocal lattice coordinates
  integer, intent(in)                ::  isort(mxddim)                   !<  G-vector corresponding to coefficient i of wavefunction

  integer, intent(in)                ::  nanl                            !<  half of number of projectors without spin
  complex(REAL64), intent(in)        ::  anlga(mxddim,mxdanl)            !<  KB projectors without spin-orbit
  real(REAL64), intent(in)           ::  xnlkb(mxdanl)                   !<  KB normalization without spin-orbit
  complex(REAL64), intent(in)        ::  danlgadrk(mxddim,mxdanl,3)      !<  d anlga / d rkpt

  logical, intent(in)                ::  lgks                            !<  generalized Kohn-Sham meta-GGA (vtau term)
  real(REAL64), intent(in)           ::  vtaumsh(mxdscr)                 !<  d (rho eps_xc) / d tau in the fft real space mesh (only used if lgks)
  integer, intent(in)                ::  kmscr(7)                        !<  max value of kgv(i,n) used for the potential fft mesh and fft mesh size

! output

  complex(REAL64), intent(out)       ::  dhdkpsi(mxddim,mxdbnd,3)        !<  (d H /d k) |Psi>

! local allocatable variables

  real(REAL64), allocatable          ::  qcontra(:,:)                    !  contravariant rkpt+kgv
  complex(REAL64), allocatable       ::  vtaupsi(:,:)                    !  vtau |psi>
  complex(REAL64), allocatable       ::  dvtaupsi(:,:,:)                 !  (d H_tau / d k) |psi>
  real(REAL64), allocatable          ::  bdot(:,:)                       !  metric in reciprocal space.  It is allocatable to keep subroutine thread-safe

! local variables

  real(REAL64)      ::  qk1, qk2, qk3

  real(REAL64)      ::  vcell

! constants

  real(REAL64), parameter     ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)
  complex(REAL64), parameter  ::  C_UM = cmplx(UM,ZERO,REAL64)

! counters

  integer    ::  j, m, n

! calculates the non local part

  call berry_dvnldk_psi(mtxd, neig, psi, dhdkpsi,                        &
    nanl, anlga, xnlkb, danlgadrk,                                       &
    mxddim, mxdbnd, mxdanl)

! kinetic energy operator derivative.  momentum operator in contravariant coordinates

  allocate(qcontra(3,mtxd))
  allocate(bdot(3,3))

  call adot_to_bdot(adot,vcell,bdot)

  do m = 1,mtxd
    qk1 = rkpt(1) + UM*(kgv(1,isort(m)))
    qk2 = rkpt(2) + UM*(kgv(2,isort(m)))
    qk3 = rkpt(3) + UM*(kgv(3,isort(m)))
    do j = 1,3
      qcontra(j,m) = bdot(j,1)*qk1 + bdot(j,2)*qk2 + bdot(j,3)*qk3
    enddo
  enddo

  do j = 1,3
    do n = 1,neig
      do m = 1,mtxd
        dhdkpsi(m,n,j) = dhdkpsi(m,n,j) + qcontra(j,m)*psi(m,n)
      enddo
    enddo
  enddo

  deallocate(qcontra)
  deallocate(bdot)

! vtau term of the generalized Kohn-Sham meta-GGA

  if(lgks) then

    allocate(vtaupsi(mxddim,mxdbnd))
    allocate(dvtaupsi(mxddim,mxdbnd,3))

    call dvtaudk_psi_c16(mtxd, neig, psi, vtaupsi, dvtaupsi,              &
        rkpt, adot, isort, kgv, vtaumsh, kmscr,                          &
        mxddim, mxdbnd, mxdgve, mxdscr)

    do j = 1,3
      do n = 1,neig
        do m = 1,mtxd
          dhdkpsi(m,n,j) = dhdkpsi(m,n,j) + dvtaupsi(m,n,j)
        enddo
      enddo
    enddo

    deallocate(vtaupsi)
    deallocate(dvtaupsi)

  endif

  return

end subroutine berry_dhdk_psi

