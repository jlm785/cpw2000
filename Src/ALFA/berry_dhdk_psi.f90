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
!>  with respect to the k-vector.
!>
!>  \author       Jose Luis Martins, Carlos Loia Reis
!>  \version      5.06
!>  \date         9 January 2023.
!>  \copyright    GNU Public License v2


subroutine berry_dhdk_psi(rkpt, adot, mtxd, neig, psi, dhdkpsi,          &
    kgv, isort,                                                          &
    nanl, anlga, xnlkb, danlgadrk,                                       &
    mxddim, mxdbnd, mxdgve, mxdanl)

! adapted from psi_vnl_psi_der, psi_p_psi and CLR phonon hk_psi_nl_lr_c16


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)
  integer, parameter          :: LMAX = 3                                !   hard coded max ang. mom.

! input

  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands
  integer, intent(in)                ::  mxdgve                          !<  array dimension of G-space vectors
  integer, intent(in)                ::  mxdanl                          !<  array dimension of number of projectors

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

! output

  complex(REAL64), intent(out)       ::  dhdkpsi(mxddim,mxdbnd,3)        !<  (d H /d k) |Psi>

! local allocatable variables


  complex(REAL64), allocatable       ::  dhd(:,:)                        !  < anl | psi >
  complex(REAL64), allocatable       ::  ddhddrk(:,:,:)                  !  < d anl / d rkpt | psi >

  real(REAL64), allocatable          ::  qcontra(:,:)                    !  contravariant rkpt+kgv
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


  allocate(dhd(nanl,neig))
  allocate(ddhddrk(nanl,neig,3))

  call zgemm('c','n', nanl, neig, mtxd, C_UM, anlga, mxddim, psi,        &
     &                mxddim, C_ZERO, dhd, nanl)

  do j = 1,3

    call zgemm('c','n', nanl, neig, mtxd, C_UM, danlgadrk(1,1,j),        &
               mxddim, psi, mxddim, C_ZERO, ddhddrk(1,1,j), nanl)

  enddo

  do n = 1,neig
    do m = 1,nanl
      dhd(m,n) = xnlkb(m)*dhd(m,n)
    enddo
  enddo

  do j = 1,3
    do n = 1,neig
      do m = 1,nanl
        ddhddrk(m,n,j) = xnlkb(m)*ddhddrk(m,n,j)
      enddo
    enddo
  enddo

! dhdkpsi =  | d anl / d k > Diag(xnl) < anl | psi > + | anl >  Diag(xnl) < d anl dk | psi >

  do j = 1,3

    call zgemm('n','n', mtxd, neig, nanl, C_UM, danlgadrk(:,:,j), mxddim,  &
               dhd, nanl, C_ZERO, dhdkpsi(:,:,j), mxddim)

    call zgemm('n','n', mtxd, neig, nanl, C_UM, anlga, mxddim,           &
                ddhddrk(:,:,j), nanl, C_UM, dhdkpsi(:,:,j), mxddim)

  enddo

  deallocate(dhd)
  deallocate(ddhddrk)

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

  return

end subroutine berry_dhdk_psi

