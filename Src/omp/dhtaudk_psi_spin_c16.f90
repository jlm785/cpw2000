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

!>  Spinor version of dhtaudk_psi_c16.  vtau_sp has the same
!>  convention as vscr_sp, components (1, sigma_z, sigma_x, sigma_y),
!>  so the vtau term is  -1/2 nabla . (vtau_1 + vtau . sigma) nabla.
!>  As the operator is linear in vtau, the result is the Pauli
!>  combination of dhtaudk_psi_c16 applied to the two spinor
!>  components (interleaved, 2*i-1 and 2*i) for each component of vtau_sp.
!>
!>  \author       claude
!>  \version      5.13
!>  \date         8 October 2026, 13 October 2026.
!>  \copyright    GNU Public License v2

subroutine dhtaudk_psi_spin_c16(mtxd, neig, psi_sp, vtaupsi_sp,          &
    dhtaupsi_sp,                                                         &
    rkpt, adot, isort, kgv, vtau_sp, kmscr, nsp,                         &
    mxddim, mxdbnd, mxdgve, mxdscr, mxdnsp)

! Written 8 October 2026. claude
! vtau_sp with spin components (1,2,4) as vscr_sp. 12 October 2026. claude
! OpenMP directives. 13 October 2026. JLM+claude
! Renamed from dvtaudk_psi_spin_c16, the derivative is of H_tau, vtau does not depend on k. 13 October 2026. JLM+claude

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands
  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxdscr                          !<  array dimension of vtau_sp
  integer, intent(in)                ::  mxdnsp                          !<  array dimension for number of spin components (1,2,4)

  integer, intent(in)                ::  mtxd                            !<  wavefunction dimension (basis size, does not include spin)
  integer, intent(in)                ::  neig                            !<  number of spin-wavefunctions
  complex(REAL64), intent(in)        ::  psi_sp(2*mxddim,mxdbnd)         !<  spin-wavevectors

  real(REAL64), intent(in)           ::  rkpt(3)                         !<  component in lattice coordinates of the k-point
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space
  integer, intent(in)                ::  isort(mxddim)                   !<  g-vector associated with row/column i of hamiltonian
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length

  real(REAL64), intent(in)           ::  vtau_sp(mxdscr,mxdnsp)          !<  d (rho eps_xc) / d tau in the fft real space mesh (1, sigma_z, sigma_x, sigma_y)
  integer, intent(in)                ::  kmscr(7)                        !<  max value of kgv(i,n) used for the potential fft mesh and fft mesh size
  integer, intent(in)                ::  nsp                             !<  number of spin components of vtau_sp (1,2,4)

! output

  complex(REAL64), intent(out)       ::  vtaupsi_sp(2*mxddim,mxdbnd)     !<  vtau |psi_sp>
  complex(REAL64), intent(out)       ::  dhtaupsi_sp(2*mxddim,mxdbnd,3)  !<  (d H_tau / d k_l) |psi_sp>, k in lattice coordinates

! local allocatable arrays

  complex(REAL64), allocatable       ::  psi(:,:)                        !  spin components as independent wave-functions (up: 1..neig, down: neig+1..2*neig)
  complex(REAL64), allocatable       ::  vtaupsi(:,:)                    !  result for one component of vtau_sp
  complex(REAL64), allocatable       ::  dhtaupsi(:,:,:)                 !  result for one component of vtau_sp

! local variables

  complex(REAL64)    ::  cup, cdn                                        !  coefficients of the up and down results in the up row (cup) and in the down row (cdn)

! constants

  real(REAL64), parameter     ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)
  complex(REAL64), parameter  ::  C_I = cmplx(ZERO,UM,REAL64)

! counters

  integer   ::   m, n, l, isp


  if(nsp /= 1 .and. nsp /= 2 .and. nsp /= 4) then
    write(6,*)
    write(6,'("   STOPPED in dhtaudk_psi_spin_c16.  nsp = ",i8,          &
           &  " is an incorrect value, should be 1,2,4")') nsp
    write(6,*)

    stop

  endif

  allocate(psi(mtxd,2*neig))
  allocate(vtaupsi(mtxd,2*neig))
  allocate(dhtaupsi(mtxd,2*neig,3))

!$omp parallel do default(shared) private(n, m)
  do n = 1,neig
    do m = 1,mtxd
      psi(m,n     ) = psi_sp(2*m-1,n)
      psi(m,n+neig) = psi_sp(2*m  ,n)
    enddo
  enddo
!$omp end parallel do

!$omp parallel do default(shared) private(n, m, l)
  do n = 1,neig
    do m = 1,2*mtxd
      vtaupsi_sp(m,n) = C_ZERO
    enddo
    do l = 1,3
      do m = 1,2*mtxd
        dhtaupsi_sp(m,n,l) = C_ZERO
      enddo
    enddo
  enddo
!$omp end parallel do

  do isp = 1,nsp

    call dhtaudk_psi_c16(mtxd, 2*neig, psi, vtaupsi, dhtaupsi,           &
        rkpt, adot, isort, kgv, vtau_sp(:,isp), kmscr,                   &
        mtxd, 2*neig, mxdgve, mxdscr)

!   (1, sigma_z): diagonal in spin, +1/+1 and +1/-1
!   (sigma_x, sigma_y): off-diagonal, up row gets 1 or -i times down, down row 1 or +i times up

    if(isp == 1 .or. isp == 2) then

      cup = C_ZERO
      cdn = C_ZERO
      if(isp == 1) cdn = cmplx(UM,ZERO,REAL64)
      if(isp == 2) cdn = cmplx(-UM,ZERO,REAL64)

!$omp parallel do default(shared) private(n, m, l)
      do n = 1,neig
        do m = 1,mtxd
          vtaupsi_sp(2*m-1,n) = vtaupsi_sp(2*m-1,n) + vtaupsi(m,n     )
          vtaupsi_sp(2*m  ,n) = vtaupsi_sp(2*m  ,n) +                    &
                                      cdn*vtaupsi(m,n+neig)
        enddo
        do l = 1,3
          do m = 1,mtxd
            dhtaupsi_sp(2*m-1,n,l) = dhtaupsi_sp(2*m-1,n,l) +            &
                                      dhtaupsi(m,n     ,l)
            dhtaupsi_sp(2*m  ,n,l) = dhtaupsi_sp(2*m  ,n,l) +            &
                                      cdn*dhtaupsi(m,n+neig,l)
          enddo
        enddo
      enddo
!$omp end parallel do

    else

      if(isp == 3) then
        cup = cmplx(UM,ZERO,REAL64)
        cdn = cmplx(UM,ZERO,REAL64)
      else
        cup = -C_I
        cdn = C_I
      endif

!$omp parallel do default(shared) private(n, m, l)
      do n = 1,neig
        do m = 1,mtxd
          vtaupsi_sp(2*m-1,n) = vtaupsi_sp(2*m-1,n) +                    &
                                      cup*vtaupsi(m,n+neig)
          vtaupsi_sp(2*m  ,n) = vtaupsi_sp(2*m  ,n) +                    &
                                      cdn*vtaupsi(m,n     )
        enddo
        do l = 1,3
          do m = 1,mtxd
            dhtaupsi_sp(2*m-1,n,l) = dhtaupsi_sp(2*m-1,n,l) +            &
                                      cup*dhtaupsi(m,n+neig,l)
            dhtaupsi_sp(2*m  ,n,l) = dhtaupsi_sp(2*m  ,n,l) +            &
                                      cdn*dhtaupsi(m,n     ,l)
          enddo
        enddo
      enddo
!$omp end parallel do

    endif

  enddo

  deallocate(psi)
  deallocate(vtaupsi)
  deallocate(dhtaupsi)

  return

end subroutine dhtaudk_psi_spin_c16
