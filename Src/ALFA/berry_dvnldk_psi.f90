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

!>  Calculates the derivatives of V_NL|Psi> for a separable non-local
!>  pseudopotential with respect to the k-vector.
!>
!>  \author       Jose Luis Martins, Carlos Loia Reis
!>  \version      5.09
!>  \date         9 January 2023.
!>  \copyright    GNU Public License v2


subroutine berry_dvnldk_psi(mtxd, neig, psi, dvnldkpsi,                  &
    nanl, anlga, xnlkb, danlgadrk,                                       &
    mxddim, mxdbnd, mxdanl)

! adapted from psi_vnl_psi_der, psi_p_psi and CLR phonon hk_psi_nl_lr_c16
! removed kinetic part from berry_dhdk_psi. 15 December 2023. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands
  integer, intent(in)                ::  mxdanl                          !<  array dimension of number of projectors

  integer, intent(in)                ::  mtxd                            !<  wavefunction dimension
  integer, intent(in)                ::  neig                            !<  wavefunction dimension

  complex(REAL64), intent(in)        ::  psi(mxddim, mxdbnd)             !<  |psi> (in principle eigen-functions)

  integer, intent(in)                ::  nanl                            !<  half of number of projectors without spin
  complex(REAL64), intent(in)        ::  anlga(mxddim,mxdanl)            !<  KB projectors without spin-orbit
  real(REAL64), intent(in)           ::  xnlkb(mxdanl)                   !<  KB normalization without spin-orbit
  complex(REAL64), intent(in)        ::  danlgadrk(mxddim,mxdanl,3)      !<  d anlga / d rkpt

! output

  complex(REAL64), intent(out)       ::  dvnldkpsi(mxddim,mxdbnd,3)      !<  (d V_NL /d k) |Psi>

! local allocatable variables


  complex(REAL64), allocatable       ::  dhd(:,:)                        !  < anl | psi >
  complex(REAL64), allocatable       ::  ddhddrk(:,:,:)                  !  < d anl / d rkpt | psi >

! constants

  real(REAL64), parameter     ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)
  complex(REAL64), parameter  ::  C_UM = cmplx(UM,ZERO,REAL64)

! counters

  integer    ::  j, m, n


  if(nanl > 0) then

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

!   dvnldkpsi =  | d anl / d k > Diag(xnl) < anl | psi > + | anl >  Diag(xnl) < d anl dk | psi >

    do j = 1,3

      call zgemm('n','n', mtxd, neig, nanl, C_UM, danlgadrk(:,:,j), mxddim,  &
                 dhd, nanl, C_ZERO, dvnldkpsi(:,:,j), mxddim)

      call zgemm('n','n', mtxd, neig, nanl, C_UM, anlga, mxddim,             &
                  ddhddrk(:,:,j), nanl, C_UM, dvnldkpsi(:,:,j), mxddim)

    enddo

    deallocate(dhd)
    deallocate(ddhddrk)

  else

     dvnldkpsi(:,:,:) = C_ZERO

  endif

  return

end subroutine berry_dvnldk_psi

