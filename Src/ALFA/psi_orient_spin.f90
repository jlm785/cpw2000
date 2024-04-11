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

!>  Orients the degenerate spin wavefunctions on a degenerate subspace
!>  with respect to spin.
!>
!>  \author       Jose Luis Martins
!>  \version      5.11
!>  \date         11 April 2024.
!>  \copyright    GNU Public License v2

subroutine psi_orient_spin(mtxd, neig, psi_sp, ei_sp,                    &
    mxddim, mxdbnd)

! Adapted form non-spin version. 11 April 2024. JLM



  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands

  integer, intent(in)                ::  mtxd                            !<  wavefunction dimension
  integer, intent(in)                ::  neig                            !<  number of eigenvalues

  real(REAL64), intent(in)           ::  ei_sp(2*mxdbnd)                 !<  eigenvalue no. i. (hartree)


! input and output

  complex(REAL64), intent(inout)     ::  psi_sp(2*mxddim, 2*mxdbnd)      !<  |psi_xyz> "aligned" on output psi within each degeneracy level

! local allocatable variables

  integer, allocatable               ::  levdeg(:)                       !  degeneracy of energy level
  integer, allocatable               ::  leveigs(:,:)                    !  states belonging to level

  complex(REAL64), allocatable       ::  aoper(:,:)                      !  hermitian operator
  complex(REAL64), allocatable       ::  aeigvec(:,:)                    !  eigenvectors
  real(REAL64), allocatable          ::  aeig(:)                         !  eigenvalues

  complex(REAL64), allocatable       ::  psi_tmp(:,:)                    !  temporary |psi>

! local variables

  integer           ::  mxdlev                          !  array dimension for number of levels
  integer           ::  mxddeg                          !  array dimension for number of levels

  integer           ::  nlevel                          !  number of energy levels

  integer           ::  maxdeg
  integer           ::  info

  real(REAL64)      ::  xmax, tmax
  integer           ::  imax
  complex(REAL64)   ::  phase

  integer           ::  neigtmp

! parameters

  real(REAL64), parameter       ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter    ::  C_UM = cmplx(UM,ZERO,REAL64)
  complex(REAL64), parameter    ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)
  real(REAL64), parameter       ::  TOL = 1.0E-6_REAL64

! counters

  integer    ::  i, n
  integer    ::  nl, nk


  neigtmp = neig

! identifies the degeneracies

  allocate(levdeg(1))
  allocate(leveigs(1,1))

  call berry_degeneracy(.TRUE., neig, neigtmp, ei_sp, TOL,               &
         nlevel, maxdeg, levdeg, leveigs,                                &
         2*mxdbnd, 1, 1)

  mxdlev = nlevel
  mxddeg = maxdeg

  deallocate(levdeg)
  deallocate(leveigs)

  allocate(levdeg(mxdlev))
  allocate(leveigs(mxdlev,mxddeg))

  call berry_degeneracy(.FALSE., neig, neigtmp, ei_sp, TOL,              &
         nlevel, maxdeg, levdeg, leveigs,                                &
         2*mxdbnd, mxdlev, mxddeg)

! prepares alignment with xyz

  allocate(aoper(mxddeg,mxddeg))
  allocate(aeigvec(mxddeg,mxddeg))
  allocate(aeig(mxddeg))

  allocate(psi_tmp(2*mxddim,mxddeg))

! loop over degenerate levels

  do nl = 1,nlevel

    if(levdeg(nl) /= 1) then

!     first diagonalizes sigma_z

      do n = 1,levdeg(nl)
        nk = leveigs(nl,1) + n - 1
        do i = 1,mtxd
          psi_tmp(2*i-1,n) =  psi_sp(2*i-1,nk)
          psi_tmp(2*i  ,n) = -psi_sp(2*i  ,nk)
        enddo
      enddo

      call zgemm('c','n', levdeg(nl), levdeg(nl), 2*mtxd,                &
           C_UM, psi_sp(1,leveigs(nl,1)), 2*mxddim, psi_tmp, 2*mxddim,   &
           C_ZERO, aoper, mxddeg)

      do n = 1,levdeg(nl)
        do i = 1,mtxd
          psi_tmp(2*i  ,n) = -psi_tmp(2*i  ,n)
        enddo
      enddo

      call diag_c16(levdeg(nl), aoper, aeig, aeigvec, mxddeg, info)
      if(info /= 0) then

        write(6,*) "   STOPPED in psi_orient_spin_xyz  info = ",info

        stop

      endif

      call zgemm('n','n', 2*mtxd, levdeg(nl), levdeg(nl),                &
           C_UM, psi_tmp, 2*mxddim, aeigvec, mxddeg,                     &
           C_ZERO,psi_sp(1,leveigs(nl,1)), 2*mxddim)

    endif

  enddo

! chooses a phase

  do nl = 1,nlevel
    do nk = 1,levdeg(nl)

      n = leveigs(nl,nk)

      xmax = real(psi_sp(1,n)*conjg(psi_sp(1,n)),REAL64)
      imax = 1

      do i = 1,mtxd
        tmax = real(psi_sp(i,n)*conjg(psi_sp(i,n)),REAL64)
        if(tmax > xmax + TOL) then
          imax = i
          xmax = tmax
        endif
      enddo

      xmax = real(psi_sp(imax,n)*conjg(psi_sp(imax,n)),REAL64)
      xmax = sqrt(xmax)
      phase = conjg(psi_sp(imax,n)) / xmax
      psi_sp(:,n) = phase*psi_sp(:,n)

    enddo
  enddo

  deallocate(aoper)
  deallocate(aeigvec)
  deallocate(aeig)

  deallocate(psi_tmp)

  return

end subroutine psi_orient_spin
