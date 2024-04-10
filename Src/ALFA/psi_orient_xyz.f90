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

!>  Orients the degenerate wavefunctions on a degenerate subspace
!>  along the canonical xyz axes.
!>  Also chooses a phase for each wave-function
!>
!>  \author       Jose Luis Martins
!>  \version      5.11
!>  \date         30 October 2023.
!>  \copyright    GNU Public License v2

subroutine psi_orient_xyz(rkpt, adot, mtxd, neig, isort,                 &
    psi, ei,                                                             &
    ng, kgv,                                                             &
    mxddim, mxdbnd, mxdgve)

! Written 30 October 2023.  JLM
! Modified TOL, coeficients, overwrites psi. 10 April 2024. JLM



  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands
  integer, intent(in)                ::  mxdgve                          !<  array dimension of G-space vectors

  real(REAL64), intent(in)           ::  rkpt(3)                         !<  k-point reciprocal lattice coordinates
  integer, intent(in)                ::  mtxd                            !<  wavefunction dimension
  integer, intent(in)                ::  neig                            !<  number of eigenvalues

  integer, intent(in)                ::  isort(mxddim)                   !<  G-vector corresponding to coefficient i of wavefunction

  real(REAL64), intent(in)           ::  ei(mxdbnd)                      !<  eigenvalue no. i. (hartree)

  integer, intent(in)                ::  ng                              !<  size of g-space
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  G-vectors in reciprocal lattice coordinates
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space


! input and output

  complex(REAL64), intent(inout)     ::  psi(mxddim, mxdbnd)             !<  |psi_xyz> "aligned" on output psi within each degeneracy level

! local allocatable variables

  integer, allocatable               ::  levdeg(:)                       !  degeneracy of energy level
  integer, allocatable               ::  leveigs(:,:)                    !  states belonging to level

  complex(REAL64), allocatable       ::  aoper(:,:)                      !  hermitian operator
  complex(REAL64), allocatable       ::  aeigvec(:,:)                    !  eigenvectors
  real(REAL64), allocatable          ::  aeig(:)                         !  eigenvalues

  real(REAL64), allocatable          ::  ekin(:)                         !  kinetic energy

  real(real64), allocatable          ::  ekpg(:)                         !  kinetic energy sorted
  complex(REAL64), allocatable       ::  psi_tmp(:,:)                    !  temporary |psi>

! local variables

  integer           ::  mxdlev                          !  array dimension for number of levels
  integer           ::  mxddeg                          !  array dimension for number of levels

  integer           ::  nlevel                          !  number of energy levels

  real(REAL64)      ::  avec(3,3), bvec(3,3)
  real(REAL64)      ::  qk(3), qcar(3)
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

  call berry_degeneracy(.TRUE., neig, neigtmp, ei, TOL,                  &
         nlevel, maxdeg, levdeg, leveigs,                                &
         mxdbnd, 1, 1)

  mxdlev = nlevel
  mxddeg = maxdeg

  deallocate(levdeg)
  deallocate(leveigs)

  allocate(levdeg(mxdlev))
  allocate(leveigs(mxdlev,mxddeg))

  call berry_degeneracy(.FALSE., neig, neigtmp, ei, TOL,                 &
         nlevel, maxdeg, levdeg, leveigs,                                &
         mxdbnd, mxdlev, mxddeg)

! alignment with xyz

  allocate(aoper(mxddeg,mxddeg))
  allocate(aeigvec(mxddeg,mxddeg))
  allocate(aeig(mxddeg))

  allocate(psi_tmp(mxddim,mxddeg))

  allocate(ekin(ng))

  allocate(ekpg(mxddim))

  call adot_to_avec_sym(adot,avec,bvec)

  do i=1,ng
    qk(1) = rkpt(1) + kgv(1,i)
    qk(2) = rkpt(2) + kgv(2,i)
    qk(3) = rkpt(3) + kgv(3,i)
    qcar(1) = bvec(1,1)*qk(1) + bvec(1,2)*qk(2) + bvec(1,3)*qk(3)
    qcar(2) = bvec(2,1)*qk(1) + bvec(2,2)*qk(2) + bvec(2,3)*qk(3)
    qcar(3) = bvec(3,1)*qk(1) + bvec(3,2)*qk(2) + bvec(3,3)*qk(3)
    ekin(i) = 1.11*qcar(1)*qcar(1) + 2.03*qcar(2)*qcar(2) + 4.05*qcar(3)*qcar(3)
    ekin(i) = ekin(i) / 2
  enddo

  do i = 1,mtxd
    ekpg(i) = ekin(isort(i))
  enddo

  do nl = 1,nlevel

    if(levdeg(nl) /= 1) then

      do n = 1,levdeg(nl)
        nk = leveigs(nl,1) + n - 1
        psi_tmp(:,n) = psi(:,nk)
      enddo

      call psi_kin_psi(mtxd, levdeg(nl), psi_tmp, aoper, ekpg,           &
           mxddim, mxddeg)

      call diag_c16(levdeg(nl), aoper, aeig, aeigvec, mxddeg, info)
      if(info /= 0) then

        write(6,*) "   STOPPED in psi_orient_xyz  info = ",info

        stop

      endif

      call zgemm('n','n', mtxd, levdeg(nl), levdeg(nl),                  &
           C_UM, psi_tmp, mxddim, aeigvec, mxddeg,                       &
           C_ZERO,psi(:,leveigs(nl,1)), mxddim)

    endif

  enddo

! chooses a phase

  do nl = 1,nlevel
    do nk = 1,levdeg(nl)

      n = leveigs(nl,nk)

      xmax = real(psi(1,n)*conjg(psi(1,n)),REAL64)
      imax = 1

      do i = 1,mtxd
        tmax = real(psi(i,n)*conjg(psi(i,n)),REAL64)
        if(tmax > xmax + TOL) then
          imax = i
          xmax = tmax
        endif
      enddo

      xmax = real(psi(imax,n)*conjg(psi(imax,n)),REAL64)
      xmax = sqrt(xmax)
      phase = conjg(psi(imax,n)) / xmax
      psi(:,n) = phase*psi(:,n)

    enddo
  enddo

  deallocate(aoper)
  deallocate(aeigvec)
  deallocate(aeig)

  deallocate(ekin)
  deallocate(ekpg)

  deallocate(psi_tmp)

  return

end subroutine psi_orient_xyz
