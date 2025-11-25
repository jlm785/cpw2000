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

!>  Computes, symmetrizes, and adds the "kinetic energy density"
!>  from the eigenvectors at a given k-point
!>
!>  \author       Carlos Loia Reis, José Luís Martins
!>  \version      5.12
!>  \date         2 October 2015, 25 November 2025.
!>  \copyright    GNU Public License v2

subroutine tau_by_fft(tauk, mtxd, neig, occp, isort, psi,                &
     rkpt, adot,                                                         &
     ng, kgv, phase, conj, ns, inds, kmax, mstar,                        &
     mxddim, mxdbnd, mxdgve, mxdnst)

! Adapted from charge_by_fft.
! written 2 October 2015. JLM, CLR.
! Modified documentation, January 2020. JLM
! Corrected openmp bug, 30 September 2022. JLM
! Remove the doubling of kinetic energy (as in libxc). 25 November 2025. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands
  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxdnst                          !<  array dimension for g-space stars

  integer, intent(in)                ::  mtxd                            !<  dimension of the hamiltonian
  integer, intent(in)                ::  neig                            !<  number of eigenvectors
  real(REAL64), intent(in)           ::  occp(mxdbnd)                    !<  fractional ocupation of level j
  integer, intent(in)                ::  isort(mxddim)                   !<  g-vector associated with row/column i of hamiltonian

  complex(REAL64), intent(in)        ::  psi(mxddim,mxdbnd)              !<  |psi>

  real(REAL64), intent(in)           ::  rkpt(3)                         !<  k-point reciprocal lattice coordinates
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space

  integer, intent(in)                ::  ng                              !<  size of g-space
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  G-vectors in reciprocal lattice coordinates
  complex(REAL64), intent(in)        ::  phase(mxdgve)                   !<  phase factor of G-vector n
  real(REAL64), intent(in)           ::  conj(mxdgve)                    !<  is -1 if one must take the complex conjugate of x*phase
  integer, intent(in)                ::  ns                              !<  number os stars with length less than gmax
  integer, intent(in)                ::  inds(mxdgve)                    !<  star to which g-vector n belongs
  integer, intent(in)                ::  kmax(3)                         !<  max value of |kgv(i,n)|
  integer, intent(in)                ::  mstar(mxdnst)                   !<  number of g-vectors in the j-th star

! output

  complex(REAL64), intent(out)       ::  tauk(mxdnst)                    !<  symmetrized density (in stars of G)

! local allocatable arrays

  complex(REAL64), allocatable       ::  tauu(:)                         !  unsymmetrized density (in G)
  complex(REAL64), allocatable       ::  taumsh(:)
  real(REAL64), allocatable          ::  wrkfft(:)
  complex(REAL64), allocatable       ::  chd(:,:)
  integer, allocatable               ::  ipoint(:)

! local variables

  integer         ::  mxdfft,mxdwrk
  integer         ::  jmin, jmax
  integer         ::  nsfft(3)
  integer         ::  n1, n2, n3, id
  integer         ::  nn1, nn2, nn3
  integer         ::  kd1, kd2, kd3
  integer         ::  k1, k2, k3
  integer         ::  ntot, it
  integer         ::  iadd
  complex(REAL64) ::  xp

  real(REAL64)    ::  qk(3),qcar(3)
  real(REAL64)    ::  avec(3,3),bvec(3,3)

! counters

  integer         ::  i, j

! parameters

  real(REAL64), parameter :: SMALL = 0.000000000001_REAL64
  real(REAL64), parameter :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)
  complex(REAL64), parameter  ::  C_I = cmplx(ZERO,UM,REAL64)


  call adot_to_avec(adot,avec,bvec)

  do i=1,ns
    tauk(i) = C_ZERO
  enddo

  if (neig > 0) then

!   find min and max band with nonzero occupation

    jmin = 0
    jmax = 0
    do i=1,neig
      if (abs(occp(i)) > SMALL .and. jmin == 0) jmin = i
      if (abs(occp(i)) > SMALL .and. jmin /= 0) jmax = i
    enddo

    if (jmin > 0) then

!     find n for fast fourier transform
!     ni is the number of points used in direction i.

      call size_fft(kmax,nsfft,mxdfft,mxdwrk)

      allocate(chd(mxdfft,3))
      allocate(wrkfft(mxdwrk))

      n1 = nsfft(1)
      n2 = nsfft(2)
      n3 = nsfft(3)
      id = n1
      ntot = id * n2 * n3
      nn1 = (n1-1) / 2
      nn2 = (n2-1) / 2
      nn3 = (n3-1) / 2
      kd1 = 0
      kd2 = 0
      kd3 = 0

!     fills array ipoint (gather-scatter index)

      allocate(ipoint(mtxd))

      do i = 1,mtxd
        it = isort(i)
        k1 = kgv(1,it)
        if (iabs(k1) > kd1) kd1 = iabs(k1)
        if (k1 < 0) k1 = n1 + k1
        k2 = kgv(2,it)
        if (iabs(k2) > kd2) kd2 = iabs(k2)
        if (k2 < 0) k2 = n2 + k2
        k3 = kgv(3,it)
        if (iabs(k3) > kd3) kd3 = iabs(k3)
        if (k3 < 0) k3 = n3 + k3
        ipoint(i) = (k3*n2 + k2)*id + k1 + 1
      enddo
      if (kd1 > nn1 .or. kd2 > nn2 .or. kd3 > nn3) then
        write(6,*)
        write(6,'("     STOPPED in tau_by_fft.  size of matrix:    ",    &
           &      i7," fft mesh ",3i5)') mtxd,n1,n2,n3

        stop

      endif

!     initialize rhomsh

      allocate(taumsh(mxdfft))

      do i = 1,ntot
        taumsh(i) = ZERO
      enddo

!     start loop over eigenvectors

      do j = jmin,jmax

!$omp parallel do default(shared) private(i)
        do i=1,ntot
          chd(i,1) = C_ZERO
          chd(i,2) = C_ZERO
          chd(i,3) = C_ZERO
        enddo
!$omp end parallel do


!$omp parallel do default(shared) private(i,qk,qcar)
        do i=1,mtxd
          qk(1) = rkpt(1) + dble(kgv(1,isort(i)))
          qk(2) = rkpt(2) + dble(kgv(2,isort(i)))
          qk(3) = rkpt(3) + dble(kgv(3,isort(i)))
          qcar(1) = bvec(1,1)*qk(1) + bvec(1,2)*qk(2) + bvec(1,3)*qk(3)
          qcar(2) = bvec(2,1)*qk(1) + bvec(2,2)*qk(2) + bvec(2,3)*qk(3)
          qcar(3) = bvec(3,1)*qk(1) + bvec(3,2)*qk(2) + bvec(3,3)*qk(3)

          chd(ipoint(i),1) = C_I*qcar(1)*psi(i,j)
          chd(ipoint(i),2) = C_I*qcar(2)*psi(i,j)
          chd(ipoint(i),3) = C_I*qcar(3)*psi(i,j)
        enddo
!$omp end parallel do

!       Fourier transform to real space

        call cfft_wf_c16(chd(:,1), id,n1,n2,n3, kd1,kd2,kd3, -1, wrkfft,mxdwrk)
        call cfft_wf_c16(chd(:,2), id,n1,n2,n3, kd1,kd2,kd3, -1, wrkfft,mxdwrk)
        call cfft_wf_c16(chd(:,3), id,n1,n2,n3, kd1,kd2,kd3, -1, wrkfft,mxdwrk)

!       Calculates the square of chd (wave-function -> charge)

!$omp parallel do default(shared) private(i,xp)
        do i=1,ntot
          xp = conjg(chd(i,1))*chd(i,1) + conjg(chd(i,2))*chd(i,2)  &
                                        + conjg(chd(i,3))*chd(i,3)
          taumsh(i) = taumsh(i) + occp(j)*real(xp,REAL64) / 2
        enddo
!$omp end parallel do

      enddo                  !  end loop over eigenvectors


!     Fourier transform to momentum space

      deallocate(ipoint)

      call rfft_c16(taumsh, id,n1,n2,n3, 1, wrkfft,mxdwrk)

      allocate(tauu(mxdgve))

      do i = 1,ng
        tauu(i) = C_ZERO
      enddo

      do i=1,ng
        k1 = kgv(1,i)
        if (iabs(k1) <= nn1) then
        if (k1 < 0) k1=n1+k1
        k2 = kgv(2,i)
        if (iabs(k2) <= nn2) then
        if (k2 < 0) k2=n2+k2
        k3 = kgv(3,i)
        if (iabs(k3) <= nn3) then
        if (k3 < 0) k3=n3+k3
          iadd = (k3*n2 + k2)*id + k1 + 1
          tauu(i) = taumsh(iadd)
        endif
        endif
        endif
      enddo

      deallocate(chd)
      deallocate(wrkfft)
      deallocate(taumsh)


!     CONJ SHOULD BE CONVERTED TO INTEGER OR LOGICAL

      call star_of_g_fold(tauk, tauu, .FALSE.,                           &
           ng, phase, conj, ns, inds, mstar,                             &
           mxdgve, mxdnst)

      deallocate(tauu)

    endif

  endif

  return

end subroutine tau_by_fft
