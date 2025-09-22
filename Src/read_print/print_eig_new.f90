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

!>  Prints the eigenvalues and if requested the eigenvectors
!>  adapted from Sverre Froyen plane wave program
!>
!>  \author       Jose Luis Martins
!>  \version      5.12
!>  \date         June 6 1987. 21 September 2025.
!>  \copyright    GNU Public License v2

subroutine print_eig_new(ipr, nspin, lrk, irk, labelk, nrka, rkpt,       &
    mtxd, icmplx, neig, psi,                                             &
    adot, ei, ekpsi, isort, kgv,                                         &
    mxddim, mxdbnd, mxdgve)

! written June 6 1987. jlm
! modified February 18 1990. jlm
! version 4.0. 17 October 93. jlm
! modified 22 March 1999. jlm
! modified 31 March 2004. jlm
! modified (f90) 14 January 2014. jlm
! Modified, documentation, August 2019.
! Modified, neig > 9999, 9 November 2020. JLM
! Placement of complex components. 10 December 2022. JLM
! Merged with/without spin-orbit versions. 21 September 2025. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands
  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors

  integer, intent(in)                ::  nspin                           !<  spin components (1:no spin or 2:spin present)

  integer, intent(in)                ::  ipr                             !<  print option, 0 no printing, 1 eigenvalues, 2 eigenv+kin energy, >2 prints eigenvalues and ipr componets of eigenvectors
  logical, intent(in)                ::  lrk                             !<  if true prints explanation of rk coordiantes

  integer, intent(in)                ::  irk                             !<  number of the reduced k point
  character(len=5), intent(in)       ::  labelk                          !<  label of the k-point (optional, see nrka)
  integer, intent(in)                ::  nrka                            !<  if nrka > 0, point has known label
  real(REAL64), intent(in)           ::  rkpt(3)                         !<  j-th component in lattice coordinates of the k-point
  integer, intent(in)                ::  mtxd                            !<  dimension of the hamiltonian
  integer, intent(in)                ::  icmplx                          !<  indicates if the structure factor is complex
  integer, intent(in)                ::  neig                            !<  number of eigenvectors required (maybe modified on output)
  complex(REAL64), intent(in)        ::  psi(nspin*mxddim,nspin*mxdbnd)  !<  component j of eigenvector i (guess on input)

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space
  real(REAL64), intent(in)           ::  ei(nspin*mxdbnd)                !<  eigenvalue no. i. (hartree)
  real(REAL64), intent(in)           ::  ekpsi(nspin*mxdbnd)             !<  kinetic energy of eigenvector i. (hartree)
  integer, intent(in)                ::  isort(mxddim)                   !<  g-vector associated with row/column i of hamiltonian
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length

! local varaibles

  real(REAL64)        ::  avec(3,3),bvec(3,3),rkcar(3)
  character(len=60)   ::  prform

! parameters

  real(REAL64), parameter  ::  ZERO = 0.0_REAL64
  real(REAL64), parameter  ::  EV = 27.21138505_REAL64
  real(REAL64), parameter  ::  BOHR = 0.5291772109_REAL64

! counters

  integer    ::  i, j, jmax, nc


  if (ipr == 0) return

  if(nspin /=1 .and. nspin /= 2) then
    write(6,*)
    write(6,*) '    STOPPED in print_eig_new, nspin = ', nspin
    write(6,*)

    STOP

  endif

! calculates coordinates

  call adot_to_avec_sym(adot,avec,bvec)

  do i=1,3
    rkcar(i) = ZERO
    do j=1,3
      rkcar(i) = rkcar(i) + bvec(i,j)*rkpt(j)
    enddo
  enddo

  if (ipr > 2) then
    write(6,*)
    write(6,*)
  endif

  if (irk == 1) then
    if(lrk) then
      write(6,*)
      write(6,*)
      write(6,'("      k   mtxd      en(k) ",73x,"r(k) cartesian",18x,     &
         &     "r(k) lattice",18x,"r(k) cart Ang-1")')
      write(6,*)
    else
      write(6,*)
      write(6,*)
      write(6,'("       k  mtxd      en(k) ")')
    endif
  endif

  do i=1,nspin*neig,8

!   eigenvalues

    jmax = i+7
    if (jmax > nspin*neig) jmax = nspin*neig
    if(nrka > 0) then
      if (i == 1)  write(6,'(4x,a5,21x,i7,3x,8f9.5)')                    &
                        labelk,mtxd,(EV*ei(j),j=1,jmax)
    else
      if (i == 1) then
        write(prform,"( '(1x,i7,1x,i7,2x,',i1,'f9.5,5x,3f7.2,3x,3f11.4,3x,3f8.3)' )")  jmax
        write(6,prform) irk,mtxd,(EV*ei(j),j=1,jmax), (rkcar(j),j=1,3),  &
               (rkpt(j),j=1,3), (rkcar(j)/BOHR,j=1,3)
        if(ipr > 2) write(6,*)
      endif
    endif
    if (i > 1) then
      if(ipr > 2) write(6,*)
      write(6,'(18x,8f9.5)') (EV*ei(j),j=i,jmax)
      if(ipr > 2) write(6,*)
    endif

!   kinetic energies

    if (ipr == 2)  write(6,'(13x,"ekin ",8f9.5)') (EV*ekpsi(j),j=i,jmax)

!   eigenvectors, ncmax first components

    if (ipr > 2) then
      do nc = 1,ipr
        if(nspin == 2) then
          write(6,'(1x,i5,2x,3i3,1x,8f9.5)') nc,(kgv(j,isort(nc)),j=1,3),    &
                                 (real(psi(2*nc-1,j)),j=i,jmax)
          write(6,'(18x,8f9.5)') (aimag(psi(2*nc-1,j)),j=i,jmax)
          write(6,'(18x,8f9.5)') (real(psi(2*nc  ,j)),j=i,jmax)
          write(6,'(18x,8f9.5)') (aimag(psi(2*nc  ,j)),j=i,jmax)
        else
          write(6,'(1x,i5,2x,3i3,1x,8f9.5)') nc, (kgv(j,isort(nc)),j=1,3),   &
                                   (real(psi(nc,j)),j=i,jmax)
          if (icmplx /= 0) then
            write(6,'(18x,8f9.5)') (aimag(psi(nc,j)),j=i,jmax)
          endif
        endif
      enddo
    endif
  enddo

  write(6,*)
  if(nspin == 2) write(6,*)

  return

end subroutine print_eig_new


subroutine print_eig_so(ipr, irk, labelk, nrka, rkpt,                    &
    mtxd, neig, psi,                                                     &
    adot, ei, ekpsi, isort, kgv,                                         &
    mxddim, mxdbnd, mxdgve)


! written 21 September 2025. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands
  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors

  integer, intent(in)                ::  ipr                             !<  print option, 0 no printing, 1 eigenvalues, 2 eigenv+kin energy, 3 eigenvalues and eigenvectors
  integer, intent(in)                ::  irk                             !<  number of the reduced k point
  character(len=5), intent(in)       ::  labelk                          !<  label of the k-point (optional, see nrka)
  integer, intent(in)                ::  nrka                            !<  if nrka > 0, point has known label
  real(REAL64), intent(in)           ::  rkpt(3)                         !<  j-th component in lattice coordinates of the k-point

  integer, intent(in)                ::  mtxd                            !<  dimension of the hamiltonian
  integer, intent(in)                ::  neig                            !<  number of eigenvectors required (maybe modified on output)
  complex(REAL64), intent(in)        ::  psi(2*mxddim,2*mxdbnd)          !<  component j of eigenvector i (guess on input)
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space
  real(REAL64), intent(in)           ::  ei(2*mxdbnd)                    !<  eigenvalue no. i. (hartree)
  real(REAL64), intent(in)           ::  ekpsi(2*mxdbnd)                 !<  kinetic energy of eigenvector i. (hartree)
  integer, intent(in)                ::  isort(mxddim)                   !<  g-vector associated with row/column i of hamiltonian
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length

  integer    ::  iprnew
  iprnew = ipr
  if(iprnew > 2) iprnew = 20

  call print_eig_new(iprnew, 2, .FALSE., irk, labelk, nrka, rkpt,        &
      mtxd, 1, neig, psi,                                                &
      adot, ei, ekpsi, isort, kgv,                                       &
      mxddim, mxdbnd, mxdgve)


  return

end subroutine print_eig_so


subroutine print_eig(ipr, irk, labelk, nrka, rkpt,                       &
    mtxd, icmplx, neig, psi,                                             &
    adot, ei, ekpsi, isort, kgv,                                         &
    mxddim,mxdbnd,mxdgve)

! written 21 September 2025. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands
  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors

  integer, intent(in)                ::  ipr                             !<  print option, 0 no printing, 1 eigenvalues, 2 eigenv+kin energy, 3 eigenvalues and eigenvectors
  integer, intent(in)                ::  irk                             !<  number of the reduced k point
  character(len=5), intent(in)       ::  labelk                          !<  label of the k-point (optional, see nrka)
  integer, intent(in)                ::  nrka                            !<  if nrka > 0, point has known label
  real(REAL64), intent(in)           ::  rkpt(3)                         !<  j-th component in lattice coordinates of the k-point

  integer, intent(in)                ::  mtxd                            !<  dimension of the hamiltonian
  integer, intent(in)                ::  icmplx                          !<  indicates if the structure factor is complex
  integer, intent(in)                ::  neig                            !<  number of eigenvectors required (maybe modified on output)
  complex(REAL64), intent(in)        ::  psi(mxddim,mxdbnd)              !<  component j of eigenvector i (guess on input)

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space
  real(REAL64), intent(in)           ::  ei(mxdbnd)                      !<  eigenvalue no. i. (hartree)
  real(REAL64), intent(in)           ::  ekpsi(mxdbnd)                   !<  kinetic energy of eigenvector i. (hartree)
  integer, intent(in)                ::  isort(mxddim)                   !<  g-vector associated with row/column i of hamiltonian
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length

  integer    ::  iprnew
  iprnew = ipr
  if(iprnew > 2) iprnew = 20

  call print_eig_new(iprnew, 1, .TRUE., irk, labelk, nrka, rkpt,         &
      mtxd, icmplx, neig, psi,                                           &
      adot, ei, ekpsi, isort, kgv,                                       &
      mxddim, mxdbnd, mxdgve)

  return

end subroutine print_eig


