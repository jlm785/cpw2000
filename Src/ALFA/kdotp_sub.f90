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

!>  Tests the kdotp output file
!>
!>  \author       Jose Luis Martins
!>  \version      5.12
!>  \date         May 3, 2014, 13 August 2025.
!>  \copyright    GNU Public License v2


  subroutine kdotp_sub(ioreplay)

! Written May 3, 2014. JLM
! Modified, documentation, name, 15 June 2020. JLM
! Modified, last question, 30 September 2020. JLM
! Modified, neigso, new API kdotp_in. 8 September 2022. JLM
! Modified, ztot in out_band_circuit_size. 26 July 2024. JLM
! Modified, rk in out_band_eref, 13 August 2025. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


  integer, intent(in)         :: ioreplay                                !<  tape number for reproducing calculations

! dimensions

  integer                            ::  mxdbnd                          !  array dimension for the number of bands

! main variables

  character(len=64)                  ::  filename                        !  file to be written
  integer                            ::  iotape                          !  tape number

  character(len=50)                  ::  title                           !  title for plots
  character(len=140)                 ::  subtitle                        !  subtitle for plots

  integer                            ::  neig                            !  number of bands
  real(REAL64)                       ::  adot(3,3)                       !  metric in direct space
  real(REAL64)                       ::  rk0(3)                          !  reference k-point (in lattice coordinates)

! allocatable arrays

  complex(REAL64), allocatable       ::  h0(:,:)                         !  <Psi|H|Psi>
  complex(REAL64), allocatable       ::  dh0drk(:,:,:)                   !  d <Psi|H|Psi> d k
  complex(REAL64), allocatable       ::  d2h0drk2(:,:,:,:)               !  d^2 <Psi|H|Psi> d k^2

  real(REAL64), allocatable          ::  ei(:)                           !  eigenvalue no. i. (hartree)

! allocatable arrays for Brillouin zone path

  integer                            ::  nlines                          !  number of lines in reciprocal space
  integer, allocatable               ::  nkstep(:)                       !  number of steps in line
  logical, allocatable               ::  ljump(:)                        !  indicates if the new line contains a jump from the preceeding
  integer                            ::  nvert                           !  number of vertical lines in plot
  real(REAL64), allocatable          ::  xcvert(:)                       !  x coordinate of vertical line
  real(REAL64), allocatable          ::  xk(:)                           !  x coordinate of k-point in plot
  real(REAL64), allocatable          ::  rk(:,:)                         !  x coordinate of k-point in plot
  real(REAL64), allocatable          ::  e_of_k(:,:)                     !  band energies of k-point in plot
  character(len=6), allocatable      ::  label(:)                        !  label of symmetry k-points
  real(REAL64), allocatable          ::  xklab(:)                        !  x coordinate of label

! local variables

  integer           ::  nrk2
  real(REAL64)      ::  rkpt(3)             !  j-th component in lattice coordinates of the k-point
  real(REAL64)      ::  eref                !  reference energy for plot
  integer           ::  nocc                !  number of occupied states (different color)
  integer           ::  nstyle              !  choice of plot style
  integer           ::  ichoice

  integer           ::  neltot              !  total number of electrons
  real(REAL64)      ::  ztot                !  total number of electrons
  real(REAL64)      ::  efermi              !  Fermi energy (not used)
  integer           ::  ispin               !  spin degeneracy

  integer           ::  neigso              !  number of eigenvalues with spin-orbit

!  character(len=1)  ::  yesno

! counters

  integer    ::  irk, j


! reads silvaco style file

  write(6,*)
  write(6,*) '  Which file do you want to use? (1,2,3)'
  write(6,*)
  write(6,*) '  1) matrix_kdotp_so.dat  (default with spin-orbit)'
  write(6,*) '  2) matrix_kdotp.dat  (default without spin-orbit)'
  write(6,*) '  3) Choose file name.'
  write(6,*)

  read(5,*) ichoice
  write(ioreplay,*) ichoice,'   default file choice'

  if(ichoice == 1) then
    filename = 'matrix_kdotp_so.dat'
  elseif(ichoice == 2) then
    filename = 'matrix_kdotp.dat'
  elseif(ichoice == 3) then
    write(6,*)
    write(6,*) ' Enter file name'
    write(6,*)
    read(5,*) filename
    write(ioreplay,*) filename,'   filename'
  else
    filename = 'matrix_kdotp_so.dat'
    write(6,*)
    write(6,*) ' Wrong answer, using  matrix_kdotp_so.dat'
    write(6,*)
  endif

  efermi = 0.0

  iotape = 23

  call kdotp_in_size(trim(filename), iotape, mxdbnd)

  allocate(h0(mxdbnd,mxdbnd))
  allocate(dh0drk(mxdbnd,mxdbnd,3))
  allocate(d2h0drk2(mxdbnd,mxdbnd,3,3))

  call kdotp_in(trim(filename), iotape, neig, adot, rk0, neltot, ispin,  &
      h0, dh0drk, d2h0drk2,                                              &
      mxdbnd)

  ztot = neltot
  iotape = 13
  call out_band_circuit_size('BAND_LINES.DAT', iotape, 1, adot, ztot,    &
                   neig, nrk2, nlines, nvert)

  allocate(xk(nrk2))
  allocate(rk(3,nrk2))
  allocate(xcvert(nvert))
  allocate(ljump(nlines))
  allocate(nkstep(nlines))
  allocate(label(nvert+nlines))
  allocate(xklab(nvert+nlines))

  call out_band_get_circuit('BAND_LINES.DAT', iotape, 1, adot,           &
                   xk, rk, xcvert, ljump, nkstep, label, xklab,          &
                   neig, nrk2, nlines, nvert)

  neigso = (2*neig)/ispin

  if(neigso > mxdbnd) then
    neigso = mxdbnd

    write(6,*)
    write(6,*) '  neig reduced to: ', neigso
    write(6,*)

  endif

  allocate(e_of_k(neigso,nrk2))

! finds mxddim, mxdbnd

  allocate(ei(mxdbnd))

  do irk=1,nrk2

!   loop over k-points

    do j=1,3
      rkpt(j) = rk(j,irk)
    enddo


    call kdotp_diag_nopsi(rkpt, mxdbnd, ei,                              &
          rk0, h0, dh0drk, d2h0drk2,                                     &
          mxdbnd)

    do j=1,neigso
      e_of_k(j,irk) = ei(j)
    enddo

  enddo

  nstyle = 2

  call out_band_eref(neig, nrk2, rk, ztot, efermi, ispin, 1, e_of_k, eref, nocc)

  call out_band_gnuplot('band_kp.gp', iotape,                            &
         neigso, nrk2, xk, e_of_k, eref,                                 &
         nvert, xcvert, nlines, ljump, nkstep, label, xklab)

  call out_band_xmgrace('band_kp.agr', iotape,                           &
         title, subtitle, nstyle,                                        &
         neigso, nrk2, xk, e_of_k, eref, nocc,                           &
         nvert, xcvert, nlines, ljump, nkstep, label, xklab)

  return

end subroutine kdotp_sub

