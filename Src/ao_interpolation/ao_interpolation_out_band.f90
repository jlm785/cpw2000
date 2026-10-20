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

!>  This subroutine calculates the band structure along a path
!>  using atomic orbital interpolation
!>
!>  \author       Carlos Loia Reis
!>  \version      5.13
!>  \date         before 2017. 20 October 2026.
!>  \copyright    GNU Public License v2

  subroutine ao_interpolation_out_band(title, subtitle,                  &
         noiData, ztot, adot)

! Written by Carlos Loia Reis before 2017, based on previous code,
! Modified, documentation, May 2020. JLM
! Modified, ztot in out_band_circuit_size, indentation. 26 July 2024. JLM
! Bug when neig greater than nband. 1 November 2024. JLM
! Increase dimension of label. 17 September 2025. JLM
! Constants updated to CODATA 2022 (HARTREE). 3 October 2026. JLM+claude
! Progress of the k-points with progress_line. 4 October 2026. JLM+claude
! Number of bands with spin-orbit limited to the spinor bands available. 19 October 2026. JLM+claude
! Old style BAND_MESH.DAT calculation removed, with its variables. 20 October 2026. JLM+claude


  use NonOrthoInterp

  implicit none

  integer, parameter                 :: REAL64 = selected_real_kind(12)

! input

  character(len=50), intent(in)      ::  title                           !<  title for plots
  character(len=140), intent(in)     ::  subtitle                        !<  subtitle for plots

  type(noiData_t)                    ::  noiData                         !<  see NonOrthoInterp

  real(REAL64), intent(in)           ::  ztot                            !<  total charge density (electrons/cell)
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space

! allocatable arrays for Brillouin zone path

  integer                            ::  nlines                          !  number of lines in reciprocal space
  integer, allocatable               ::  nkstep(:)                       !  number of steps in line
  logical, allocatable               ::  ljump(:)                        !  indicates if the new line contains a jump from the preceeding
  integer                            ::  nvert                           !  number of vertical lines in plot
  real(REAL64), allocatable          ::  xcvert(:)                       !  x coordinate of vertical line
  real(REAL64), allocatable          ::  xk(:)                           !  x coordinate of k-point in plot
  real(REAL64), allocatable          ::  rk(:,:)                         !  x coordinate of k-point in plot
  real(REAL64), allocatable          ::  e_of_k(:,:)                     !  band energies of k-point in plot
  real(REAL64), allocatable          ::  e_of_k_so(:,:)                  !  spin-orbit band energies of k-point in plot
  character(len=10), allocatable     ::  label(:)                        !  label of symmetry k-points
  real(REAL64), allocatable          ::  xklab(:)                        !  x coordinate of label

! local variables


  integer           ::  neig           !  number of eigenvectors required (maybe modified on output)
  real(REAL64)      ::  rkpt(3)        !  j-th component in lattice coordinates of the k-point

  real(REAL64)      ::  eref           !  reference energy for plot
  integer           ::  nocc           !  number of occupied states (different color)
  integer           ::  nstyle         !  choice of plot style

  integer           ::  irk
  integer           ::  iotape
  integer           ::  nrk2

! counters

  integer    ::  j, n

  real(REAL64), allocatable :: ev_interp(:)


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

! checks maximum of number of bands (with spin-orbit nband counts the spinor bands)

  if(noiData%lso == 1) then
    if(2*neig > noiData%nband) neig = noiData%nband/2
  else
    if(neig > noiData%nband) neig = noiData%nband
  endif

  allocate(e_of_k(neig,nrk2))
  allocate(e_of_k_so(2*neig,nrk2))

!  mxdbnd = neig

  allocate(ev_interp(noiData%nband))

  do irk=1,nrk2

!   loop over k-points

    do j=1,3
      rkpt(j) = rk(j,irk)
    enddo

    call progress_line(irk, nrk2, 'k-points:', 6)

    call NonOrthoInterpRun(noiData,rkpt,ev_interp)

    if(noiData%lso==1) then
     do j=1,2*neig
       e_of_k_so(j,irk) = ev_interp(j)
     enddo
    else
     do j=1,neig
       e_of_k(j,irk) = ev_interp(j)
     enddo

    endif

  enddo

  iotape = 15
  nstyle = 1

! writes the output files for xmgrace and gnuplot

  if(noiData%lso==1) then

     n = min(nint(ztot + 0.01),2*neig)
     eref = e_of_k_so(n,1)
     do irk = 1,nrk2
     do j=1,n
       if(e_of_k_so(j,irk) > eref) eref = e_of_k_so(j,irk)
     enddo
     enddo

     nocc = n

    call out_band_xmgrace('band_interp_so.agr', iotape,                  &
        title, subtitle, nstyle,                                         &
        2*neig, nrk2, xk, e_of_k_so, eref, nocc,                         &
        nvert, xcvert, nlines, ljump, nkstep, label, xklab)

    call out_band_dots_xmgrace('dots_band_interp_so.agr', iotape,        &
        nint(ztot), 2*neig, nrk2, xk, e_of_k_so,                         &
        eref, nvert, xcvert, nlines, ljump, nkstep, label, xklab)

  else

     n = min(nint(0.5*ztot + 0.01),neig)
     eref = e_of_k(n,1)
     do irk = 1,nrk2
     do j=1,n
       if(e_of_k(j,irk) > eref) eref = e_of_k(j,irk)
     enddo
     enddo

     nocc = n

     call out_band_xmgrace('band_interp.agr', iotape,                    &
         title, subtitle, nstyle,                                        &
         neig, nrk2, xk, e_of_k, eref, nocc,                             &
         nvert, xcvert, nlines, ljump, nkstep, label, xklab)

     call out_band_dots_xmgrace('dots_band_interp.agr', iotape,          &
         nint(ztot)/2, neig, nrk2, xk, e_of_k,                           &
         eref, nvert, xcvert, nlines, ljump, nkstep, label, xklab)

  endif


! Cleanup

  deallocate(xk)
  deallocate(rk)
  deallocate(xcvert)
  deallocate(ljump)
  deallocate(nkstep)
  deallocate(label)
  deallocate(xklab)

  deallocate(e_of_k)
  deallocate(e_of_k_so)

  deallocate(ev_interp)

  return

end subroutine ao_interpolation_out_band
