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

!>generates the data for a 3D contourplot using xcrysden or vesta

subroutine plot_contour3D(ioreplay, func,                                &
     adot, ntype, natom, nameat, rat,                                    &
     ng, kgv,                                                            &
     mxdtyp, mxdatm, mxdgve)

! rewriten february 24,1995.jlm
! modified 28 april 2004
! modified, f90, subroutine, 27 May 2014. JLM
! Written 30 May 2014, from 2D and LightSO code. JLM
! Modified, documentation, 11 June 2020. JLM
! Modified to use calls to write_xsf, etc.., 1-4 February 2021. JLM
! copyright  Jose Luis Martins, Carlos Loia Reis/INESC-MN

! version 4.99


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                :: ioreplay                         !<  tape number for reproducing calculations

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors

  complex(REAL64), intent(in)        ::  func(mxdgve)                    !<  quantity to be plotted for g-vector j

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space
  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i
  character(len=2), intent(in )      ::  nameat(mxdtyp)                  !<  chemical symbol for the type i
  real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)            !<  k-th component (in lattice coordinates) of the position of the n-th atom of type i

  integer, intent(in)                ::  ng                              !<  total number of g-vectors with length less than gmax
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length

! allocatable arrays

  complex(REAL64), allocatable       :: chd(:,:,:)                       !  function in the fft grid
  complex(REAL64), allocatable       :: wrk(:)                           !  work array.

! other variables

  integer      ::  id, n1,n2,n3
  integer      ::  k1,k2,k3, kd
  integer      ::  mxdwrk

  integer      ::  kmscr(3), nsfft(3)
  integer      ::  mfft, mwrk
  integer      ::  ktmp(3)
  
  integer      ::  ierr
  integer      ::  iotape

  real(REAL64) ::  dmin, dmax, cmax, small, abschd
  
  real(REAL64) ::  avec(3,3), bvec(3,3)
  logical                ::  lvesta
  
  character(len=1)  ::  yesno

! constants

  real(REAL64), parameter ::  ZERO = 0.0_REAL64
  real(REAL64), parameter ::  BOHR = 0.5291772109_REAL64

! counters

  integer      ::  i, j, k


  call adot_to_avec_sym(adot,avec,bvec)

! finds basic fft grid 

  do j=1,3
    kmscr(j) = 0
  enddo

  do i=1,ng
    do j=1,3
      if(abs(kgv(j,i)) > kmscr(j)) kmscr(j) = kgv(j,i)
    enddo
  enddo

  call size_fft(kmscr, nsfft, mfft, mwrk)

  write(6,*) 
  write(6,'("  The basic fft grid is: ",3i8)') (nsfft(j),j=1,3)
  write(6,*)
  write(6,*) '  For a final smooth plot you may need to double' 
  write(6,*) '  these values.' 
  write(6,*)
  write(6,*) '  Do you want a finer grid? (y/n)'
  write(6,*)
  
  read(5,*) yesno
  write(ioreplay,'(2x,a1,"   new fft grid")') yesno

  if(yesno == 'y' .or. yesno == 'Y') then

    write(6,*)
    write(6,*) '  Enter the new grid size (n1,n2,n3)'
    write(6,*) '  These values will be overridden for fast FFT!'
    write(6,*)
    
    read(5,*) (ktmp(j),j=1,3)
    write(ioreplay,'(3(2x,i8),"   fft grid")') (ktmp(j),j=1,3)

    do j=1,3
      ktmp(j) = (ktmp(j)-1)/2
      kmscr(j) = max(kmscr(j),ktmp(j))
    enddo

    call size_fft(kmscr, nsfft, mfft, mwrk)

    write(6,*) 
    write(6,'("  The fft grid will be: ",3i8)') (nsfft(j),j=1,3)
    write(6,*)

  endif

! initialize the array

  write(6,*)
  write(6,'("  The integral of the function is: ",g12.4)') real(func(1))
  write(6,*)


  n1 = nsfft(1)
  n2 = nsfft(2)
  n3 = nsfft(3)

  id = n1 + 1

! initialize charge density array and enter symmetrized
! charge.

  allocate(chd(id,n2,n3))

  do k=1,n3
  do j=1,n2
  do i=1,id
    chd(i,j,k) = cmplx(ZERO,ZERO,REAL64)
  enddo
  enddo
  enddo

! wraparound for charge (should not occur here )

  do i=1,ng
    k1 = kgv(1,i)
    kd = n1*(k1/n1)
    if (k1 < 0) kd = kd - n1
    k1 = k1 - kd + 1
    k2 = kgv(2,i)
    kd = n2*(k2/n2)
    if (k2 < 0) kd = kd - n2
    k2 = k2 - kd + 1
    k3 = kgv(3,i)
    kd = n3*(k3/n3)
    if (k3 < 0) kd = kd - n3
    k3 = k3 - kd + 1

    chd(k1,k2,k3) = chd(k1,k2,k3) + func(i)
  enddo

  mxdwrk = 2*max(n1*n2,n1*n3,n2*n3)

  allocate(wrk(mxdwrk))

! fourier transform to real space

  call cfft_c16(chd, id, n1,n2,n3, -1, wrk, mxdwrk)

  dmax = real(chd(1,1,1))
  dmin = dmax
  cmax = ZERO
  ierr = 0
  small=1.0d-8
  do k=1,n3
  do j=1,n2
  do i=1,n1
    if (real(chd(i,j,k)) > dmax) dmax = real(chd(i,j,k))
    if (real(chd(i,j,k)) < dmin) dmin = real(chd(i,j,k))
    abschd = abs(aimag(chd(i,j,k)))
    if (abschd > cmax) cmax = abschd
    if (abschd > small) ierr = ierr+1
    chd(i,j,k) = cmplx(real(chd(i,j,k)),ZERO,REAL64)
  enddo
  enddo
  enddo

  write(6,'("  max and min of function ",3f12.4)') dmax,dmin,cmax
  if(ierr > 0) then

    write(6,'("  stopped:  complex function ")')

    stop

  endif
  
  iotape = 47
  lvesta = .TRUE.

  open(unit=iotape,file='rho3D.xsf',form='formatted')

! writes crystal structure information

  call plot_xsf_crys(iotape, lvesta,                                &
      adot, ntype, natom, nameat, rat,                              &
      mxdtyp, mxdatm)

! writes the function on the grid

  call plot_xsf_data(iotape, lvesta, adot, chd, id, n1,n2,n3)


  close(unit=iotape)
  
  write(6,*)
  write(6,*) 'The file rho3D.xsf was written'
  write(6,*)
  write(6,*) 'You can see the plot with VESTA or xcrysden'
  write(6,*)

  deallocate(chd)
  deallocate(wrk)

  return
end subroutine plot_contour3D
