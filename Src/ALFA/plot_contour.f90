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

!>  Generates a contourplot of a function
!>  in a plane determined by one corner c0(i) and two vectors
!>  vx(i) and vy(i) giving two sides. these are all given in units
!>  of the lattice basis vectors.

subroutine plot_contour(ioreplay, func, adot, ng, kgv, mxdgve)

! rewriten february 24,1995.jlm
! modified 28 april 2004
! modified, f90, subroutine, 27 May 2014. JLM
! Documentation, merge rho/psi. 4 February 2021. JLM
! copyright  Jose Luis Martins/INESC-MN

! version 4.99


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                :: ioreplay                         !<  tape number for reproducing calculations

  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors

  complex(REAL64), intent(in)        ::  func(mxdgve)                    !<  function for g-vector j

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space

  integer, intent(in)                ::  ng                              !<  total number of g-vectors with length less than gmax
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length

! allocatable arrays

  complex(REAL64), allocatable       :: rhopl(:,:)                       !  function on planar grid
  real(REAL64), allocatable          :: ro(:,:)                          !  function on planar grid
  complex(REAL64), allocatable       :: chd(:,:,:)                       !  function in the fft grid
  complex(REAL64), allocatable       :: wrk(:)                           !  work array.

! other variables

  integer      :: id, n1,n2,n3, k1,k2,k3, kd
  integer      :: mxdwrk

  integer      ::  kmscr(3), nsfft(3)
  integer      ::  mfft, mwrk
  integer      ::  ktmp(3)
  
  integer      ::  nx,ny
  integer      ::  ierr
  
  integer      ::  io

  real(REAL64) ::  dmin,dmax,cmax,small,abschd
  
  character(len=1)  ::  yesno
  real(REAL64)         ::  c0(3)                       !  corner (origin) of plane
  real(REAL64)         ::  dx(3),dy(3)                 !  step vectors that define the plane (lattice coordinates)
  real(REAL64)         ::  xscale,yscale               !  aspect ratio of plot

! constants

  real(REAL64), parameter :: ZERO = 0.0_REAL64

! counters

  integer      ::  i, j, k

! finds basic fft grid 

  do j = 1,3
    kmscr(j) = 0
  enddo

  do i=1,ng
    do j=1,3
      if(abs(kgv(j,i)) > kmscr(j)) kmscr(j) = kgv(j,i)
    enddo
  enddo

  call size_fft(kmscr,nsfft,mfft,mwrk)

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

    do j = 1,3
      ktmp(j) = (ktmp(j)-1)/2
      kmscr(j) = max(kmscr(j),ktmp(j))
    enddo

    call size_fft(kmscr,nsfft,mfft,mwrk)

    write(6,*) 
    write(6,'("  The fft grid will be: ",3i8)') (nsfft(j),j=1,3)
    write(6,*)

  endif

! initialize function array 

  write(6,*)
  write(6,'("  The integral of the function is: ",g12.4)') real(func(1))
  write(6,*)


  n1 = nsfft(1)
  n2 = nsfft(2)
  n3 = nsfft(3)

  id = n1 + 1

  allocate(chd(id,n2,n3))

  do k=1,n3
  do j=1,n2
  do i=1,id
    chd(i,j,k) = cmplx(ZERO,ZERO,REAL64)
  enddo
  enddo
  enddo

! wraparound for function (should not occur here )

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

  deallocate(wrk)

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


! gets the details of plot

  write(6,*)
  write(6,*) '  The program will now generate a contour plot '
  write(6,*) '  of the function '
  write(6,*)
  write(6,*)
  write(6,*) '  enter number of points in surface nx,ny '
  write(6,*) '  (suggestion:  101 101 for nice plot) '
  write(6,*)

  read(5,*) nx,ny
  write(ioreplay,'(2(2x,i8),"   2D plot grid")') nx,ny

  allocate(rhopl(nx,ny))
  allocate(ro(nx,ny))

  call plot_get_plane(ioreplay, adot, nx,ny, c0, dx,dy, xscale,yscale)

! quadratic interpolation in the fft mesh

  call plot_mat_fft(rhopl, nx,ny, chd, id, n1,n2,n3, c0,dx,dy)

  do j=1,ny
  do i=1,nx
    ro(i,j) = real(rhopl(i,j))
  enddo
  enddo

  deallocate(rhopl)

! plots

  write(6,*)
  write(6,*) '  Do you want to see a rough ASCII plot? (y/n)'
  write(6,*)
  
  read(5,*) yesno
  write(ioreplay,'(2x,a1,"    ASCII plot")') yesno

  if(yesno == 'y' .or. yesno == 'Y') then

    call plot_alf_plt(ro, nx,ny)

  endif

  io = 15
  if(io == ioreplay) io = 16

  call plot_2D_gnuplot(ioreplay, 'rho_2D.gp', io, ro, nx,ny, xscale,yscale)

  deallocate(ro)
  deallocate(chd)

  return

  end subroutine plot_contour
