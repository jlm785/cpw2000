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

!>  Estimates the radiative recombination rate from the
!>  optical functions.  master subroutine

subroutine opt_rad_sub(ioreplay)

! Written 15-30 October 2020. JLM
! based on optical subroutine by Carlos Loia Reis. July 2020
! copyright  Carlos Loia Reis/Jose Luis Martins/INESC-MN

! version 4.98

  implicit none

  integer, parameter          ::  REAL64 = selected_real_kind(12)
  integer, parameter          ::  REAL32 = selected_real_kind(6)

! input

  integer, intent(in)                 ::  ioreplay                       !<  tape number for reproducing calculations

! local allocatable variables
  
  real(REAL64), allocatable           ::  el(:,:)                        !  eigenvalues in Hartree in irreducible wedge

  integer, allocatable                ::  kmap(:,:,:,:)                  !  kmap(1,...) corresponding k-point, kmap(2,...) = 1 additional inversion, kmap(3,...) symmetry operation
  integer, allocatable                ::  indk(:,:)                      !  index of the six k-points neighbouring k-point i
  real(REAL64),allocatable            ::  rk(:,:)                        !  component in lattice coordinates of the k-point in the mesh
  real(REAL64),allocatable            ::  wgk(:)                         !  weight in the integration of k-point

  integer, allocatable                ::  nband(:)                       !  number of bands for each k

  real(REAL64), allocatable           ::  ev_grid(:,:,:,:)               !  eigenvalue in grid

  real(REAL64), allocatable           ::  ehist(:)                       !  histogram energy (DOS)
  real(REAL64), allocatable           ::  dhist(:)                       !  histogram dos
  real(REAL64), allocatable           ::  chist(:)                       !  histogram intergrated dos)

  real(REAL64), allocatable           ::  e_re(:), e_im(:)               !  real and imaginary parts of the dielectric function

  real(REAL64), allocatable        ::  e_of_k(:,:)                       !  band energies of k-point in plot
  real(REAL64), allocatable        ::  e_of_k_so(:,:)                    !  spin-orbit band energies of k-point in plot

  real(REAL64), allocatable           ::  ehopt(:)                       !  histogram energy (Optical)

  real(REAL64), allocatable           ::  ehrad(:)                       !  histogram energy (Optical)
  real(REAL64), allocatable           ::  e_rad_re(:), e_rad_im(:)       !  real and imaginary parts of the dielectric function
  real(REAL64), allocatable           ::  rsrad(:)                       !  function t be integrated to give radiative coefficient

! local variables read from file

  integer                             ::  mxdbnd                         !  array dimension for the number of bands

  integer                             ::  nx, ny ,nz                     !  grid size
  integer                             ::  nkgrid(3)                      !  nkgrid(1) = nx,...

  integer                             ::  ntrans                         !  number of symmetry operations in the factor group
  integer                             ::  mtrx(3,3,48)                   !  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the factor group

  integer                             ::  nrk                            !  number of k-points for integration in the irreducible wedge of the brillouin zone

  real(REAL64)                        ::  ztot                           !  total number of electrons
  integer                             ::  neig                           !  number of bands
  integer                             ::  nval                           !  number of valence bands
  integer                             ::  ncond                          !  number of conduction bands
  integer                             ::  nvtc                           !  nval*ncond or...
  
  real(REAL64)                        ::  adot(3,3)                      !  metric in real space

  integer                             :: ispin

  character(len=50)                   ::  title                          !  title for plots
  character(len=140)                  ::  subtitle                       !  subtitle for plots

  logical               ::  lscl                                         !  write for scalar (no-spin-orbit)
  logical               ::  lso                                          !  also write for spin-orbit

  character(len = 1) ::  yesno
  integer            ::  io, iotape

! other local variables

  integer            ::  nhtarg                                          !  target number of points in  histogram

  integer            ::  nvbm, ncbm                                      !  approximate indices for conduction and valence band
  real(REAL64)       ::  evbb                                            !  energy of bottom of valence band
  real(REAL64)       ::  evbm                                            !  energy of maximum of valence band
  real(REAL64)       ::  ecbm                                            !  energy of minimum of conduction band
  logical            ::  linsul                                          !  true if it appears to be an insulator
  real(REAL64)       ::  ezero                                           !  zero of energy
  real(REAL64)       ::  egap                                            !  band gap (underestimate)

  integer            ::  ninter                                          !  number of interpolated points for big grid (~4)
  integer            ::  nint_rad                                         !  number of interpolated points for optical grid (~4-16)
  integer            ::  nhist                                           !  number of points in histogram
  real(REAL64)       ::  emin                                            !  energy minimum in histogram
  real(REAL64)       ::  deltae                                          !  histogram step
  logical            ::  lunif                                           !  true if uniform histogram

  real(REAL64)       ::  tempk                                           !  temperature in K
  real(REAL64)       ::  tau                                             !  temperature in au

  real(REAL64)       ::  xnp                                             !  n*p * exp(E_gap/k_B T)
  real(REAL64)       ::  xni, xnisi                                      !  intrincic carrier concentration

  real(REAL64)       ::  xmu                                             !  electron chemical potential, AKA Fermi level
  real(REAL64)       ::  xmdosv, xmdosc                                  !  dos effective masses for valence and conduction
  real(REAL64)       ::  xnimass                                         !  intrinsic carrier

  logical            ::  lidos                                           !  true if integrated density of states is to be computed.
  real(REAL64)       ::  vcell                                           ! cell volume
  real(REAL64)       ::  bdot(3,3)                                       ! metric in reciprocal space (contravariant components)

  real(REAL64)       ::  egapopt                                         !  optical gap
  real(REAL64)       ::  pcvsq                                           !  pcv^2

  integer            ::  nhopt                                           !  number of histogram points for optical response
  real(REAL64)       ::  ezopt                                           !  zero of energy for optical response

  
  integer            ::  nhrad                                           !  number of histogram points for radiative coefficient
  integer            ::  nhradmin
  integer            ::  nrange                                          !  number of band pairs within erange

  integer            ::  nordp1                                          !  order of Lagrange interpolation
  real(REAL64)       ::  dymax                                           !  estimate of error of Lagrange interpolation

  real(REAL64)       ::  xlth                                            !  thermal wavelength
  real(REAL64)       ::  pref                                            !  prefactor
  real(REAL64)       ::  bcoef                                           !  coefficient B of radiative recombination rate
  real(REAL64)       ::  rsint                                           !  Integral for Rs
  real(REAL64)       ::  abscoef                                         !  absorption coefficient
  real(REAL64)       ::  e_mod                                           !  epsilon**2
  real(REAL64)       ::  xn_im                                           !  imaginary part (extiction coefficient) of the refractive index
  real(REAL64)       ::  xn_re                                           !  refractive index (real part)
  
  character(len=16)                   ::  filename                       !  input file name

  integer                            ::  identif                         !  identifier that is almost 

  character(len=40)  ::  filename_out
  character(len=50)  ::  func                                            !  name of function

  integer                            ::  io_unsym
  integer                            ::  io_dhdrk
  integer                            ::  io_tmp
  integer                            ::  io_grid

  character(len=16)                  ::  filedhdrk 
  character(len=16)                  ::  fileunsym 
  character(len=16)                  ::  filetmp 
  character(len=16)                  ::  filegrid 

! constants

  real(REAL64), parameter ::  ZERO = 0.0_REAL64 , UM = 1.0_REAL64
  real(REAL64), parameter ::  PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter ::  EPS = 1.0E-10_REAL64

  real(REAL64), parameter ::  HARTREE = 27.211386246_REAL64
  real(REAL64), parameter ::  TAUTOK = 11604.9_REAL64 * HARTREE
  real(REAL64), parameter ::  BOHR = 0.5291772109E-10_REAL64 
  real(REAL64), parameter ::  ALPHA = UM / 137.036
  real(REAL64), parameter ::  AUT = 2.4188843266E-17_REAL64
  real(REAL64), parameter ::  ELMASS = 9.10938356E-31_REAL64

! counters

  integer            :: irk, i


! parameters for the calculation

  filename = "dos_file.dat"
  filetmp = "tmp_tmp.dat"
  fileunsym = "tmp_unsym.dat"
  filegrid = "tmp_grid.dat"

  io = 21
  iotape = 15

  io_dhdrk = 80
  io_unsym = 81
  io_tmp = 82
  io_grid = 83

! these parameters are hardcoded.  Change them if you need more precision

  nhtarg = 2000
  nordp1 = 4

  write(6,*)
  write(6,*) '   ESTIMATE of rhe radiative recombination rate'
  write(6,*) 

  write(6,*) '   Enter temperature in K'
  write(6,*) 

  read(5,*) tempk
  write(ioreplay,*) tempk,'     temperature in K'
  write(6,*) 

  if(tempk < 10.0 .or. tempk > 10000.0) then
    write(6,*) '   strange temperature, expect problems....'
    write(6,*) 
  endif

  tau = tempk / TAUTOK

  write(6,*) '   do you want a calculation including Spin-Orbit (y/n) ?'
  write(6,*) 

  read(5,*) yesno
  write(ioreplay,*) yesno,'     spin-orbit for IE'

  if(yesno == 'y' .or. yesno == 'Y') then
    ispin = 1
  else
    ispin = 2
  endif


  write(6,*)
  write(6,*) '   How many interpolation points do you want for optical response (~4)?'
  write(6,*) '   Enter 1 for no interpolation.'
  write(6,*) '   Do not use interpolation for perfect superlattices.'
  read(5,*) ninter
  write(ioreplay,*) ninter,'     number of superquadratic interpolation points for optical'

  if(ninter < 1) then
    ninter = 1
    write(6,*)
    write(6,*) '   wrong value, ninter set to ', ninter
    write(6,*)
  endif
 
  if(ninter > 10) then
    ninter = 10
    write(6,*)
    write(6,*) '   too large value, ninter set to ', ninter
    write(6,*)
  endif


  write(6,*)
  write(6,*) '   How many interpolation points do you want for radiative rates (~4 to 16)?'
  write(6,*) '   Enter 1 for no interpolation.'
  write(6,*) '   Do not use interpolation for perfect superlattices.'
  read(5,*) nint_rad
  write(ioreplay,*) nint_rad,'     number of superquadratic interpolation points'

  if(nint_rad < 1) then
    nint_rad = 1
    write(6,*)
    write(6,*) '   wrong value, nint_rad set to ', nint_rad
    write(6,*)
  endif
 
  if(nint_rad > 32) then
    nint_rad = 32
    write(6,*)
    write(6,*) '   too large value, nint_rad set to ', nint_rad
    write(6,*)
  endif

! reads data file

  do i = 1,140
    subtitle(i:i) = ' '
  enddo
  do i = 1,50
    title(i:i) = ' '
  enddo

! reads data file


  call dos_read_size(trim(filename), io, nrk, mxdbnd, nkgrid)

  nx = nkgrid(1)
  ny = nkgrid(2)
  nz = nkgrid(3)

  allocate(e_of_k(mxdbnd,nrk))
  allocate(e_of_k_so(2*mxdbnd,nrk))
  allocate(nband(nrk))
  allocate(wgk(nrk))
  allocate(rk(3,nrk))
  allocate(indk(6,nrk))
  allocate(kmap(3,nx,ny,nz))

  call dos_read_data(trim(filename), io, title, subtitle,                &
    lscl, lso, identif,                                                  &
    nrk, nx, ny, nz, ztot, adot, ntrans, mtrx,                           &
    nband, rk, wgk, indk, kmap, e_of_k, e_of_k_so,                       &
    nrk, mxdbnd)

  call adot_to_bdot(adot,vcell,bdot)

  neig = nband(1)
  do irk = 1,nrk
    if(neig > nband(irk)) neig = nband(irk)
  enddo

  write(6,*)
  write(6,*) '   The title and subtitle associated with the data file are:' 
  write(6,*)
  write(6,*) '   ',title
  write(6,*) '   ',subtitle
  write(6,*)


  if(ispin == 1) then
    if(.not. lso) then
      write(6,*) '    returning in opt_sub:  file does not contain spin-orbit info'

      return

    endif
    mxdbnd = 2*mxdbnd
    neig = 2*neig
    nval = nint(ztot)
    ncond = neig - nval
    allocate(el(mxdbnd,nrk))
    el(:,:) = e_of_k_so(:,:)
    do irk = 1,nrk
      nband(irk) = 2*nband(irk)
    enddo

    filedhdrk = 'opt_dhdrk_so.dat'

    write(6,*)
    write(6,*) '     Optical functions calculated with spin-orbit'
    write(6,*)

  else
    if(.not. lscl) then
      write(6,*) '    returning in dos_sub:  file does not contain scalar info'

      return

    endif

    nval = nint(0.5*ztot)
    ncond = neig - nval

    allocate(el(mxdbnd,nrk))
    el(:,:) = e_of_k(:,:)

    filedhdrk = 'opt_dhdrk.dat'

    write(6,*)
    write(6,*) '     Optical functions calculated WITHOUT spin-orbit'
    write(6,*)

  endif

  deallocate(e_of_k)
  deallocate(e_of_k_so)
  deallocate(rk)
  deallocate(indk)
  deallocate(wgk)

! generates unsymmetrized grid

  call opt_grid_unsym(el, filedhdrk, io_dhdrk, fileunsym, io_unsym, kmap, mtrx, &
    neig, nx,ny,nz, nrk, mxdbnd)


  allocate(ev_grid(nx,ny,nz,mxdbnd))

  call dos_el_to_egrid(el, nrk, nband, nkgrid, kmap, ev_grid, nkgrid, mxdbnd)

  call dos_bands_rough(el, nrk, nband, ztot, ispin, evbb, evbm, ecbm, linsul, mxdbnd)

  deallocate(kmap)
  deallocate(nband)
  
  ezero = evbm

  if(linsul) then

    deltae = 0.001_REAL64/HARTREE
    emin = evbm - ezero - 0.5/HARTREE
    emin = deltae*real(int(emin/deltae))
    nhist = int((ecbm-evbm + 1.0/HARTREE) / deltae)

  else

    write(6,*) '  The system is not an insulator.  Rewrite the code for metals'

    stop

  endif

  allocate(ehist(nhist))
  allocate(chist(nhist))
  allocate(dhist(nhist))

  do i=1,nhist
    ehist(i) = emin + deltae*(i-1)
    dhist(i) = ZERO
    chist(i) = ZERO
  enddo


  lunif = .TRUE.
  lidos = .FALSE.

  call dos_quad(nkgrid, ev_grid,ninter, nhist, ehist, dhist, chist, ezero, ispin, &
    lunif, .TRUE., lidos, mxdbnd)

  nvtc = nval*ncond


  call opt_set_in_grid(0,0, adot, fileunsym, io_unsym, filetmp, io_tmp,       &
                       filegrid, io_grid, neig, nval, ncond, nx, ny, nz)

  call execute_command_line("rm " // fileunsym // " 2> /dev/null ")
  call execute_command_line("rm " // filetmp // " 2> /dev/null ")

! Calculate the optical response in a broad range

  ezopt = ZERO


  call opt_plot_range(el, neig, nval, nhtarg, deltae, nhopt, nrk, mxdbnd)

  deallocate(el)

  allocate(ehopt(nhopt))
  allocate(e_im(nhopt))
  allocate(e_re(nhopt))

  do i = 1,nhopt
    ehopt(i) = deltae*i
!    ehopt(i) = UM/(10*HARTREE) + deltae*(i-1)
  enddo

  call opt_calc_quad(ispin, adot, filegrid, io_grid, ninter, ehopt,    &
                         e_re, e_im, nx,ny,nz, nhopt, nvtc)

  call opt_write(0,0, iotape, title, subtitle,                    &
                     nhopt, ehopt, e_re, e_im, vcell, ztot, ispin)

!   filename_out = 'optical_epsilon_im.dat'
!   func = 'Imaginary Part of Dielectric Function             '
! 
!   call opt_write_file(trim(filename_out), iotape, title, subtitle, func, &
!                       nhopt, ehopt, e_im, vcell, ztot)
! 
!   filename_out = 'optical_epsilon_re.dat'
!   func = 'Real Part of Dielectric Function                  '
! 
!   call opt_write_file(trim(filename_out), iotape, title, subtitle, func, &
!                       nhopt, ehopt, e_re, vcell, ztot)

! range for radiative coefficient integration

  call opt_rad_plot_range(filegrid, io_grid, nhtarg, 0.01D0, nrange,     &
          emin, deltae, nhrad, nx,ny,nz, nvtc)

  allocate(ehrad(nhrad))
  allocate(e_rad_im(nhrad))
  allocate(e_rad_re(nhrad))
  allocate(rsrad(nhrad))

  do i=1,nhrad
    ehrad(i) = emin + deltae*(i-1)
    e_rad_im(i) = ZERO
    e_rad_re(i) = ZERO
    rsrad(i) = ZERO
  enddo


  call opt_calc_quad(ispin, adot, filegrid, io_grid, nint_rad, ehrad,      &
                      e_rad_re, e_rad_im, nx,ny,nz, nhrad, nrange)


! as the range is small e_rad_re is completely wrong, but doesn't change much on that range

  call grid_interp(ehopt, e_re, nhopt, ehrad(1), ehrad(nhrad), nhrad, e_rad_re, nordp1, dymax)

  if(dymax/e_rad_re(1) > 0.001_REAL64) then
    write(6,*)
    write(6,'("   WARNING:  estimate of error in epsilon_re is ",e12.3)') dymax
    write(6,*)
  endif

  filename_out = 'optical_rad_epsilon_im.dat'
  func = 'Imaginary Part of Dielectric Function             '

  call opt_write_file(trim(filename_out), iotape, title, subtitle, func, &
                      nhrad, ehrad, e_rad_im, vcell, ztot)


! Estimate based on JAP 86,3241 (1999)

  write(6,*) '   Do you want an ESTIMATE based on DOS effective masses?'
  read(5,*) yesno
  write(ioreplay,*) yesno,'     based on DOS'
  write(6,*)

  if(yesno == 'y' .or. yesno == 'Y') then

    write(6,*)
    write(6,*) "   Results based on DOS.  They are very rough!!!"
    write(6,*)


    call opt_rad_rough_pcvsq(filegrid, io_grid, egapopt, pcvsq, nvtc, nx,ny,nz)

    write(6,*)
    write(6,'("  rough estimate of optical gap: ",f12.6," eV")') egapopt*HARTREE
    write(6,'("  rough estimate of Pcv*Pcv:     ",e14.3," au")') pcvsq
    write(6,'("  rough estimate of Pcv:         ",e14.3," SI")') ELMASS*BOHR*sqrt(pcvsq)/AUT
    write(6,*)

    call dos_f_level(tempk,nint(ztot),nhist,ehist,dhist,chist,           &
      lidos,(evbm+ecbm)/2-ezero,vcell, nvbm, ncbm, xmu, xmdosv, xmdosc, xnimass)

    egap = ehist(ncbm) - ehist(nvbm)

    if(egapopt - egap > 2*tau) then

      write(6,*)
      write(6,*) '   WARNING the semiconductor is indirect. Expect large error'
      pref = exp((egapopt - egap) / tau)
      write(6,'("    could be wrong by a factor of ",e12.3)') pref
      write(6,*)

    endif

    do i = 1,nhopt
      xn_re = e_re(i)
      if(ehopt(i) > egapopt) exit
    enddo

    xlth = sqrt(2*PI/tau)
    pref = ALPHA * xlth * sqrt(UM / (xmdosv + xmdosc))
    pref = pref*pref*pref
    bcoef = pref * xn_re * pcvsq * egap * (UM + 3*tau/(2*egap))



    write(6,*)
    write(6,'("  Estimate of recombination coefficient B = ",f12.6," au")') bcoef
    write(6,'("  Estimate of recombination coefficient B = ",e11.3,               &
              & " m^3 / s  (SI)")') bcoef*BOHR*BOHR*BOHR/AUT
    write(6,*)
    write(6,*)
    write(6,'("  Estimate of recombination coefficient B = ",e11.3,               &
              & " cm^3 / s  (cgs)")') 1.0E+6*bcoef*BOHR*BOHR*BOHR/AUT
    write(6,*)
    write(6,*)

  endif

! calculates the intrisic carrier concentration

  call opt_rad_np(xnp, egap, tau, adot, (evbm+ecbm)/2-ezero, ehist, dhist, nhist)

  nhradmin = 1
  do i = 1,nhrad
    if(e_rad_im(i) > EPS) exit
    nhradmin = i
  enddo
  
  do i = nhradmin,nhrad
    e_mod = sqrt(e_rad_re(i)*e_rad_re(i) + e_rad_im(i)*e_rad_im(i))
    xn_im = sqrt( (e_mod - e_rad_re(i)) / 2 )
    xn_re = sqrt( (e_mod + e_rad_re(i)) / 2 )
    abscoef = 2*ehrad(i)*xn_im*ALPHA
    rsrad(i) = ehrad(i)*ehrad(i) * xn_re*xn_re * abscoef * exp(-(ehrad(i)-egap)/tau)
  enddo

! for use in the future

  i = nhradmin
  e_mod = sqrt(e_rad_re(i)*e_rad_re(i) + e_rad_im(i)*e_rad_im(i))
  xn_re = sqrt( (e_mod + e_rad_re(i)) / 2 )

  rsint = ZERO
  do i = nhradmin,nhrad-1
    rsint = rsint + (rsrad(i)+rsrad(i+1)) * (ehrad(i+1)-ehrad(i)) / 2
  enddo

  pref = ALPHA / PI
  pref = pref*pref
  rsint = pref*rsint
  
  write(6,'("  Energy gap ",f12.3," precision",f12.5)') egap*HARTREE, deltae*HARTREE
  xni = sqrt(xnp*exp(-egap/tau))
  xnisi = xni/(BOHR**3)
  write(6,'("  intrinsic carrier concentration ",e12.3," au",5x,e12.3," SI",5x,e12.3,         &
               &  " per cubic cm" )') xni, xnisi, 1.0E-6*xnisi
  write(6,*)

  write(6,'("  radiative recombination rate  ",e12.3," au")') rsint
  write(6,'("  radiative recombination coefficient B  ",e12.3," au")') rsint/xnp
  write(6,'("  Radiative recombination coefficient B = ",e11.3,               &
             " SI (cubic meter per second)")') (rsint/xnp)*BOHR*BOHR*BOHR/AUT
  write(6,*)
  write(6,*)
  write(6,'("  Radiative recombination coefficient B = ",e11.3," cm^3 / s  (cgs) ")')       &
                  1.0E+6*(rsint/xnp)*BOHR*BOHR*BOHR/AUT
  write(6,*)
  write(6,*)


  call execute_command_line("rm " // filegrid // " 2> /dev/null ")

  deallocate(ev_grid)

  deallocate(ehist)
  deallocate(chist)
  deallocate(dhist)

  deallocate(ehopt)
  deallocate(e_im)
  deallocate(e_re)
  deallocate(ehrad)
  deallocate(e_rad_im)
  deallocate(e_rad_re)
  deallocate(rsrad)


end subroutine opt_rad_sub

