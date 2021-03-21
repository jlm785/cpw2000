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

!>  Calculates the density of states using the
!>  information on the PW_DOS_xyz files from cpw_postprocess

subroutine dos_sub(ioreplay)

! rewritten 12 may 2004.
! modified November 16, 2013. jlm
! Modified, documentation, 19 September 2020. JLM
! Modified, egrid, 16-20 October 2020. JLM
! Modified, new unformatted input file. 13 December 2020. JLM
! copyright  J.L.Martins, INESC-MN.

! version 4.99 of cpw

  implicit none

  integer, parameter  :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)              :: ioreplay                           !<  tape number for reproducing calculations

! dimensions of arrays

  integer                          ::  mxdbnd                            !  size of number of bands array
  integer                          ::  nrk                               !  number of k-points
  integer                          ::  nx(3)                             !  number of k-points in each direction in regular grid

! allocatable arrays

  real(REAL64), allocatable        ::  el(:,:)                           !  eigenvalues in Hartree in irreducible wedge
  integer, allocatable             ::  nband(:)                          !  number of bands for each k-points
  real(REAL64), allocatable        ::  wgk(:)                            !  weight for each k-point
  integer, allocatable             ::  indk(:,:)                         !  index of 6 neighboring points
  integer, allocatable             ::  kmap(:,:,:,:)                     !  kmap(1,...) corresponding k-point, kmap(2,...) = 1 additional inversion, kmap(3,...) symmetry operation
  real(REAL64), allocatable        ::  egrid(:,:,:,:)                    !  eigenvalues in Hartree in regular grid
  real(REAL64), allocatable        ::  ehist(:)                          !  histogram energy
  real(REAL64), allocatable        ::  dhist(:)                          !  histogram dos
  real(REAL64), allocatable        ::  chist(:)                          !  histogram integrated dos

  real(REAL64), allocatable        ::  e_of_k(:,:)                       !  band energies of k-point in plot
  real(REAL64), allocatable        ::  e_of_k_so(:,:)                    !  spin-orbit band energies of k-point in plot

  real(REAL64), allocatable        ::  rk(:,:)                           !  component in lattice coordinates of the k-point in the mesh

! scalar variables

  character(len=50)     ::  title                                        !  title for plots
  character(len=140)    ::  subtitle                                     !  subtitle for plots

  logical               ::  lscl                                         !  write for scalar (no-spin-orbit)
  logical               ::  lso                                          !  also write for spin-orbit
  integer               ::  identif                                      !  identifier that is almost unique to the calculation

  real(REAL64)          ::  adot(3,3)                                    !  metric in direct space
  integer               ::  ntrans                                       !  number of symmetry operations in the factor group
  integer               ::  mtrx(3,3,48)                                 !  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the factor group


  integer            ::  nhist                                           !  number of points in histogram
  real(REAL64)       ::  emin, deltae, vol
  logical            ::  lunif                                           !  true if grid is uniform
  real(REAL64)       ::  ztot                                            !  total charge density (electrons/cell)
  real(REAL64)       ::  vcell, bdot(3,3)                                !  unit cell volume (in atomic units)
  integer            ::  nvbm, ncbm                                      !  approximate indices for conduction and valence band
  real(REAL64)       ::  evbb                                            !  energy of bottom of valence band
  real(REAL64)       ::  evbm                                            !  energy of maximum of valence band
  real(REAL64)       ::  ecbm                                            !  energy of minimum of conduction band
  logical            ::  linsul                                          !  true if it appears to be an insulator

  real(REAL64)       ::  tempk                                           !  temperature in K

  character(len = 1) ::  yesno
  integer            ::  ilevel

  real(REAL64)       ::  ezero                                           !  zero of energy
  character(len=40)  ::  filename                                        !  file with band data
  integer            ::  io                                              !  tape number

  integer            ::  ispin                                           !  spin degeneracy (2 for non-spin-polarized, 1 for spin-orbit)
  logical            ::  lidos                                           !  true if integrated density of states is to be computed.
  integer            ::  ichoice
  logical            ::  lask

  real(REAL64)       ::  xmu                                             !  electron chemical potential, AKA Fermi level
  real(REAL64)       ::  xmdosv, xmdosc                                  !  dos effective masses for valence and conduction
  real(REAL64)       ::  xni                                             !  intrinsic carrier concentration

  integer            ::  nhtarg                                          !  target number of points in histogram
  integer            ::  ninter                                          !  number of interpolation points for superquadratic algorithm

! counters

  integer   ::  i, irk

! constants

  real(REAL64), parameter :: HARTREE = 27.21138386_REAL64
  real(REAL64), parameter :: ZERO = 0.0_REAL64 , UM = 1.0_REAL64

  
  io = 40

  write(6,*)
  write(6,*) '   Density of States (DOS) analysis '
  write(6,*)

  write(6,*)
  write(6,*) '   Do you want defaults or full control? (1,2,3)'
  write(6,*)
  write(6,*) '   1) Full range DOS and IDOS '
  write(6,*) '   2) DOS near the gap '
  write(6,*) '   3) Full control of plots and options '
  write(6,*)

  read(5,*) ichoice
  write(ioreplay,*) ichoice,'   dos plot control'

  if(ichoice < 1 .or. ichoice > 3) then

    write(6,*)
    write(6,*) '  Wrong value, enter your choice again (1,2,3) '
    write(6,*)

    read(5,*) ichoice
    write(ioreplay,*) ichoice,'   dos plot control'

    if(ichoice < 1 .or. ichoice > 3) stop

  endif

  filename = 'dos_file.dat'

  if(ichoice == 3) then

    write(6,*) '   What file do you want to use? (1,2)'
    write(6,*)
    write(6,*) ' 1) dos_file.dat (default)'
    write(6,*) ' 2) Choose filename'
    write(6,*)

    read(5,*) i
    write(ioreplay,*) i,'    default file choice'


    if(i /= 1) then

      write(6,*)
      write(6,*) '   Enter desired filename '
      write(6,*)

      read(5,*) filename
      write(ioreplay,*) filename,'   filename'

    endif

  endif

  write(6,*)
  write(6,*) '   Do you want results with spin-orbit (y/n)'
  write(6,*)

  read(5,*) yesno
  write(ioreplay,*) yesno,'   results with spin-orbit'

  if(yesno == 'y' .or. yesno == 'Y') then
    ispin = 1
  else
    ispin = 2
  endif

  call dos_read_size(trim(filename), io, nrk, mxdbnd, nx)

  allocate(e_of_k(mxdbnd,nrk))
  allocate(e_of_k_so(2*mxdbnd,nrk))
  allocate(nband(nrk))
  allocate(wgk(nrk))
  allocate(rk(3,nrk))
  allocate(indk(6,nrk))
  allocate(kmap(3,nx(1),nx(2),nx(3)))

  call dos_read_data(trim(filename), io, title, subtitle,                &
    lscl, lso, identif,                                                  &
    nrk, nx(1), nx(2), nx(3), ztot, adot, ntrans, mtrx,                  &
    nband, rk, wgk, indk, kmap, e_of_k, e_of_k_so,                       &
    nrk, mxdbnd)

  call adot_to_bdot(adot,vcell,bdot)

  if(ispin == 1) then
    if(.not. lso) then
      write(6,*) '    returning in dos_sub:  file does not contain spin-orbit info'

      return

    endif
    mxdbnd = 2*mxdbnd
    allocate(el(mxdbnd,nrk))
    el(:,:) = e_of_k_so(:,:)
    do irk = 1,nrk
      nband(irk) = 2*nband(irk)
    enddo

    write(6,*)
    write(6,*) '     DOS calculated with spin-orbit'
    write(6,*)

  else
    if(.not. lscl) then
      write(6,*) '    returning in dos_sub:  file does not contain scalar info'

      return

    endif
    allocate(el(mxdbnd,nrk))
    el(:,:) = e_of_k(:,:)

    write(6,*)
    write(6,*) '     DOS calculated WITHOUT spin-orbit'
    write(6,*)

  endif

  deallocate(e_of_k)
  deallocate(e_of_k_so)
  deallocate(rk)

! following lines do not depend on change of data file format

  allocate(egrid(nx(1),nx(2),nx(3),mxdbnd))

  call dos_el_to_egrid(el,nrk,nband,nx,kmap,egrid,nx,mxdbnd)

  call dos_bands_rough(el,nrk,nband,ztot,ispin,evbb,evbm,ecbm,linsul,mxdbnd)

  if(ichoice == 1 .or. ichoice == 2) then

    ezero = evbm
    nhtarg = 2000
    ninter = 4

  else

    write(6,*)
    write(6,*)  '  Where do you want the zero of energy?'
    write(6,*)
    write(6,*)  '  0)  Leave it at the average internal potential'
    write(6,*)  '  1)  Use the estimate of valence band maximum'
    write(6,*)  '  2)  Use the estimate of conduction band minimum'
    write(6,*)  '  3)  Use the bottom of valence band'
    write(6,*)  '  4)  Use the estimate of Fermi energy'
    write(6,*)

    read(5,*)   ilevel
    write(ioreplay,*) ilevel,'   zero energy choice'

    ezero = ZERO
    if(ilevel < 0 .or. ilevel > 4) then
      ilevel = 0
      write(6,*) '  wrong value: using default (average potential)'
    elseif(ilevel == 1) then
      ezero = evbm
    elseif(ilevel == 2) then
      ezero = ecbm
    elseif(ilevel == 3) then
      ezero = evbb
    elseif(ilevel == 4) then
      ezero = (evbm + ecbm) / 2
    endif

    write(6,*)
    write(6,*)  '  Enter approximate number of points you want in the plot.'
    write(6,*)  '  Actual number will differ to have nice energy begin, end and step in eV.'
    write(6,*)

    read(5,*)   nhtarg
    write(ioreplay,*) nhtarg,'   target number of points'

    if(nhtarg < 10 .or. nhtarg > 100000) then

      write(6,*)
      write(6,*)  '  The number of points is unreasonable, using 2000 for target.'
      write(6,*)

      nhtarg = 2000

    endif


    write(6,*)
    write(6,*)  '  Enter number of interpolation points for superquadratic algorithm.'
    write(6,*)

    read(5,*)   ninter
    write(ioreplay,*) ninter,'   number of interpolation points'

    if(ninter < 2 .or. ninter > 10) then

      write(6,*)
      write(6,*)  '  The number of points is unreasonable, using 4.'
      write(6,*)

      ninter = 4

    endif

  endif

  write(6,*)
  write(6,'("  The new reference energy (new zero) for the ",            &
      "density of states is at",f10.3," eV")') ezero*HARTREE
  write(6,'("  with respect to the average potential")')
  write(6,*)

  if(ichoice == 1) then
    lidos = .TRUE.
  elseif(ichoice == 2) then
    lidos = .FALSE.
  else

    write(6,*)
    write(6,*) '   Do you want an integrated DOS? (y/n) '
    write(6,*)

    read(5,*) yesno
    write(ioreplay,*) yesno,'   IDOS'

    if(yesno == 'y' .or. yesno == 'Y') then
      lidos = .TRUE.
    else
      lidos = .FALSE.
    endif

  endif

  if(ichoice == 3) then

    write(6,*)
    write(6,*)  '  Do you want a linear interpolation plot? (y/n)'
    write(6,*)  '  This has lower resolution...'
    write(6,*)

    read(5,*) yesno
    write(ioreplay,*) yesno,'   linear interpolation'

    if(yesno == 'y' .or. yesno == 'Y') then
 
      call dos_plot_range(nx,egrid,nhtarg,emin,deltae,nhist,ezero,nx,mxdbnd)

      allocate(ehist(nhist))
      allocate(chist(nhist))
      allocate(dhist(nhist))

      do i=1,nhist
        ehist(i) = emin + deltae*(i-1)
        dhist(i) = ZERO
        chist(i) = ZERO
      enddo
      lunif = .TRUE.

      vol = UM / (nx(1)*nx(2)*nx(3))

      call dos_lin(nx,egrid, nhist,ehist,dhist,chist,ezero,ispin,        &
               lunif,.TRUE.,lidos,vol,mxdbnd)

      write(6,*) ' Do you want to see a rough DOS now?(y/n)'
      write(6,*)

      read(5,*) yesno
      write(ioreplay,*) yesno,'   rough plot'

      if(yesno == 'y' .or. yesno == 'Y') then

        call dos_print_ascii(emin,deltae,nhist,dhist,chist,lidos,20,60)

      endif

      filename = 'DOS_LIN.DAT'
      call dos_print_file(trim(filename),nhist,ehist,dhist,chist,        &
             lidos,ezero,(evbm+ecbm)/2-ezero,vcell,ztot)

      write(6,*)
      write(6,*) '  The file ',filename,' has been written.'
      write(6,*)

      deallocate(ehist)
      deallocate(chist)
      deallocate(dhist)

    endif

    write(6,*) ' Do you want to see a gaussian DOS now?(y/n)'
    write(6,*)

    read(5,*) yesno
    write(ioreplay,*) yesno,'   gaussian broadening'

    if(yesno == 'y' .or. yesno == 'Y') then

      call dos_gau(el,nrk,nband,wgk,ezero,lidos,100,mxdbnd)

    endif

  endif

! begin of the main calculations

  write(6,*)
  write(6,*) '  Calculating the DOS with a "quadratic" interpolation.'
  write(6,*)

  if(ichoice == 2) then

    if(linsul) then

      deltae = 0.001_REAL64/HARTREE
      emin = evbm - ezero - 0.5/HARTREE
      emin = deltae*real(int(emin/deltae))
      nhist = int((ecbm-evbm + 1.0/HARTREE) / deltae)

    else

      deltae = 0.001_REAL64/HARTREE
      emin = ezero - 1.0/HARTREE
      emin = deltae*real(int(emin/deltae)-5)
      nhist = int(2.0 / deltae) + 10

    endif

  else

    call dos_plot_range(nx,egrid,nhtarg,emin,deltae,nhist,ezero,nx,mxdbnd)

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

  call dos_quad(nx,egrid,4, nhist,ehist,dhist,chist,ezero,ispin,         &
        lunif,.TRUE.,lidos,mxdbnd)

  if(ichoice == 3) then

    write(6,*) ' Do you want to see a rough DOS now?(y/n)'
    write(6,*)

    read(5,*) yesno
    write(ioreplay,*) yesno,'   rough plot'

    if(yesno == 'y' .or. yesno == 'Y') then

      call dos_print_ascii(emin,deltae,nhist,dhist,chist,lidos,100,60)

    endif

  endif

  filename = 'DOS_QUAD.DAT'

  call dos_print_file(trim(filename),nhist,ehist,dhist,chist,       &
         lidos,ezero,(evbm+ecbm)/2-ezero,vcell,ztot)

  write(6,*)
  write(6,*) '  The file ',filename,' has been written.'
  write(6,*)

  filename = 'dos_quad.gp'

  if(ichoice == 3) then
    lask = .TRUE.
  else
    lask = .FALSE.
  endif

  call dos_out_gnuplot(ioreplay,trim(filename),                     &
         nhist,ehist,dhist,chist,lidos,lask)

! Fermi level and effective masses

  if(ichoice == 3) then

    do i=1,100
      write(6,*)
      write(6,'("  Enter desired temperature (in K) for carrier ",  &
        " concentration and Fermi level. (Negative exits)")')
      read(5,*) tempk
      write(ioreplay,*) tempk,'   temperature (K)'

      if(tempk < 0.0) exit

      call dos_f_level(tempk,nint(ztot),nhist,ehist,dhist,chist,    &
      lidos,(evbm+ecbm)/2-ezero,vcell, nvbm, ncbm, xmu, xmdosv, xmdosc, xni)

    enddo

  else

    tempk = 300.0
    call dos_f_level(tempk,nint(ztot),nhist,ehist,dhist,chist,      &
      lidos,(evbm+ecbm)/2-ezero,vcell, nvbm, ncbm, xmu, xmdosv, xmdosc, xni)

  endif

  deallocate(ehist)
  deallocate(chist)
  deallocate(dhist)

  deallocate(el)
  deallocate(nband)
  deallocate(wgk)
  deallocate(indk)
  deallocate(kmap)
  deallocate(egrid)

  return
end subroutine dos_sub






