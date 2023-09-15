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

!>  Calculates the layer average and double average of
!>  the charge density and electrostatic potential.
!>
!>  For the double average see  PRL 61, 734 (1988).
!>
!>  \author       Jose Luis Martins
!>  \version      5.06
!>  \date         September 5, 2012. 2 April 2023.
!>  \copyright    GNU Public License v2


subroutine plot_rho_v_average(ioreplay,                                  &
         adot, ntype, natom, nameat, rat, zv,                            &
         ng, kgv, phase, conj, ns, mstar,                                &
         den,                                                            &
         mxdtyp, mxdatm, mxdgve, mxdnst)


! writen September 5, 2012.jlm
! Modified, split code, complex variables, 26 May 2014. JLM
! Documentation, merge psi_plot, 5 February 2021. JLM
! Double average by material.  March-April 2023. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                :: ioreplay                         !<  tape number for reproducing calculations

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxdnst                          !<  array dimension for g-space stars

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space
  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i
  character(len=2), intent(in)       ::  nameat(mxdtyp)                  !<  chemical symbol for the type i
  real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)            !<  k-th component (in lattice coordinates) of the position of the n-th atom of type i
  real(REAL64), intent(in)           ::  zv(mxdtyp)                      !<  valence of atom with type i

  integer, intent(in)                ::  ng                              !<  total number of g-vectors with length less than gmax
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length
  complex(REAL64), intent(in)        ::  phase(mxdgve)                   !<   phase factor of G-vector n
  real(REAL64), intent(in)           ::  conj(mxdgve)                    !<  is -1 if one must take the complex conjugate of x*phase
  integer, intent(in)                ::  ns                              !<  number os stars with length less than gmax
  integer, intent(in)                ::  mstar(mxdnst)                   !<  number of g-vectors in the j-th star

  complex(REAL64), intent(in)        ::  den(mxdnst)                     !<  valence charge density for the prototype g-vector in star j

! allocatable arrays

  complex(REAL64), allocatable       ::  rho(:)                          !  charge density for G-vector j
  complex(REAL64), allocatable       ::  rhogau(:)                       !  atom centered gaussian charge density for G-vector j
  complex(REAL64), allocatable       ::  gtmp(:)                         !  temporary array for G-vector quantities
  real(REAL64), allocatable          ::  width(:)                        !  width for the double average
  real(REAL64), allocatable          ::  widthgeom(:)                    !  width for the double average from material info
  real(REAL64), allocatable          ::  widthrho(:)                     !  width for the double average from charge density
  integer, allocatable               ::  izval(:)                        !  valence of atom of type i
  real(REAL64), allocatable          ::  xpeak(:)                        !  position of first peak of auto-correlation for atom of type i

  real(REAL64), allocatable   ::  ave(:)                                 !  layer average of charge density in the 3d direction of fft grid
  real(REAL64), allocatable   ::  dave(:,:)                              !  double average of charge density in the 3d direction of fft grid
  real(REAL64), allocatable   ::  gave(:)                                !  average of the nuclear "gaussian" charge density in the 3d direction of fft grid
  real(REAL64), allocatable   ::  conv(:)                                !  repeated average of charge density by convolution in the 3d direction of fft grid

  real(REAL64), allocatable   ::  autocorr(:)                            !  autocorrelation functiom

  integer, allocatable        ::  iptype(:)                              !  type of atom with n-th lower z coordinate
  integer, allocatable        ::  ipnatom(:)                             !  atom of iptupe with n-th lower z coordinate

  integer, allocatable        ::  nrepeat(:)                             !  number of repeat units for each material
  real(REAL64), allocatable   ::  rbottom(:)                             !  bottom lattice coordinate of the material
  real(REAL64), allocatable   ::  yave(:)                                !  average of the function in the center of the material.

! main variables

  integer                     ::  nmat                                   !  number of materials
  integer                     ::  nptot                                  !  total number of repeat units

  real(REAL64)                ::  xave                                   !  distance over which the average is made

  logical                     ::  linter                                 !  if .TRUE. asks if the figure should be shown
  character(len=1)            ::  yesno

  integer                     ::  ntot                                   !  total number of atoms

! parameters that control the quality of the plots

  integer                     ::  nmult                                  !  increases point density in z direction
  real(REAL64)                ::  sigref                                 !  controls the width of the core gaussian

! other variables

  integer             ::  kmscr(3), nsfft(3), ktmp(3)                    !  fft sizes
  integer             ::  mfft, mwrk                                     !  fft stuff

  integer             ::  istar, istop
  real(REAL64)        ::  xfirst                                         !  autocorrelation peak

  real(REAL64)        ::  gsquare, vcell, bdot(3,3), fac
  integer             ::  id, n1,n2,n3, nn
  integer             ::  nwidth
  integer             ::  mxdwid

  real(REAL64)        ::  sigma                                          !  width of gaussian (in real space)

  integer             ::  iotape
  character(len=40)   ::  filename

  real(REAL64)        ::  height                                         !  height of the cell

  logical             ::  lfound, lopt, lfoundall

! constants

  real(REAL64), parameter  ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  real(REAL64), parameter  ::  PI=3.141592653589793_REAL64
  real(REAL64), parameter  ::  HARTREE = 27.21138386_REAL64
  real(REAL64), parameter  ::  BOHR = 0.5291772109_REAL64

! counters

  integer      :: i, j, k, nt

! These constants may be changed to get nicer plots.

  nmult = 3
  sigref = 0.5

! finds basic fft grid

  write(6,*)
  write(6,'("  Do you want to see the plots interactively? (y/n)")')
  write(6,*)

  read(5,*) yesno
  write(ioreplay,'(2x,a1,"   interactive plot")') yesno

  if(yesno == 'y' .or. yesno == 'Y') then

    linter = .TRUE.

    write(6,*)
    write(6,'("  Program will generate files for later plotting ",       &
            & "with gnuplot")')
    write(6,'("  BEWARE: plots may hide below each other")')
    write(6,*)

  else

    linter = .FALSE.

    write(6,*)
    write(6,'("  Program will generate files for later plotting ",       &
            & "with xmgrace")')
    write(6,*)

  endif


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
  write(6,'("  the basic fft grid is: ",3i8)') (nsfft(j),j=1,3)
  write(6,*) '  you can change it by changing the energy cutoff'
  write(6,*)

  allocate(iptype(mxdatm*mxdtyp),ipnatom(mxdatm*mxdtyp))

  call plot_z1D_print_ordered(1, 6, ntype, natom, nameat, rat,           &
         ntot, iptype, ipnatom,                                          &
         mxdtyp, mxdatm)

  write(6,*)
  write(6,'("  Enter number of different materials")')
  write(6,'("  Must be between 1 and ",i5)') ntot
  write(6,*)

  read(5,*) nmat
  write(ioreplay,'(2x,i8,"   number of materials")') nmat

  if(nmat <1 .or. nmat > ntot) then
    write(6,*)
    write(6,'("  Wrong value enter it again")')
    write(6,*)

    read(5,*) nmat
    write(ioreplay,'(2x,i8,"   number of materials again")') nmat

    if(nmat <1 .or. nmat > ntot) then
      write(6,*)
      write(6,'("  STOPPING  unreasonable value")')
      write(6,*)

      stop

    endif
  endif

  if(nmat == 1) then
    write(6,*)
    write(6,'("  Only one type of material...")')
    write(6,*)
  endif

  allocate(nrepeat(nmat))
  allocate(rbottom(nmat))
  allocate(yave(nmat))

  allocate(widthgeom(nmat))

  nwidth = nmat

  call plot_z1D_material_width(ioreplay,                                 &
       ntot, iptype, ipnatom, rat, adot,                                 &
       height, nptot, nrepeat, rbottom, widthgeom,                       &
       nmat, mxdtyp, mxdatm)

! make it more dense on the averaging direction

  ktmp(1) = kmscr(1)
  ktmp(2) = kmscr(2)
  ktmp(3) = nmult*kmscr(3)

! tries to get a multiple of nmat (if nmat is a prime number it just increases the number of points)

  do i = 1,5
    call size_fft(ktmp,nsfft,mfft,mwrk)

    if(mod(nsfft(3),nmat) == 0) exit

    ktmp(3) = nsfft(3) / 2
  enddo

!  kave = nsfft(3) / nmat

  write(6,*)
  write(6,'("  the new fft grid is: ",3i8)') (nsfft(j),j=1,3)
  write(6,*)

! unfold the charge density and store in phase array

  allocate(rho(ng))
  allocate(rhogau(ng))
  allocate(gtmp(ng))
  allocate(izval(mxdtyp))
  allocate(xpeak(mxdtyp))

! if inds was available should call charge_unfold...

  rho = ZERO

  istop = 0
  do i=1,ns
    istar = istop+1
    istop = istar+mstar(i)-1
    do j=istar,istop
      rho(j) = den(i)*conjg(phase(j))
      if(conj(j) < ZERO) rho(j) = conjg(rho(j))
    enddo
  enddo

! choice of width for gaussian broadening

! sigma = 0.02*height/nmat
  sigma = sigref

  n1 = nsfft(1)
  n2 = nsfft(2)
  n3 = nsfft(3)

  nn = nsfft(3)

  id = n1 + 1

! initialize charge density array and enter symmetrized
! charge.

  allocate(ave(nn))
  allocate(dave(nn,1))
  allocate(gave(nn))

  allocate(autocorr(0:nn/2))

  iotape = 11
  if(iotape == ioreplay) iotape = iotape + 1


  do nt = 1,ntype

    do j=1,ntype
      izval(j) = 0
    enddo

    izval(nt) = nint(zv(nt))

    call plot_gauss(sigma, rhogau,                                       &
         adot, ntype, natom, rat, izval,                                 &
         ng, kgv,                                                        &
         mxdtyp, mxdatm, mxdgve)


    call plot_zave1D(gave, rhogau, nptot, id, n1,n2,n3, ng, kgv)

    if(linter) then
      filename = 'rho_gauss_' // adjustl(trim(nameat(nt)))//'.gp'
      call plot_z1D_gnuplot(ioreplay, iotape, gave, dave, 0, n3, height, &
            adjustl(trim(filename)),                                     &
            'Broadened Nuclear Density '//nameat(nt),                    &
            '{/Symbol r} (1/cell)',linter)
    else
      filename = 'rho_gauss_'// adjustl(trim(nameat(nt)))//'.agr'
      call plot_z1D_xmgr(iotape, gave, dave, 0, n3, height,              &
            adjustl(trim(filename)),                                     &
           'Broadened Nuclear Density '//nameat(nt),                     &
           '\f{Symbol} r\f{} (1/cell)')
    endif

!   auto-correlation function.

    call plot_z1D_auto_corr(nn, gave, nn/nptot, autocorr, xpeak(nt), lfound)

    if(.not. lfound) then
      write(6,*)
      write(6,*) '  peak for atom ',nameat(nt),' may not be accurate'
      write(6,*)
    endif

    if(linter) then
      filename = 'rho_nucl_' // adjustl(trim(nameat(nt))) //             &
                     '_auto_corr.gp'
      call plot_z1D_gnuplot(ioreplay, iotape, autocorr,dave,             &
            0, n3/2, height/2,  adjustl(trim(filename)),                 &
            nameat(nt)//' Nuclear Density Auto-correlation',             &
            '  ',linter)
    else
      filename = 'rho_nucl_' // adjustl(trim(nameat(nt))) //             &
                     '_auto_corr.agr'
      call plot_z1D_xmgr(iotape, autocorr, dave, 0, n3/2, height/2,      &
            adjustl(trim(filename)),                                     &
            nameat(nt)//' Nuclear Density Auto-correlation', '  ')
    endif

  enddo


  call plot_zave1D(ave, rho, nptot, id, n1,n2,n3, ng, kgv)

  if(linter) then
    call plot_z1D_gnuplot(ioreplay, iotape, ave, dave, 0, n3,            &
          height, 'rho_ave.gp','Electron Density Average',               &
          '{/Symbol r} (e/cell)', linter)
  else
    call plot_z1D_xmgr(iotape, ave, dave ,0, n3, height,                 &
          'rho_ave.agr',                                                 &
          'Electron Density Average', '\f{Symbol} r\f{} (e/cell)')
  endif


! quick and dirty auto-correlation function. See Wiener–Khinchin theorem

  call plot_z1D_auto_corr(nn, ave, nn/nptot, autocorr, xfirst, lfound)

  if(linter) then
    call plot_z1D_gnuplot(ioreplay, iotape, autocorr, dave,              &
          0, n3/2, height/2, 'rho_auto_corr.gp',                         &
          'Electron Density Auto-correlation', '  ',linter)
  else
    call plot_z1D_xmgr(iotape, autocorr, dave, 0, n3/2, height/2,        &
         'rho_auto_corr.agr', 'Electron Density Auto-correlation',       &
            '  ')
  endif

! find widths from auto-correlation of localy projected electron density

  allocate(widthrho(nmat))

  call plot_z1D_local_corr(nn, ave, nmat, height, rbottom, nrepeat, widthrho)


! chooses convolution widths from several options

  mxdwid = max(nmat,4)
  allocate(width(mxdwid))

  call plot_z1d_choice_width(ioreplay, nwidth, width, lopt,              &
         ntype, natom, nameat,                                           &
         height, nmat, nrepeat, widthgeom, widthrho,                     &
         nn, xfirst, xpeak,                                              &
         mxdtyp, mxdwid)

! uses convolution for the double average for electron density,
! electron plus gaussian broadened ion density and electrostatic potential

  deallocate(dave)
  allocate(dave(nn,nwidth))

  allocate(conv(nn))

! Electron density

  do k = 1,nwidth

    xave = width(k)/height
    call plot_convol_opt(n3, xave, ave, dave(:,k), yave(k), lfound, lopt,  &
         nrepeat(k), rbottom(k), rbottom(mod(k,nmat)+1))

  enddo

  if(linter) then
    call plot_z1D_gnuplot(ioreplay, iotape, ave, dave, nwidth,           &
          n3,height, 'rho_dave.gp',                                      &
          'Electon Density','{/Symbol r} (e/cell)', linter)
  else
    call plot_z1D_xmgr(iotape, ave, dave, nwidth, n3, height,            &
          'rho_dave.agr',                                                &
          'Electon Density','\f{Symbol} r\f{} (e/cell)')
  endif

! repeats for the electron+ion charge density
! use smaller sigma

  sigma = 2.0*height/nsfft(3)

  do nt = 1,ntype

    izval(nt) = nint(zv(nt))

  enddo

    call plot_gauss(sigma, rhogau,                                       &
         adot, ntype, natom, rat, izval,                                 &
         ng, kgv,                                                        &
         mxdtyp, mxdatm, mxdgve)


  do i=1,ng
    gtmp(i) = rho(i) - rhogau(i)
  enddo

  call plot_zave1D(ave, gtmp, nptot, id, n1,n2,n3, ng, kgv)


  if(linter) then
    call plot_z1D_gnuplot(ioreplay, iotape, ave, dave, 0, n3, height,    &
         'rho_total_ave.gp', 'Total Charge Density Average',             &
         '{/Symbol r} (e/cell)', linter)
  else
    call plot_z1D_xmgr(iotape, ave, dave, 0, n3, height,                 &
          'rho_total_ave.agr', 'Total Charge Density Average',           &
            '\f{Symbol} r \f{} (e/cell)')
  endif

! each width convoluted separately
! uses convolution for the double average

  do k=1,nwidth

    xave = width(k)/height
    call plot_convol_opt(n3, xave, ave, dave(:,k), yave(k), lfound, lopt,  &
         nrepeat(k), rbottom(k), rbottom(mod(k,nmat)+1))

  enddo


  if(linter) then
    call plot_z1D_gnuplot(ioreplay, iotape, ave, dave, nwidth, n3,       &
          height, 'rho_total_dave.gp', 'Total Charge Density',           &
          '{/Symbol r} (e/cell)', linter)
  else
    call plot_z1D_xmgr(iotape, ave, dave, nwidth, n3, height,            &
          'rho_total_dave.agr','Total Charge Density',                   &
          '\f{Symbol} r\f{} (e/cell)')
  endif


! repeats for the Hartree potential

  call adot_to_bdot(adot, vcell, bdot)

  fac = 4*PI / vcell

  do i=2,ng
    gsquare = ZERO
    do j=1,3
    do k=1,3
      gsquare = gsquare + kgv(j,i)*bdot(j,k)*kgv(k,i)
    enddo
    enddo
    gtmp(i) = fac*(rho(i) - rhogau(i))/gsquare
  enddo

  gtmp(1) = cmplx(ZERO,ZERO,REAL64)

  call plot_zave1D(ave, gtmp, 1, id, n1,n2,n3, ng, kgv)

  do i=1,n3
    ave(i) = HARTREE*ave(i)
  enddo


  if(linter) then
    call plot_z1D_gnuplot(ioreplay, iotape, ave, dave, 0, n3,            &
          height, 'pot_ave.gp', 'Electrostatic Potential Average',       &
         'V (eV)', linter)
  else
    call plot_z1D_xmgr(iotape, ave, dave, 0, n3, height,                 &
          'pot_ave.agr', 'Electrostatic Potential Average', 'V (eV)')
  endif

  do i=1,n3
    ave(i) = ave(i) / HARTREE
  enddo

! each width convoluted separately
! uses convolution for the double average

  lfoundall = .TRUE.
  do k=1,nwidth

    xave = width(k)/height
    call plot_convol_opt(n3, xave, ave, dave(:,k), yave(k), lfound, lopt,  &
         nrepeat(k), rbottom(k), rbottom(mod(k,nmat)+1))

    if(.not. lfound) lfoundall = .FALSE.

  enddo

  do i=1,n3
    ave(i) = ave(i) * HARTREE
    do k = 1,nwidth
      dave(i,k) = dave(i,k) * HARTREE
    enddo
  enddo

  if(linter) then
    call plot_z1D_gnuplot(ioreplay, iotape, ave, dave, nwidth, n3,       &
          height, 'pot_dave.gp', 'Electrostatic Potential',              &
          'V (eV)', linter)
  else
    call plot_z1D_xmgr(iotape, ave, dave, nwidth, n3, height,            &
          'pot_dave.agr', 'Electrostatic Potential', 'V (eV)')
  endif

! prints the information for use in band alignments

  if(nmat > 1) then

    write(6,*)
    write(6,*) '   Electrostatic potential steps in eV'
    write(6,*) '   Check the figures to see it they are accurate!'
    write(6,*)
    if(.not. lfoundall .and. lopt) write(6,*) '   WARNING    optimization did not work!!!'
    write(6,*)

    do j = 1,nmat-1
    do k = j+1,nmat
      write(6,'("    Step between material ",i3," and ",i3," is ",       &
          &  f12.3)')   k,j,(yave(k)-yave(j))*HARTREE
    enddo
    enddo

  endif


! deallocates the stuff

  deallocate(iptype,ipnatom)
  deallocate(nrepeat)
  deallocate(rbottom)
  deallocate(yave)

  deallocate(ave)
  deallocate(dave)
  deallocate(gave)

  deallocate(rho)
  deallocate(rhogau)

  deallocate(gtmp)
  deallocate(izval)
  deallocate(width)
  deallocate(widthgeom)
  deallocate(xpeak)

  deallocate(autocorr)
  deallocate(conv)

  return

end subroutine plot_rho_v_average
