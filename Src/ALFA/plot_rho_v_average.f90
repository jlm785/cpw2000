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
  real(REAL64), allocatable          ::  widthold(:)                     !  width for the double average in earlier code
  integer, allocatable               ::  izval(:)                        !  valence of atom of type i
  real(REAL64), allocatable          ::  xpeak(:)                        !  position of first peak of auto-correlation for atom of type i

  real(REAL64), allocatable   ::  ave(:)                                 !  layer average of charge density in the 3d direction of fft grid
  real(REAL64), allocatable   ::  dave(:,:)                              !  double average of charge density in the 3d direction of fft grid
  real(REAL64), allocatable   ::  gave(:)                                !  average of the nuclear "gaussian" charge density in the 3d direction of fft grid
  real(REAL64), allocatable   ::  conv(:)                                !  repeated average of charge density by convolution in the 3d direction of fft grid

  integer, allocatable        ::  indx(:)                                !  index of species with most atoms

  real(REAL64), allocatable   ::  autocorr(:)                            ! autocorrelation functiom
  real(REAL64), allocatable   ::  convtmp(:)                             !  repeated average of charge

  integer, allocatable        ::  iptype(:)
  integer, allocatable        ::  ipnatom(:)
  integer, allocatable        ::  nrepeat(:)
  real(REAL64), allocatable   ::  rleft(:)

! main variables

  integer                     ::  nmat                                   !  number of materials
  integer                     ::  nptot                                  !  total number of repeat units

  real(REAL64)                ::  xave                                   !  distance over which the average is made

  logical                     ::  linter                                 !  if .TRUE. asks if the figure should be shown
  character(len=1)            ::  yesno
  integer                     ::  nab, ntotal, nwd

  integer                     ::  ntot

! parameters that control the quality of the plots

  integer                     ::  nmult                                  !  increases point density in z direction
  real(REAL64)                ::  sigref                                 !  controls the width of the core gaussian

! other variables

  integer             ::  kmscr(3), nsfft(3)
  integer             ::  istar, istop
  real(REAL64)        ::  xfirst

  real(REAL64)        ::  gsquare, vcell, bdot(3,3), fac
  integer             ::  id, n1,n2,n3, nn
  integer             ::  ku, kd
  integer             ::  nwidth
  integer             ::  nwidthold
  real(REAL64)        ::  sigma                      !  width of gaussian (in real space)

  integer             ::  iotape
  character(len=40)   ::  filename
  integer             ::  mfft, mwrk
  integer             ::  ktmp(3)

  real(REAL64)        ::  rright, height

  logical             ::  lfound
  integer             ::  ichoice

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
  write(6,'("  Enter number of different materials  ")')
  write(6,*)

  read(5,*) nmat
  write(ioreplay,'(2x,i8,"   number of materials")') nmat

  allocate(nrepeat(nmat))
  allocate(rleft(nmat))

  allocate(widthgeom(nmat))

  nwidth = nmat

  call plot_z1D_material_width(ioreplay,                                 &
       ntot, iptype, ipnatom, rat, adot,                                 &
       height, nptot, nrepeat, rleft, widthgeom,                         &
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
!  ntot = id * n2 * n3

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
      filename = 'rho_nucl_' // adjustl(trim(nameat(nt))) //        &
                     '_auto_corr.gp'
      call plot_z1D_gnuplot(ioreplay, iotape, autocorr,dave,        &
            0, n3/2, height/2,  adjustl(trim(filename)),            &
            nameat(nt)//' Nuclear Density Auto-correlation',        &
            '  ',linter)
    else
      filename = 'rho_nucl_' // adjustl(trim(nameat(nt))) //        &
                     '_auto_corr.agr'
      call plot_z1D_xmgr(iotape, autocorr, dave, 0, n3/2, height/2, &
            adjustl(trim(filename)),                                &
            nameat(nt)//' Nuclear Density Auto-correlation', '  ')
    endif

  enddo



  call plot_zave1D(ave, rho, nptot, id, n1,n2,n3, ng, kgv)

  if(linter) then
    call plot_z1D_gnuplot(ioreplay, iotape, ave, dave, 0, n3,       &
          height, 'rho_ave.gp','Electron Density Average',          &
          '{/Symbol r} (e/cell)', linter)
  else
    call plot_z1D_xmgr(iotape, ave, dave ,0, n3, height,            &
          'rho_ave.agr',                                            &
          'Electron Density Average', '\f{Symbol} r\f{} (e/cell)')
  endif


! quick and dirty auto-correlation function. See Wiener–Khinchin theorem

  call plot_z1D_auto_corr(nn, ave, nn/nptot, autocorr, xfirst, lfound)

! find widths from auto-correlation of localy projected electron density

  allocate(widthrho(nmat))

  call plot_z1D_local_corr(nn, ave, nmat, height, rleft, nrepeat, widthrho)

  if(linter) then
    call plot_z1D_gnuplot(ioreplay, iotape, autocorr, dave,         &
          0, n3/2, height/2, 'rho_auto_corr.gp',                    &
          'Electron Density Auto-correlation', '  ',linter)
  else
    call plot_z1D_xmgr(iotape, autocorr, dave, 0, n3/2, height/2,   &
         'rho_auto_corr.agr', 'Electron Density Auto-correlation',  &
            '  ')
  endif

! identifies most abundant species and finds default widths

  allocate(indx(mxdtyp))
  allocate(widthold(mxdtyp))

  call isort(ntype, natom, indx)

  nab = 0
  ntotal = 0
  do nt = 1,ntype
    ntotal = ntotal + natom(nt)
  enddo
  do nt = 1,ntype
    if(4*natom(nt) > ntotal) nab = nab + 1
  enddo

  if(nab == 0) then
    nwidthold = 1
    widthold(1) = height/nptot
  else
    nwidthold = nab
    do j=1,nab
      widthold(j) = xpeak(indx(ntype-j+1))*height/nn
    enddo
  endif

  write(6,*)
  write(6,*)
  write(6,*)
  write(6,*) ' You have to choose both the number and'
  write(6,*) ' width of the square well averages'
  write(6,*)
  write(6,'("  Width from total number of planes: ",g14.6)') height/nptot

  write(6,*)
  write(6,'("  Width from electron density autocorrelation: ",g14.6)')   &
         xfirst*height/nn
  write(6,*)

  do nt = 1,ntype
    write(6,'("  Width from ",a2," atomic autocorrelation: ",g14.6)')    &
         nameat(indx(nt)),xpeak(indx(nt))*height/nn
    write(6,*)
  enddo

  write(6,*)
  write(6,*)
  write(6,*)
  write(6,'("  Old algorithm suggest ",i2," widths with values: ",       &
        &     99f12.5)') nwidthold,(widthold(j),j=1,nwidthold)
  write(6,*)

  write(6,*)
  write(6,'("  Materials description (layer size) suggest ",i2,          &
        &     " widths with values: ",99f12.5)')                         &
               nwidth,(widthgeom(j),j=1,nwidth)
  write(6,*)

  write(6,*)
  write(6,'("  Density auto-correlation for each material suggest ",i2,  &
        &     " widths with values: ",99f12.5)')                         &
               nwidth,(widthrho(j),j=1,nwidth)
  write(6,*)


  write(6,*)
  write(6,*) "  Choose which values you want"
  write(6,*) "  1) Materials description"
  write(6,*) "  2) Density autocorrelation"
  write(6,*) "  3) Old algorithm"
  write(6,*) "  4) Enter other values of your choice"
  write(6,*)
  write(6,*) "  Enter your choice (1--4)"


  read(5,*) ichoice
  write(ioreplay,'(2x,i5,"   choice of widths")') ichoice

  if(ichoice < 1 .or. ichoice > 4) then
    write(6,*)
    write(6,*) "  Wrong value, enter again"
    read(5,*) ichoice
    write(ioreplay,'(2x,i5,"   choice of widths")') ichoice
    if(ichoice < 1 .or. ichoice > 4) then
      write(6,*)
      write(6,*) "  Wrong value, using materials description"
      write(6,*)
      ichoice = 1
    endif
  endif

  allocate(width(nmat))

  if(ichoice == 1) then
    do i = 1,nwidth
      width(i) = widthgeom(i)
    enddo
  elseif(ichoice == 2) then
    do i = 1,nwidth
      width(i) = widthrho(i)
    enddo
  elseif(ichoice == 3) then
    do i = 1,nwidth
      width(i) = widthold(i)
    enddo
  else


    write(6,'(" Enter number of averages ")')
    write(6,*)

    read(5,*)  nwd
    write(ioreplay,'(2x,i8,"   number of average widths")') nwd

    if(nwidth < 0 .or. nwidth > mxdtyp) then

      write(6,*) "  Invalid answer using defaults "

    else

      nwidth = nwd
      write(6,'(" Enter the ",i5," widths ")') nwidth
      read(5,*) (width(i),i=1,nwidth)
      write(ioreplay,'(2x,20f15.5)') (width(i),i=1,nwidth)

    endif

  endif

! After this point widths are set.

  do k=1,nwidth

    if(width(k) <= 0 .or. width(k) > height) then

      write(6,*) "  Invalid widths  "

      stop

    endif

  enddo

! uses convolution for the double average

  deallocate(dave)
  allocate(dave(nn,nwidth))

  allocate(conv(nn))
  allocate(convtmp(nn))

  do k = 1,nwidth

    xave = width(k)/height
    call plot_convol(n3, xave, ave, dave(:,k))

  enddo

  if(linter) then
    call plot_z1D_gnuplot(ioreplay, iotape, ave, dave, nwidth,      &
          n3,height, 'rho_dave.gp',                                 &
          'Electon Density','{/Symbol r} (e/cell)', linter)
  else
    call plot_z1D_xmgr(iotape, ave, dave, nwidth, n3, height,       &
          'rho_dave.agr',                                           &
          'Electon Density','\f{Symbol} r\f{} (e/cell)')
  endif

! repeats for the electron+ion charge density
! use smaller sigma

  sigma = 2.0*height/nsfft(3)

  do nt = 1,ntype

    izval(nt) = nint(zv(nt))

  enddo

    call plot_gauss(sigma, rhogau,                                  &
         adot, ntype, natom, rat, izval,                            &
         ng, kgv,                                                   &
         mxdtyp, mxdatm, mxdgve)


  do i=1,ng
    gtmp(i) = rho(i) - rhogau(i)
  enddo

  call plot_zave1D(ave, gtmp, nptot, id, n1,n2,n3, ng, kgv)


  if(linter) then
    call plot_z1D_gnuplot(ioreplay, iotape, ave, dave, 0, n3, height,  &
         'rho_total_ave.gp', 'Total Charge Density Average',        &
         '{/Symbol r} (e/cell)', linter)
  else
    call plot_z1D_xmgr(iotape, ave, dave, 0, n3, height,            &
          'rho_total_ave.agr', 'Total Charge Density Average',      &
            '\f{Symbol} r \f{} (e/cell)')
  endif


! uses convolution for the double average

  do k=1,nwidth

    xave = width(k)/height
    call plot_convol(n3, xave, ave, dave(:,k))

  enddo


  if(linter) then
    call plot_z1D_gnuplot(ioreplay, iotape, ave, dave, nwidth, n3,  &
          height, 'rho_total_dave.gp', 'Total Charge Density',      &
          '{/Symbol r} (e/cell)', linter)
  else
    call plot_z1D_xmgr(iotape, ave, dave, nwidth, n3, height,       &
          'rho_total_dave.agr','Total Charge Density',              &
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
    call plot_z1D_gnuplot(ioreplay, iotape, ave, dave, 0, n3,       &
          height, 'pot_ave.gp', 'Electrostatic Potential Average',  &
         'V (eV)', linter)
  else
    call plot_z1D_xmgr(iotape, ave, dave, 0, n3, height,            &
          'pot_ave.agr', 'Electrostatic Potential Average', 'V (eV)')
  endif

  do i=1,n3
    ave(i) = ave(i) / HARTREE
  enddo

! uses convolution for the double average

  do k=1,nwidth

    xave = width(k)/height
    call plot_convol(n3, xave, ave, dave(:,k))

  enddo

  do i=1,n3
    ave(i) = ave(i) * HARTREE
    do k = 1,nwidth
      dave(i,k) = dave(i,k) * HARTREE
    enddo
  enddo


  if(linter) then
    call plot_z1D_gnuplot(ioreplay, iotape, ave, dave, nwidth, n3,    &
          height, 'pot_dave.gp', 'Electrostatic Potential',           &
          'V (eV)', linter)
  else
    call plot_z1D_xmgr(iotape, ave, dave, nwidth, n3, height,            &
          'pot_dave.agr', 'Electrostatic Potential', 'V (eV)')
  endif


! deallocates the stuff


  deallocate(iptype,ipnatom)
  deallocate(nrepeat)
  deallocate(rleft)


  deallocate(ave)
  deallocate(dave)
  deallocate(gave)

  deallocate(rho)
  deallocate(rhogau)

  deallocate(gtmp)
  deallocate(izval)
  deallocate(width)
  deallocate(widthold)
  deallocate(widthgeom)
  deallocate(xpeak)

  deallocate(indx)

  deallocate(autocorr)
  deallocate(conv)
  deallocate(convtmp)


  return

end subroutine plot_rho_v_average
