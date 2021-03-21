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

!>  Calculates the layer average of a scalar periodic quantity.
!>  parallel to the z-direction.  Includes the atom positions...

subroutine plot_average_simple(ioreplay, func,                           &
         adot, ntype, natom, nameat, rat, zv,                            &
         ng, kgv,                                                        &
         mxdtyp, mxdatm, mxdgve)

! writen September 5, 2012.jlm
! Modified, split code, complex variables, 26 May 2014. JLM
! Documentation, merged psi with rho_v. 5 February 2021. JLM
! copyright  Jose Luis Martins/INESC-MN

  implicit none

! version 4.99

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                :: ioreplay                         !<  tape number for reproducing calculations

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors

  complex(REAL64), intent(in)        ::  func(mxdgve)                    !<  function to be plotted for g-vector j

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space
  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i
  character(len=2), intent(in)       ::  nameat(mxdtyp)                  !<  chemical symbol for the type i
  real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)            !<  k-th component (in lattice coordinates) of the position of the n-th atom of type i
  real(REAL64), intent(in)           ::  zv(mxdtyp)                      !<  valence of atom with type i

  integer, intent(in)                ::  ng                              !<  total number of g-vectors with length less than gmax
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length

! allocatable arrays

  complex(REAL64), allocatable       ::  rhogau(:)                  !  atom centered gaussian charge density for G-vector j
  integer, allocatable               ::  izval(:)                   !  valence of atom of type i

  real(REAL64), allocatable          ::  ave(:)                     !  layer average of function in the 3d direction of fft grid
  real(REAL64), allocatable          ::  dave(:)                    !  double average of function in the 3d direction of fft grid
  real(REAL64), allocatable          ::  gave(:)                    !  average of the nuclear "gaussian" charge density in the 3d direction of fft grid

! main variables

  integer                       ::  nplane                       !  number of lattice planes

  logical                       ::  linter                       !  if .TRUE. asks if the figure should be shown 
  character(len=1)              ::  yesno
  real(REAL64)                  ::  sigma                        !  width of gaussian (in real space)

! parameters that control the quality of the plots 

  integer                       ::  nmult                        !  increases point density in z direction 
  real(REAL64)                  ::  sigref                       !  controls the width of the core gaussian

! other variables

  integer             ::  kmscr(3), nsfft(3)
  real(REAL64)        ::  height
  integer             ::  id, n1, n2, n3, nn

  integer             ::  iotape
  character(len=40)   ::  filename
  integer             ::  mfft, mwrk
  integer             ::  ktmp(3)
  
! constants

  real(REAL64), parameter  :: BOHR = 0.5291772109

! counters

  integer      :: i, j, nt

! These constants may be changed to get nicer plots.

  nmult = 3
  sigref = 0.5

! finds basic fft grid 

  write(6,*) 
  write(6,'("  Do you want to see the plots interactively?   (y/n)")') 
  write(6,*)

  read(5,*) yesno
  write(ioreplay,'(2x,a1,"   interactive plot")') yesno

  if(yesno == 'y' .or. yesno == 'Y') then

    linter = .TRUE.

    write(6,*) 
    write(6,'("  Program will generate files for later plotting with gnuplot")')
    write(6,'("  BEWARE: plots may hide below each other")') 
    write(6,*)

  else

    linter = .FALSE.

    write(6,*) 
    write(6,'("  Program will generate files for later plotting with xmgrace")')
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
  write(6,*)

  nplane = 1

! rewrite code if you want to enable other direction!!!!
! As it is it averages in the plane of first two lattice vectors,
! that is, in the direction of the third reciprocal lattice vector.


  height = adot(3,3) - adot(2,3)*adot(3,2)/adot(2,2)                     &
                     - adot(1,3)*adot(3,1)/adot(1,1)                     &
           + 2*adot(1,2)*adot(2,3)*adot(3,1)/(adot(2,2)*adot(1,1))
  height = sqrt(height)

  write(6,*)
  write(6,'("  the height of the cell is: ",f12.4," Ang.")') height*BOHR
  write(6,*)

! make it more dense on the averaging direction

  ktmp(1) = kmscr(1)
  ktmp(2) = kmscr(2)
  ktmp(3) = nmult*kmscr(3)

! tries to get a multiple of nplane (if nplane is a prime number it just increases the number of points)

  do i = 1,5
    call size_fft(ktmp, nsfft, mfft, mwrk)
    
    if(mod(nsfft(3),nplane) == 0) exit

    ktmp(3) = nsfft(3) / 2
  enddo

!  kave = nsfft(3) / nplane
 
  write(6,*) 
  write(6,'("  the new fft grid is: ",3i8)') (nsfft(j),j=1,3)
  write(6,*) 

  allocate(rhogau(ng))
  allocate(izval(mxdtyp))


! choice of width for gaussian broadening

! sigma = 0.02*height/nplane
  sigma = sigref
  
  n1 = nsfft(1)
  n2 = nsfft(2)
  n3 = nsfft(3)

  nn = nsfft(3)

  id = n1 + 1
!  ntot = id * n2 * n3

! initialize function arrays

  allocate(ave(nn))
  allocate(dave(nn))
  allocate(gave(nn))
  
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


    call plot_zave1D(gave, rhogau, nplane, id,n1,n2,n3, ng, kgv)

    if(linter) then
      filename = 'rho_gauss_' // adjustl(trim(nameat(nt)))//'.gp'
      call plot_z1D_gnuplot(ioreplay, iotape, gave, dave, 0, n3, height, &
              adjustl(trim(filename)),                                   &
              'Broadened Nuclear Density '//nameat(nt),                  &
              '{/Symbol r} (1/cell)',linter)
    else
      filename = 'rho_gauss_'// adjustl(trim(nameat(nt)))//'.agr'
      call plot_z1D_xmgr(iotape, gave, dave, 0, n3, height,              &
              adjustl(trim(filename)),                                   &
             'Broadened Nuclear Density '//nameat(nt),                   &
             '\f{Symbol} r\f{} (1/cell)')
    endif

  enddo



  call plot_zave1D(ave, func, nplane, id,n1,n2,n3, ng, kgv)

  if(linter) then
    call plot_z1D_gnuplot(ioreplay, iotape, ave, dave, 0, n3, height,    &
            'func_ave.gp','Function Average',                            &
            '{/Symbol r} (1/cell)', linter)
  else
    call plot_z1D_xmgr(11, ave, dave, 0, n3, height,'func_ave.agr',      &
           'Function Average', '\f{Symbol} r\f{} (1/cell)')
  endif


! deallocates the stuff

  deallocate(ave)
  deallocate(dave)
  deallocate(gave)

  deallocate(rhogau)

  deallocate(izval)
  
  return

  end subroutine plot_average_simple
