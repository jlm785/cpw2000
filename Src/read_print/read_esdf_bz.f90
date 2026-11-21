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

!>  Reads the plane-wave cutoff and BZ integration mesh from a file
!>  opened by a previous esdf_init
!>
!>  \author       Jose Luis Martins
!>  \version      5.13
!>  \date         16 June 2017. 19 November 2026.
!>  \copyright    GNU Public License v2

subroutine read_esdf_bz(ipr,                                             &
   emax, nbandin, mpgrid, mpshift, lbz)

! Written June 16, 2017. jlm
! Modified, documentation, August 10 2019. JLM
! Documentation, one argument per declaration. 28 September 2026. JLM+claude
! Indentation, 28 September 2026. JLM
! Off-diagonal elements of kgrid_Monkhorst_Pack, mpgrid, mpshift. 19 November 2026. JLM+claude WARNING NEW API



  use esdf

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  integer, intent(in)                ::  ipr                             !<  print control: ipr < 3 no printing, except warnings.


! output


  real(REAL64), intent(out)          ::  emax                            !<  kinetic energy cutoff of plane wave expansion (Hartree).

  integer, intent(out)               ::  mpgrid(3,3)                     !<  supercell of the Monkhorst-Pack mesh,  A_i = sum_j mpgrid(i,j) a_j  (line i of kgrid_Monkhorst_Pack)
  real(REAL64), intent(out)          ::  mpshift(3)                      !<  shift of the mesh in units of the reciprocal vectors of the supercell (usually 0.5)

  integer, intent(out)               ::  nbandin                         !<  target for number of bands
  logical, intent(out)               ::  lbz                             !<  indicates if Brillouin Zone data was successfully read.

! local variables

  integer                    ::  nlines
  integer                    ::  ioerr
  logical                    ::  ldiag
  integer                    ::  idet

! parameters

  real(REAL64), parameter ::  ZERO = 0.0_REAL64

! counters

  integer   ::  i, j


! plane-wave cutoff energy

  lbz = .TRUE.

  if(esdf_defined('PWEnergyCutoff')) then

    emax = ZERO
    emax = esdf_physical('PWEnergyCutoff',emax,'hartree')

    if(emax < 1.0) then
      write(6,*)
      write(6,*) '   WARNING'
      write(6,*)
      write(6,*) '   Plane-wave cut-off looks very small'
      write(6,'("    The value of emax is: ",e12.3)') emax
    endif

    if(emax > 500.0) then
      write(6,*)
      write(6,*) '   WARNING'
      write(6,*)
      write(6,*) '   Plane-wave cut-off looks looks very large'
      write(6,'("   The value of emax is: ",e12.3)') emax
    endif

  else
    lbz = .FALSE.
    emax = 5.0
    write(6,*)
    write(6,*) '   WARNING'
    write(6,*)
    write(6,*) '   Using default values for emax'
  endif


! Number of eigenstates

  if(esdf_defined('NumberOfEigenStates')) then

    nbandin = 0
    nbandin = esdf_integer('NumberOfEigenStates',nbandin)

    if(nbandin < 1) then
      nbandin = 1
      write(6,*)
      write(6,*) '   WARNING'
      write(6,*)
      write(6,*) '   Number of bands was 0 or negative. Set to 1.'
    endif

    if(nbandin > 5000) then
      write(6,*)
      write(6,*) '   WARNING'
      write(6,*)
      write(6,*) '   Number of bands looks very large'
      write(6,'("    The value of nbandin is: ",i8)') nbandin
    endif

  else
    lbz = .FALSE.
    nbandin = 10
    write(6,*)
    write(6,*) '   WARNING'
    write(6,*)
    write(6,*) '   Using default values for nbandin'
  endif


! BZ integration


  if(esdf_block('kgrid_Monkhorst_Pack',nlines)) then
    if(nlines /= 3) then
      write(6,*)
      write(6,*) '  esdf_read_bz:  Wrong Number of Lines'
      write(6,*) '  in kgrid_Monkhorst_Pack!'

      stop

    endif

    ioerr = 0
    do i=1,nlines

      read(block_data(i),*,iostat=ioerr) (mpgrid(i,j),j=1,3),mpshift(i)

      if(ioerr /= 0) then
        write(6,*)
        write(6,*) '    Stopped in read_esdf_bz'
        write(6,*) '    error reading kgrid_Monkhorst_Pack'

        stop

      endif

    enddo

!   line i is the supercell vector i in units of the lattice vectors (as in siesta)

    ldiag = .TRUE.
    do i = 1,3
    do j = 1,3
      if(i /= j .and. mpgrid(i,j) /= 0) ldiag = .FALSE.
    enddo
    enddo

    if(ldiag) then

      do i = 1,3
        if(mpgrid(i,i) < 1) then
          mpgrid(i,i) = 1
          write(6,*)
          write(6,*) '   WARNING'
          write(6,*)
          write(6,'("    diagonal element ",i1," of kgrid_Monkhorst_Pack ",  &
             &      "was 0 or negative. Set to 1.")') i
        endif
      enddo

    else

      idet = mpgrid(1,1)*(mpgrid(2,2)*mpgrid(3,3) - mpgrid(2,3)*mpgrid(3,2))  &
           - mpgrid(1,2)*(mpgrid(2,1)*mpgrid(3,3) - mpgrid(2,3)*mpgrid(3,1))  &
           + mpgrid(1,3)*(mpgrid(2,1)*mpgrid(3,2) - mpgrid(2,2)*mpgrid(3,1))

      if(idet == 0) then
        write(6,*)
        write(6,*) '    Stopped in read_esdf_bz'
        write(6,*) '    the kgrid_Monkhorst_Pack matrix is singular'
        write(6,*)
        write(6,'(3i6)') (mpgrid(1,j),j=1,3)
        write(6,'(3i6)') (mpgrid(2,j),j=1,3)
        write(6,'(3i6)') (mpgrid(3,j),j=1,3)
        write(6,*)

        stop

      endif

    endif

  else
    lbz = .FALSE.

    write(6,*)
    write(6,*) '   WARNING'
    write(6,*)
    write(6,*) '   Using default values for BZ integration'

    do j = 1,3
    do i = 1,3
      mpgrid(i,j) = 0
    enddo
    mpgrid(j,j) = 4
    mpshift(j) = 0.5
    enddo
  endif


! prints the results

  if(ipr > 2) then
    write(6,*)
    write(6,*)
    write(6,*)  '  The values set by read_esdf_bz are:'
    write(6,*)
    write(6,'("   The value of emax is: ",e12.3)') emax
    write(6,'("   The value of nbandin is: ",i8)') nbandin
    write(6,*) '  The kgrid_Monkhorst_Pack supercell and shift are:'
    do i = 1,3
      write(6,'(5x,3i6,f10.3)') (mpgrid(i,j),j=1,3), mpshift(i)
    enddo
    write(6,*)
    write(6,*)
  endif

  return

end subroutine read_esdf_bz
