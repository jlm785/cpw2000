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

!>  Reads the Fourier pseudo potentials and the
!>  atomic core and valence charge densities from
!>  files whose name depend on the chemical symbol.
!>  After writes it to iotape which must be opened.
!>  Reads both the old format and the 2026 format (more than one
!>  atomic basis set and core kinetic energy density).  The output
!>  is always in the new format.
!>
!>  \author       Jose Luis Martins
!>  \version      5.13
!>  \date         January 30 2008, 7 October 2026.
!>  \copyright    GNU Public License v2

  subroutine read_write_pseudo(iotape, ntype, nameat,                 &
       pseudo_path, pseudo_suffix, itape_pseudo,                      &
       mxdtyp)

! Written April 16, 2014. jlm
! adapted from Sverre Froyen plane wave program
! adapted from version 4.36 of pseukb.
! Written January 30 2008. jlm
! modified for f90, 16 June 2012. jlm
! style modifications, 7 January 2014. jlm
! modified, vkb dimensions, March 31, 2014. jlm
! modified, so pseudos for non-so file. April 12 2014. JLM
! modified, documentation, August 2019.
! Modified, ititle -> psdtitle, indentation. 20 February 2025. JLM
! Modified, filenames for pseudos. 10 October 2025. JLM
! Modified, several basis sets and core tau (2026 format), new output format. 7 October 2026. JLM+claude



  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms

  integer, intent(in)                ::  iotape                          !<  number of tape to which the pseudo is added.

  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  character(len=2), intent(in)       ::  nameat(mxdtyp)                  !<  chemical symbol for the type i

  character(len=200), intent(in)     ::  pseudo_path                     !<  path to pseudopotentials
  character(len=50), intent(in)      ::  pseudo_suffix                   !<  suffix for the pseudopotentials
  integer, intent(in)                ::  itape_pseudo                    !<  tape number to read pseudo

! local allocatable arrays

  real(REAL64), allocatable   ::  vloc(:)                                !  local pseudopotential for atom k (hartree)
  real(REAL64), allocatable   ::  vkbraw(:,:)                            !  (1/q**l) * kb nonlocal pseudo. for atom k, ang. mom. l. (non normalized to vcell, hartree)
  real(REAL64), allocatable   ::  dcor(:)                                !  core charge density for atom k
  real(REAL64), allocatable   ::  dval(:)                                !  valence charge density for atom k
  real(REAL64), allocatable   ::  wvfraw(:)                              !  wavefunction for atom k, ang. mom. l
  real(REAL64), allocatable   ::  tauc(:)                                !  core kinetic energy density for atom k

! local variables

  character(len=255)       ::  fnam
  integer                  ::  it                                        !  tape number
  character(len=2)         ::  namel, icorrt
  character(len=3)         ::  irel
  character(len=4)         ::  icore
  character(len=60)        ::  iray
  character(len=10)        ::  psdtitle(20)
  integer                  ::  izv, nql, nqnl
  integer                  ::  norb(-1:1), lo(4,-1:1)
  real(REAL64)             ::  delql, vql0
  integer                  ::  nkb(0:3,-1:1)                             !   kb pseudo.  normalization for atom k, ang. mom. l
  real(REAL64)             ::  eorb(0:3,-1:1)
  real(REAL64)             ::  eorbwv
  integer                  ::  nqwf                                      !  number of points for wavefunction interpolation for atom k
  real(REAL64)             ::  delqwf                                    !  step used in the wavefunction interpolation for atom k
  integer                  ::  norbat                                    !  number of atomic orbitals for atom k
  integer                  ::  lorb                                      !  angular momentum of orbital n of atom k
  integer                  ::  n_bsets                                   !  number of atomic basis sets for atom k
  logical                  ::  l2026                                     !  extended pseudo format from 2026
  character(len=512)       ::  line

  integer                  ::  ioerror

! constants

  real(REAL64), parameter  ::  ZERO = 0.0_REAL64

! counters

  integer                  ::  nt, n, j, l, nb

! start loop over atomic types

  do nt = 1,ntype

    it = itape_pseudo + nt
    if(it == iotape) it = 50 + nt

!   open file

    call read_pseudo_get_path(nameat(nt), fnam, pseudo_path, pseudo_suffix)

    open(unit=it,file=fnam,status='old',form='formatted')

!   read heading

    psdtitle(1:20) = '          '

    read(it,'(1x,a2,1x,a2,1x,a3,1x,a4,1x,a60,1x,20a10)',iostat=ioerror)  &
         namel, icorrt, irel, icore, iray, psdtitle

    if(ioerror /= 0) then
      backspace(it)
      read(it,'(1x,a2,1x,a2,1x,a3,1x,a4,1x,a60,1x,7a10)')                &
         namel, icorrt, irel, icore, iray, psdtitle(1:7)
    endif

    write(iotape) namel, icorrt, irel, icore, iray, psdtitle

    read(it,*) izv, nql, delql, vql0
    write(iotape) izv, nql, delql, vql0

!   read pseudopotentials

    nqnl = nql

    if(irel == 'rel') then
      read(it,*) norb(0), norb(-1), norb(1)
      write(iotape) norb(0), norb(-1), norb(1)
    else
      read(it,*) norb(0)
      write(iotape) norb(0)
    endif

    if(irel == 'rel') then
      read(it,*) (lo(j,0),j=1,norb(0)), (lo(j,-1),j=1,norb(-1)),         &
                 (lo(j,1),j=1,norb(1))
      write(iotape) (lo(j,0),j=1,norb(0)), (lo(j,-1),j=1,norb(-1)),      &
                 (lo(j,1),j=1,norb(1))
    else
      read(it,*) (lo(j,0),j=1,norb(0))
      write(iotape) (lo(j,0),j=1,norb(0))
    endif


    if(irel == 'rel') then
      read(it,*) (nkb(lo(j,0),0),j=1,norb(0)),                           &
                 (nkb(lo(j,-1),-1),j=1,norb(-1)),                        &
                 (nkb(lo(j,1),1),j=1,norb(1))
      write(iotape) (nkb(lo(j,0),0),j=1,norb(0)),                        &
                 (nkb(lo(j,-1),-1),j=1,norb(-1)),                        &
                 (nkb(lo(j,1),1),j=1,norb(1))
    else
      read(it,*) (nkb(lo(j,0),0),j=1,norb(0))
      write(iotape) (nkb(lo(j,0),0),j=1,norb(0))
    endif

    if(irel == 'rel') then
      read(it,*) (eorb(lo(j,0),0),j=1,norb(0)),                          &
                 (eorb(lo(j,-1),-1),j=1,norb(-1)),                       &
                 (eorb(lo(j,1),1),j=1,norb(1))
      write(iotape) (eorb(lo(j,0),0),j=1,norb(0)),                       &
                 (eorb(lo(j,-1),-1),j=1,norb(-1)),                       &
                 (eorb(lo(j,1),1),j=1,norb(1))
    else
      read(it,*) (eorb(lo(j,0),0),j=1,norb(0))
      write(iotape) (eorb(lo(j,0),0),j=1,norb(0))
    endif

!   reads the local potential vloc(-1 or 0,nt) should not be used!

    allocate(vloc(nql))
    do j = 1,nql
      read(it,*) vloc(j)
      write(iotape) vloc(j)
    enddo
    deallocate(vloc)

!   reads the non-local pseudopotential

    allocate(vkbraw(0:nqnl,0:1))
    do n = 1,norb(0)
      l = lo(n,0)
      if(irel == 'rel') then
        if(l == 0) then
          do j = 0,nqnl
            read(it,*) vkbraw(j,0)
            write(iotape) vkbraw(j,0)
          enddo
        else
          do j = 0,nqnl
            read(it,*) vkbraw(j,0),vkbraw(j,1)
            write(iotape) vkbraw(j,0),vkbraw(j,1)
          enddo
        endif
      else
        do j = 0,nqnl
          read(it,*) vkbraw(j,0)
          write(iotape) vkbraw(j,0)
        enddo

      endif
    enddo
    deallocate(vkbraw)

!   read core charge

    allocate(dcor(nql))
    do j=1,nql
      read(it,*) dcor(j)
      write(iotape) dcor(j)
    enddo
    deallocate(dcor)

!   read valence charge

    allocate(dval(nql))
    do j=1,nql
      read(it,*) dval(j)
      write(iotape) dval(j)
    enddo
    deallocate(dval)

!   reads the fourier transforms of the wavefunctions.
!   The 2026 format has the number of basis sets in the first line
!   and the number of orbitals in the first line of each basis set.

    line(1:512) = ' '
    read(it,'(a)') line

    read(line,*,iostat=ioerror) nqwf, delqwf, norbat, n_bsets

    l2026 = .TRUE.
    if(ioerror /= 0) then
      read(line,*) nqwf, delqwf, norbat
      n_bsets = 1
      l2026 = .FALSE.
    endif

    write(iotape) nqwf, delqwf, norbat, n_bsets

    allocate(wvfraw(0:nqwf))

    do nb = 1,n_bsets

      if(l2026) read(it,*) lorb, eorbwv, norbat
      write(iotape) norbat

      do n = 1,norbat

        if(.not. l2026 .or. n /= 1) read(it,*) lorb, eorbwv
        write(iotape) lorb, eorbwv

        do j = 0,nqwf-1
          read(it,*) wvfraw(j)
          write(iotape) wvfraw(j)
        enddo

      enddo

    enddo

    deallocate(wvfraw)

!   core kinetic energy density (zero in the old format)

    allocate(tauc(0:nql))

    if(l2026) then
      do j = 0,nql
        read(it,*) tauc(j)
      enddo
    else
      tauc(:) = ZERO
    endif

    write(iotape) (tauc(j),j=0,nql)

    deallocate(tauc)

    close (unit=it)


  enddo

! end loop over atomic types

  return

end subroutine read_write_pseudo
