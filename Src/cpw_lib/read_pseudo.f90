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

!>  reads the Fourier pseudo potentials and the
!>  atomic core and valence charge densities from
!>  files whose name depend on the chemical symbol.
!>  it then computes the same quantities on the prototype
!>  g-vectors.
!>
!>  \author       Jose Luis Martins
!>  \version      5.11
!>  \date         1980s, 20 February 2025.
!>  \copyright    GNU Public License v2

subroutine read_pseudo(ipr, author, ealraw,                              &
      nqnl, delqnl, vkbraw, nkb, vloc, dcor, dval,                       &
      norbat, nqwf, delqwf, wvfao, lorb, latorb,                         &
      ntype, natom, nameat, zv, ztot,                                    &
      mxdtyp, mxdlqp, mxdlao)

! adapted from Sverre Froyen plane wave program
! adapted from version 4.36 of pseukb.
! Written January 30 2008. jlm
! modified for f90, 16 June 2012. jlm
! style modifications, 7 January 2014. jlm
! modified, vkb dimensions, March 31, 2014. jlm
! modified, so pseudos for non-so file. April 12 2014. JLM
! modified, delqwf in wvfao. September 13, 2015. JLM
! Modified, documentation, January 2020. JLM
! Modified, eorbwv, 29 November 2021. JLM
! Modified, polarization orbitals of f not processed. 2 December 2021. JLM
! Modified, Perdew-Wang (1992) not flagged as unsupported. 12 January 2024. JLM
! Modified, ititle -> psdtitle, useless but for consistency. 20 February 2025. JLM


  implicit none
  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdlqp                          !<  array dimension for local potential
  integer, intent(in)                ::  mxdlao                          !<  array dimension of orbital per atom type

  integer, intent(in)                ::  ipr                             !<  should be equal to one if information is to be printed.
  character(len=*), intent(in)       ::  author                          !<  type of correlation
  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i
  character(len=2), intent(in)       ::  nameat(mxdtyp)                  !<  chemical symbol for the type i

! output

  real(REAL64), intent(out)          ::  ealraw                          !<  G=0 contrib. to the total energy. (non norm. to vcell,hartree)
  integer, intent(out)               ::  nqnl(mxdtyp)                    !<  number of points for pseudo interpolation for atom k
  real(REAL64), intent(out)          ::  delqnl(mxdtyp)                  !<  step used in the pseudo interpolation for atom k
  real(REAL64), intent(out)       ::  vkbraw(-2:mxdlqp,0:3,-1:1,mxdtyp)  !<  (1/q**l) * kb nonlocal pseudo. for atom k, ang. mom. l. (non normalized to vcell, hartree)
  integer, intent(out)               ::  nkb(0:3,-1:1,mxdtyp)            !<   kb pseudo.  normalization for atom k, ang. mom. l
  real(REAL64), intent(out)          ::  vloc(-1:mxdlqp,mxdtyp)          !<  local pseudopotential for atom k (hartree)
  real(REAL64), intent(out)          ::  dcor(-1:mxdlqp,mxdtyp)          !<  core charge density for atom k
  real(REAL64), intent(out)          ::  dval(-1:mxdlqp,mxdtyp)          !<  valence charge density for atom k
  integer, intent(out)               ::  norbat(mxdtyp)                  !<  number of atomic orbitals for atom k
  integer, intent(out)               ::  nqwf(mxdtyp)                    !<  number of points for wavefunction interpolation for atom k
  real(REAL64), intent(out)          ::  delqwf(mxdtyp)                  !<  step used in the wavefunction interpolation for atom k
  integer, intent(out)               ::  lorb(mxdlao,mxdtyp)             !<  angular momentum of orbital n of atom k
  real(REAL64), intent(out)          ::  wvfao(-2:mxdlqp,mxdlao,mxdtyp)  !<  (1/q**l) * wavefunction for atom k, ang. mom. l (non normalized to vcell)
  logical, intent(out)               ::  latorb                          !<  indicates if all atoms have information about atomic orbitals
  real(REAL64), intent(out)          ::  zv(mxdtyp)                      !<  valence of atom with type i
  real(REAL64), intent(out)          ::  ztot                            !<  total charge density (electrons/cell)

! vloc,dcor,dval contains the fourier transform for q=i*delq
! vkbraw contains (1/q**l) * the fourier transform for q=i*delq
! wvfao contains (1/q**l) * the fourier transform for q=i*delqwf

! local variables

  character(len=12)        :: tfile
  character(len=14)        :: fnam
  integer                  :: it                !  tape number
  character(len=2)         :: namel, icorr, icorrt
  character(len=3)         :: irel
  character(len=4)         :: icore
  character(len=60)        :: iray
  character(len=10)        :: psdtitle(20)
  integer                  :: izv,nql,norb(-1:1),lo(4,-1:1)
  real(REAL64)             :: ealpha,delql,vql0,fac
  real(REAL64)             :: eorb(0:3,-1:1)
  real(REAL64)             :: eorbwv(0:mxdlao)

  integer                  :: nskip

  integer                  :: ioerror

! constants

  real(REAL64), parameter  :: ZERO = 0.0_REAL64, UM = 1.0_REAL64

! counters

  integer                  :: nt,n,j,m,l,mmax,nc


! write heading

  if (ipr == 1) write(6,'(//," potentials :",/," ------------")')

  latorb = .TRUE.

  icorr = author(1:2)

! start loop over atomic types

  ztot = ZERO
  ealpha = ZERO

  do nt=1,ntype

    it = 40 + nt
    tfile = '_POTKB_F.DAT'

!   open file

    if(nameat(nt)(1:1) /= ' ' .and. nameat(nt)(2:2) /= ' ') then
      fnam = nameat(nt)//tfile
    else if(nameat(nt)(1:1) == ' ') then
      fnam = nameat(nt)(2:2)//tfile
    else
      fnam = nameat(nt)(1:1)//tfile
    endif

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

    read(it,*) izv,nql,delql,vql0
    zv(nt) = izv*UM
    delqnl(nt) =delql

    if (ipr == 1) write(6,'(/,1x,a2,2x,a2,2x,a3,2x,a4,/,1x,a60,/,        &
      &     1x,20a10,/," nql=",i4," delql=",f9.4)')                      &
            namel, icorrt, irel, icore, iray, psdtitle, nql, delql

    if (nql > mxdlqp) then
      write(6,'("  Stopped in read_pseudo  reading data for ",a2,        &
      &      "  increase mxdlqp (potential) to ",i7)') nameat(nt),nql

      stop

    endif

    if (namel /= nameat(nt)) write(6,'("  *** warning in read_",         &
      &     "pseudo   chemical symbols do not match")')

!     if(icorrt == 'ca') icorrt='CA'
!     if(icorrt == 'xa') icorrt='XA'
!     if(icorrt == 'wi') icorrt='WI'
!     if(icorrt == 'hl') icorrt='HL'
!     if(icorrt == 'pb') icorrt='PB'
!     if(icorrt == 'pw') icorrt='PW'

    call chrcap(icorrt,2)
    if(icorrt /= icorr) write(6,'("  *** warning in read_pseudo",        &
      &   "  correlation potential does not match")')

!   sum up vql0 for the alpha energy
!   and zv for the total valence charge

    ealpha = ealpha + vql0*natom(nt)
    ztot = ztot + zv(nt)*natom(nt)

!   read pseudopotentials

    nqnl(nt) = nql

    do m = -1,1
      norb(m) = 0
    enddo
    if(irel == 'rel') then
      read(it,*) norb(0),norb(-1),norb(1)
      if(norb(1) > norb(0) .or. norb(-1) > norb(0)) then
        write(6,'("  WARNING in read_pseudo --- normal orbitals:",       &
          &      i5,"  spin orbitals:",2i5)') norb(0),norb(-1),norb(1)
      endif
    else
      read(it,*) norb(0)
    endif

    if(norb(0) > 4) then
      write(6,'("  stopped in read_pseudo  reading data for ",a2,        &
      &   "  program cannot accept more than 4 orbitals")') nameat(nt)

      stop

    endif

    if(irel == 'rel') then
      read(it,*) (lo(j,0),j=1,norb(0)),(lo(j,-1),j=1,norb(-1)),          &
                 (lo(j,1),j=1,norb(1))
    else
      read(it,*) (lo(j,0),j=1,norb(0))
    endif

    if(irel /= 'rel') then
      norb(1) = norb(0)
      do j=1,norb(0)
        lo(j,1) = lo(j,0)
      enddo
      nc = 0
      do j=1,norb(0)
        if(lo(j,0) /= 0) then
          nc = nc + 1
          lo(nc, -1) = lo(j,0)
        endif
        norb(-1) = nc
      enddo
    endif

    if(irel == 'rel') then
      mmax = 1
    else
      mmax = 0
    endif
    do m = -mmax,mmax
    do j = 1,norb(m)

      if(lo(j,m) > 3) then
        write(6,'("  stopped in read_pseudo  reading data for ",         &
          &    a2,"  program does noes not accept l=",i4,                &
          &    " in potential")')  nameat(nt),lo(j,m)

        stop

      endif
      if(lo(j,m) < 0) then
        write(6,'("  stopped in read_pseudo  reading data for ",         &
          &    a2,"  program does noes not accept negative l in",        &
          &    " potential")') nameat(nt)

        stop

      endif
    enddo
    enddo

    do m = -1,1
    do l = 0,3
      nkb(l,m,nt) = 0
    enddo
    enddo


    if(irel == 'rel') then
      read(it,*) (nkb(lo(j,0),0,nt),j=1,norb(0)),                        &
                 (nkb(lo(j,-1),-1,nt),j=1,norb(-1)),                     &
                 (nkb(lo(j,1),1,nt),j=1,norb(1))
    else
      read(it,*) (nkb(lo(j,0),0,nt),j=1,norb(0))
    endif

    if(irel /= 'rel') then
      do j=1,norb(0)
        nkb(lo(j,1),1,nt) = nkb(lo(j,0),0,nt)
      enddo
      nc = 0
      do j=1,norb(0)
        if(lo(j,0) /= 0) then
          nc = nc + 1
          nkb(lo(nc,-1),-1,nt) = nkb(lo(j,0),0,nt)
        endif
      enddo
    endif


    do m=-1,1
    do j=0,3
      eorb(j,m) = ZERO
    enddo
    enddo
    if(irel == 'rel') then
      read(it,*) (eorb(lo(j,0),0),j=1,norb(0)),                          &
                 (eorb(lo(j,-1),-1),j=1,norb(-1)),                       &
                 (eorb(lo(j,1),1),j=1,norb(1))
    else
      read(it,*) (eorb(lo(j,0),0),j=1,norb(0))
    endif

!   reads the local potential vloc(-1 or 0,nt) should not be used!

    do j = 1,nql
      read(it,*) vloc(j,nt)
    enddo
    vloc(-1,nt) = vloc(1,nt)
    vloc(0,nt) = vloc(1,nt)

!   reads the non-local pseudopotential

    do l = 0,3
    do j = -2,mxdlqp
      vkbraw(j,l,-1,nt) = ZERO
      vkbraw(j,l, 0,nt) = ZERO
      vkbraw(j,l, 1,nt) = ZERO
    enddo
    enddo
    do n = 1,norb(0)
      l = lo(n,0)
      if(irel == 'rel') then
        if(l == 0) then
          do j = 0,nqnl(nt)
            read(it,*) vkbraw(j,l,0,nt)
            vkbraw(j,l,-1,nt) = vkbraw(j,l,0,nt)
            vkbraw(j,l, 1,nt) = vkbraw(j,l,0,nt)
          enddo
        else
          do j = 0,nqnl(nt)
            read(it,*) vkbraw(j,l,0,nt),vkbraw(j,l,1,nt)
            vkbraw(j,l,-1,nt) = vkbraw(j,l,0,nt) - ((l+1)*vkbraw(j,l,1,nt))/2
            vkbraw(j,l, 1,nt) = vkbraw(j,l,0,nt) + (l*vkbraw(j,l,1,nt))/2
          enddo
        endif
      else
        do j = 0,nqnl(nt)
          read(it,*) vkbraw(j,l,0,nt)
          vkbraw(j,l,-1,nt) = vkbraw(j,l,0,nt)
          vkbraw(j,l, 1,nt) = vkbraw(j,l,0,nt)
        enddo

      endif

!     divides by 1/q**l

      if(l > 0) then
        do j = 1,nqnl(nt)
          fac = (j*delqnl(nt))**l
          vkbraw(j,l,-1,nt) = vkbraw(j,l,-1,nt) / fac
          vkbraw(j,l, 0,nt) = vkbraw(j,l, 0,nt) / fac
          vkbraw(j,l, 1,nt) = vkbraw(j,l, 1,nt) / fac
        enddo

        vkbraw(0,l,-1,nt) = (4*vkbraw(1,l,-1,nt) - vkbraw(2,l,-1,nt)) / 3
        vkbraw(0,l, 0,nt) = (4*vkbraw(1,l, 0,nt) - vkbraw(2,l, 0,nt)) / 3
        vkbraw(0,l, 1,nt) = (4*vkbraw(1,l, 1,nt) - vkbraw(2,l, 1,nt)) / 3

      endif

      vkbraw(-1,l,-1,nt) = vkbraw(1,l,-1,nt)
      vkbraw(-1,l, 0,nt) = vkbraw(1,l, 0,nt)
      vkbraw(-1,l, 1,nt) = vkbraw(1,l, 1,nt)
      vkbraw(-2,l,-1,nt) = vkbraw(2,l,-1,nt)
      vkbraw(-2,l, 0,nt) = vkbraw(2,l, 0,nt)
      vkbraw(-2,l, 1,nt) = vkbraw(2,l, 1,nt)

    enddo
!   read core charge

    do j=1,nql
      read(it,*) dcor(j,nt)
    enddo
    dcor(-1,nt) = dcor(1,nt)
    dcor(0,nt) = dcor(1,nt) - (dcor(2,nt)-dcor(1,nt))/3

!   read valence charge

    do j=1,nql
      read(it,*) dval(j,nt)
    enddo
    dval(-1,nt) = dval(1,nt)
    dval(0,nt) = zv(nt)

!   reads the fourier transforms of the wavefunctions

    do n = 1,mxdlao
    do j = -2,mxdlqp
      wvfao(j,n,nt) = ZERO
    enddo
    enddo
    read(it,*) nqwf(nt),delqwf(nt),norbat(nt)

    if(nqwf(nt)-1 > mxdlqp) then
      write(6,'("  stopped in read_pseudo  reading data for ",           &
        &    a2,"  increase mxdlqp (wavefunction) to ",i7)')             &
          nameat(nt),nqwf(nt)-1

      stop

    endif
    if(norbat(nt) > mxdlao) then
      write(6,'("  stopped in read_pseudo  reading data for ",           &
        &    a2,"  program does not accept ",i4,                         &
        &    " orbitals (max",i4,")")') nameat(nt),norbat(nt),mxdlao

      stop

    endif

    nskip = 0
    do n=1,norbat(nt)

      read(it,*) lorb(n-nskip,nt),eorbwv(n-nskip)

      l = lorb(n-nskip,nt)

      if(l < 0) then
        write(6,'("  stopped in read_pseudo  reading data for ",         &
          &    a2,"  program does noes not accept negative l in",        &
          &    " wavefunctions")') nameat(nt)

        stop

      endif

      if(l > 3) then
        write(6,'("  WARNING in read_pseudo  reading data for ",         &
          &    a2,"  code not written for l= ",i4,                       &
          &    " in wavefunctions")') nameat(nt),l
        write(6,'("  skipping atomic orbital")')
        do j = 0,nqwf(nt)-1
          read(it,*) wvfao(j,n-nskip,nt)
        enddo
        nskip = nskip+1
      else
        do j = 0,nqwf(nt)-1
          read(it,*) wvfao(j,n-nskip,nt)
        enddo

!       divides by 1/q**l

        if(l > 0) then
          do j = 1,nqwf(nt)-1
            fac = (j*delqwf(nt))**l
            wvfao(j,n-nskip,nt) = wvfao(j,n-nskip,nt) / fac
          enddo

          wvfao(0,n-nskip,nt) = (4*wvfao(1,n-nskip,nt) - wvfao(2,n-nskip,nt)) / 3

        endif

        wvfao(-1,n-nskip,nt) = wvfao(1,n-nskip,nt)
        wvfao(-2,n-nskip,nt) = wvfao(2,n-nskip,nt)

      endif

    enddo
    norbat(nt) = norbat(nt) - nskip

    close (unit=it)

!   converts to hartree

    do j = -1,nql
      vloc(j,nt) = vloc(j,nt) / 2
    enddo
    do n = 0,3
      do j = -2,nqnl(nt)
        vkbraw(j,n,-1,nt) = vkbraw(j,n,-1,nt) / sqrt(2*UM)
        vkbraw(j,n, 0,nt) = vkbraw(j,n, 0,nt) / sqrt(2*UM)
        vkbraw(j,n, 1,nt) = vkbraw(j,n, 1,nt) / sqrt(2*UM)
      enddo
    enddo

  enddo

! end loop over atomic types

! compute the alpha energy. this is the contribution
! to the total energy coming from G=0.

  ealraw = ztot*ealpha

! converts to hartree

  ealraw = ealraw / 2

  return

end subroutine read_pseudo
