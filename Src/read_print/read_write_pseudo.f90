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

!>     reads the fourier pseudo potentials and the
!>     atomic core and valence charge densities from
!>     files whose name depend on the chemical symbol.
!>     and then writes it to iotape which must be opened.

       subroutine read_write_pseudo(iotape,ntype,nameat,                 &
     & mxdtyp)

!      Written April 16, 2014. jlm
!      adapted from Sverre Froyen plane wave program
!      adapted from version 4.36 of pseukb.
!      Written January 30 2008. jlm
!      modified for f90, 16 June 2012. jlm
!      style modifications, 7 January 2014. jlm
!      modified, vkb dimensions, March 31, 2014. jlm
!      modified, so pseudos for non-so file. April 12 2014. JLM
!      modified, documentation, August 2019.
!      Copyright INESC-MN/Jose Luis Martins

!      version 4.94



       implicit none
       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  mxdtyp                     !<  array dimension of types of atoms

       integer, intent(in)                ::  iotape                     !<  number of tape to which the pseudo is added.

       integer, intent(in)                ::  ntype                      !<  number of types of atoms
       character(len=2), intent(in)       ::  nameat(mxdtyp)             !<  chemical symbol for the type i


!      local allocatable arrays

       real(REAL64), allocatable   ::  vloc(:)       !  local pseudopotential for atom k (hartree)
       real(REAL64), allocatable   ::  vkbraw(:,:)   !  (1/q**l) * kb nonlocal pseudo. for atom k, ang. mom. l. (non normalized to vcell, hartree)
       real(REAL64), allocatable   ::  dcor(:)       !  core charge density for atom k
       real(REAL64), allocatable   ::  dval(:)       !  valence charge density for atom k
       real(REAL64), allocatable   ::  wvfraw(:)     !  wavefunction for atom k, ang. mom. l

!      local variables

       character(len=12)        :: tfile
       character(len=14)        :: fnam
       integer                  :: it                !  tape number
       character(len=2)         :: namel, icorrt
       character(len=3)         :: irel
       character(len=4)         :: icore
       character(len=60)        :: iray
       character(len=70)        :: ititle
       integer                  :: izv, nql, nqnl
       integer                  :: norb(-1:1), lo(4,-1:1)
       real(REAL64)             :: delql, vql0
       integer                  :: nkb(0:3,-1:1)       !   kb pseudo.  normalization for atom k, ang. mom. l
       real(REAL64)             :: eorb(0:3,-1:1)
       real(REAL64)             :: eorbwv
       integer                  ::  nqwf               !  number of points for wavefunction interpolation for atom k
       real(REAL64)             ::  delqwf             !  step used in the wavefunction interpolation for atom k
       integer                  ::  norbat             !  number of atomic orbitals for atom k
       integer                  ::  lorb            !  angular momentum of orbital n of atom k

!      counters

       integer                  ::  nt, n, j, l

!      start loop over atomic types

       do nt=1,ntype
 
         it = 40 + nt
         if(it == iotape) it = 50 + nt
         tfile = '_POTKB_F.DAT'

!        open file

         if(nameat(nt)(1:1) /= ' ' .and. nameat(nt)(2:2) /= ' ') then
           fnam = nameat(nt)//tfile
         else if(nameat(nt)(1:1) == ' ') then
           fnam = nameat(nt)(2:2)//tfile
         else
           fnam = nameat(nt)(1:1)//tfile
         endif

         open(unit=it,file=fnam,status='old',form='formatted')

!        read heading

         read(it,'(1x,a2,1x,a2,1x,a3,1x,a4,1x,a60,1x,a70)')             &
     &        namel,icorrt,irel,icore,iray,ititle
         write(iotape) namel,icorrt,irel,icore,iray,ititle

         read(it,*) izv,nql,delql,vql0
         write(iotape) izv,nql,delql,vql0

!        read pseudopotentials

         nqnl = nql

         if(irel == 'rel') then
           read(it,*) norb(0),norb(-1),norb(1)
           write(iotape) norb(0),norb(-1),norb(1)
         else
           read(it,*) norb(0)
           write(iotape) norb(0)
         endif

         if(irel == 'rel') then
           read(it,*) (lo(j,0),j=1,norb(0)),(lo(j,-1),j=1,norb(-1)),     &
     &                (lo(j,1),j=1,norb(1))
           write(iotape) (lo(j,0),j=1,norb(0)),(lo(j,-1),j=1,norb(-1)),  &
     &                (lo(j,1),j=1,norb(1))
         else
           read(it,*) (lo(j,0),j=1,norb(0))
           write(iotape) (lo(j,0),j=1,norb(0))
         endif


         if(irel == 'rel') then
           read(it,*) (nkb(lo(j,0),0),j=1,norb(0)),                      &
     &                (nkb(lo(j,-1),-1),j=1,norb(-1)),                   &
     &                (nkb(lo(j,1),1),j=1,norb(1))
           write(iotape) (nkb(lo(j,0),0),j=1,norb(0)),                   &
     &                (nkb(lo(j,-1),-1),j=1,norb(-1)),                   &
     &                (nkb(lo(j,1),1),j=1,norb(1))
         else
           read(it,*) (nkb(lo(j,0),0),j=1,norb(0))
           write(iotape) (nkb(lo(j,0),0),j=1,norb(0))
         endif

         if(irel == 'rel') then
           read(it,*) (eorb(lo(j,0),0),j=1,norb(0)),                    &
     &                (eorb(lo(j,-1),-1),j=1,norb(-1)),                 &
     &                (eorb(lo(j,1),1),j=1,norb(1))
           write(iotape) (eorb(lo(j,0),0),j=1,norb(0)),                 &
     &                (eorb(lo(j,-1),-1),j=1,norb(-1)),                 &
     &                (eorb(lo(j,1),1),j=1,norb(1))
         else
           read(it,*) (eorb(lo(j,0),0),j=1,norb(0))
           write(iotape) (eorb(lo(j,0),0),j=1,norb(0))
         endif

!        reads the local potential vloc(-1 or 0,nt) should not be used!

         allocate(vloc(nql))
         do j = 1,nql
           read(it,*) vloc(j)
           write(iotape) vloc(j)
         enddo
         deallocate(vloc)

!        reads the non-local pseudopotential

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

!        read core charge

         allocate(dcor(nql))
         do j=1,nql
           read(it,*) dcor(j)
           write(iotape) dcor(j)
         enddo
         deallocate(dcor)

!        read valence charge

         allocate(dval(nql))
         do j=1,nql
           read(it,*) dval(j)
           write(iotape) dval(j)
         enddo
         deallocate(dval)

!        reads the fourier transforms of the wavefunctions

         read(it,*) nqwf,delqwf,norbat
         write(iotape) nqwf,delqwf,norbat

         allocate(wvfraw(0:nqwf))
         do n=1,norbat
           read(it,*) lorb,eorbwv
           write(iotape) lorb,eorbwv

           do j = 0,nqwf-1
             read(it,*) wvfraw(j)
             write(iotape) wvfraw(j)
           enddo

         enddo
         deallocate(wvfraw)

         close (unit=it)


       enddo

!      end loop over atomic types

       return

       end subroutine read_write_pseudo
