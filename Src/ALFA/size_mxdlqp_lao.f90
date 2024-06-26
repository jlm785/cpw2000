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

!>     reads the fourier pseudo potentials and determines mxdlqp, mxdlao

       subroutine size_mxdlqp_lao(ntype,nameat,                          &
     &        mxdtyp,mxdlqp,mxdlao)

!      written, 16 June 2012. jlm
!      modified, vkb dimensions, March 31, 2014. jlm
!      Modified, documentation, January 2020. JLM
!      Modified, takes into account that some compilers may insert 
!      line wraps in free format writes, 12 July 2020. JLM
!      Copyright INESC-MN/Jose Luis Martins

!      version 4.983


       implicit none
       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  mxdtyp                     !<  array dimension of types of atoms
       integer, intent(in)                ::  ntype                      !<  number of types of atoms
       character(len=2), intent(in)       ::  nameat(mxdtyp)             !<  chemical symbol for the type i

!      output

       integer, intent(out)               ::  mxdlqp                     !<  array dimension for local potential
       integer, intent(out)               ::  mxdlao                     !<  array dimension of orbital per atom type

!      local variables

       character(len=12)        :: tfile
       character(len=14)        :: fnam
       integer                  :: it                                    !  tape number
       integer                  :: norb(-1:1)
       integer                  :: nql,nqwf,lo(4,-1:1),nkbloc(0:3,-1:1)
       real(REAL64)             :: eorb(0:3,-1:1)
       integer                  :: norbat
       character(len=2)         :: namel,icorrt
       character(len=3)         :: irel

       integer                  :: nt,n,j
       integer                  :: idummy
       real(REAL64)             :: dummy

       tfile = '_POTKB_F.DAT'

!      start loop over atomic types

       mxdlqp = 1
       mxdlao = 1

       do nt=1,ntype
         it = 40 + nt

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

         read(it,'(1x,a2,1x,a2,1x,a3)') namel,icorrt,irel

         read(it,*) idummy,nql,dummy,dummy

         mxdlqp = max(mxdlqp,nql)

!        This is more complicate than it should because of the line wrap in free format!!!

         do n = -1,1
           norb(n) = 0
         enddo

         if(irel == 'rel') then
           read(it,*) norb(0),norb(-1),norb(1)
           if(norb(1) > norb(0) .or. norb(-1) > norb(0)) then
             write(6,'("  WARNING in read_pseudo --- normal orbitals:",  &
     &         i5,"  spin orbitals:",2i5)') norb(0),norb(-1),norb(1)
           endif
         else
           read(it,*) norb(0)
         endif

         if(norb(0) > 4) then
           write(6,'("  stopped in read_pseudo_size  reading data for ",&
     &   a2,"  program cannot accept more than 4 orbitals")') nameat(nt)

           stop

         endif


         if(irel == 'rel') then
           read(it,*) (lo(j,0),j=1,norb(0)),(lo(j,-1),j=1,norb(-1)),     &
     &                (lo(j,1),j=1,norb(1))
         else
           read(it,*) (lo(j,0),j=1,norb(0))
         endif


         if(irel == 'rel') then
           read(it,*) (nkbloc(lo(j,0),0),j=1,norb(0)),                   &
     &                (nkbloc(lo(j,-1),-1),j=1,norb(-1)),                &
     &                (nkbloc(lo(j,1),1),j=1,norb(1))
         else
           read(it,*) (nkbloc(lo(j,0),0),j=1,norb(0))
         endif

         if(irel == 'rel') then
           read(it,*) (eorb(lo(j,0),0),j=1,norb(0)),                     &
     &                (eorb(lo(j,-1),-1),j=1,norb(-1)),                  &
     &                (eorb(lo(j,1),1),j=1,norb(1))
         else
           read(it,*) (eorb(lo(j,0),0),j=1,norb(0))
         endif

!        reads the local potential

         do j = 1,nql
           read(it,*) dummy
         enddo

!        reads the non-local pseudopotential

         do n = 1,norb(0)
           do j = 0,nql
             read(it,*) dummy
           enddo
         enddo

!        read core charge

         do j=1,nql
           read(it,*) dummy
         enddo

!        read valence charge

         do j=1,nql
           read(it,*) dummy
         enddo

!        reads the fourier transforms of the wavefunctions

!        q=0 is for j=1

         read(it,*) nqwf,dummy,norbat

         mxdlqp = max(mxdlqp,nqwf-1)
         mxdlao = max(mxdlao,norbat)

         close(unit=it)

       enddo

       return

       end subroutine size_mxdlqp_lao
