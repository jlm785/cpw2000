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

!>     interface for minimization with the lbfgs subroutine
!>     adapted from dav_call

!>     WARNING   WARNING   WARNING  :

!>     This is a reverse communication interface. It should be initialized with
!>     istmin = 0. In subsequent calls only istmin,energy,rat,force,adot, and stress
!>     can be changed.  Changes to other variables will lead to disaster.
!>     The subroutine does not check for those changes.

       subroutine move_lbfgs_call(energy,rat,force,adot,stress,         &
     & strext,press,istmin,iconv,minrat,minstr,                         &
     & ntype,natom,                                                     &
     & mkeep,pgtol,dxmax,                                               &
     & mxdtyp,mxdatm)

!      adapted 29 september 2006. JLM
!      modified to include epitaxial relaxation  January 2012. JLM
!      Modified to stabilize minstr=2. 8 June 2014. JLM
!      Modified, documentation. 8 August 2019. JLM
!      copyright INESC-MN/Jose Luis Martins

!      version 4.94
!      version 1.51 of md

       implicit none
       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input:

       integer, intent(in)             ::  mxdtyp                        !<  array dimension of types of atoms
       integer, intent(in)             ::  mxdatm                        !<  array dimension of number of atoms of a given type
       real(REAL64), intent(in)        ::  energy                        !<  energy in hartree
       real(REAL64), intent(in)        ::  force(3,mxdatm,mxdtyp)        !<  component (in lattice coordinates) of the force of the n-th atom of type i (hartree/au)
       real(REAL64), intent(in)        ::  stress(3,3)                   !<  stress tensor in lattice coordinates * vcell (hartree/au)
       real(REAL64), intent(in)        ::  strext(3,3)                   !<  external (applied) stress tensor in lattice coordinates * vcell (hartree/au)
       real(REAL64), intent(in)        ::  press                         !<  applied pressure (hartree/au)
       integer, intent(in)             ::  istmin                        !<  equal to 0 in first step of minimization
       integer, intent(in)             ::  minrat                        !<  if =1 minimize with respect to atomic positions
       integer, intent(in)             ::  minstr                        !<  if =1 minimize with respect to all adot variables, if =2 minimize with respect to adot(3,3)
       integer, intent(in)             ::  ntype                         !<  number of types of atoms
       integer, intent(in)             ::  natom(mxdtyp)                 !<  number of atoms of type i
       integer, intent(in)             ::  mkeep                         !<  number of iterations remembered by lbfgs
       real(REAL64), intent(in)        ::  pgtol                         !<  tolerance for minimization convergence
       real(REAL64), intent(in)        ::  dxmax                         !<  maximum step size

!      input and output:

       real(REAL64), intent(inout)     ::  adot(3,3)                     !<  metric in direct space (atomic units)
       real(REAL64), intent(inout)     ::  rat(3,mxdatm,mxdtyp)          !<  component (in lattice coordinates) of the position of the n-th atom of type i

!      output:

       integer, intent(out)            ::  iconv                         !<  status of minimization: success=0, failure=-1,tryagain=1

!      scratch

       real(REAL64),allocatable, save  :: xdav(:)
       real(REAL64),allocatable, save  :: gdav(:)
       real(REAL64),allocatable, save  :: wa(:)
       integer,allocatable, save       :: iwa(:) 
      
!      constants

       real(REAL64), parameter :: PI = 3.14159265358979323846_REAL64
       real(REAL64), parameter :: TWOPI = 2.0_REAL64*PI       
       real(REAL64), parameter :: ZERO = 0.0_REAL64, DOIS = 2.0_REAL64      

!      local variables

       real(REAL64)         :: adinv(3,3),vcell,entpy,fcov(3)
       integer              ::  ndav
       real(REAL64)         ::  strimb(3,3),temp(3,3)
       real(REAL64)         :: adotnew(3,3)
       logical              :: lmetric
       real(REAL64), save   :: cntrol
       real(REAL64), save   :: a0(3,3),b0(3,3),ad33

!      counters

       integer jc, nt, i, j, k, l

!      allocates work space and saves information about first iteration
       
       
       if(istmin == 0) then

         ad33 = adot(3,3)
         call adot_to_avec(adot,a0,b0)
         do i=1,3
         do j=1,3
           b0(j,i) = b0(j,i) / TWOPI
         enddo
         enddo
         
         if(allocated(xdav)) deallocate(xdav)
         if(allocated(gdav)) deallocate(gdav)
         allocate(xdav(3*mxdatm*mxdtyp+9))
         allocate(gdav(3*mxdatm*mxdtyp+9))

         if(allocated(wa)) deallocate(wa)
         if(allocated(iwa)) deallocate(iwa)
         allocate(wa( (2*mkeep+9)*(3*mxdtyp*mxdatm+9)                   & 
     &                    + 11*mkeep*mkeep+8*mkeep ))
         allocate(iwa(4*(3*mxdtyp*mxdatm+9)))

       endif


       call adot_to_bdot(adot,vcell,adinv)

       do i=1,3
       do j=1,3
         adinv(i,j) = adinv(i,j) / (TWOPI*TWOPI)
       enddo
       enddo


!      converts to lbfgs format in quasi-cartesian coordinates
!      notice that a0 from first iteration is used for all iterations

!      enthalpy (no external stress or pressure for epilbf)

       if (minstr == 2) then

         entpy = energy

       else

         entpy = energy + press*vcell  +                                 &
     &          (strext(1,1)*adot(1,1) +                                 &
     &           strext(1,2)*adot(2,1) +                                 &
     &           strext(1,3)*adot(3,1) +                                 &
     &           strext(2,1)*adot(1,2) +                                 &
     &           strext(2,2)*adot(2,2) +                                 &
     &           strext(2,3)*adot(3,2) +                                 &
     &           strext(3,1)*adot(1,3) +                                 &
     &           strext(3,2)*adot(2,3) +                                 &
     &           strext(3,3)*adot(3,3)) / DOIS

       endif

!      atomic coordinates

       jc = 0
       if(minrat == 1) then

         do nt=1,ntype
           do j=1,natom(nt)

             do k=1,3
               xdav(jc+k) = ZERO
               do l=1,3
                 xdav(jc+k) = xdav(jc+k) + a0(k,l)*rat(l,j,nt)
               enddo
             enddo

             do k=1,3
               fcov(k) = ZERO
               do l=1,3
                 fcov(k) = fcov(k) + adot(k,l)*force(l,j,nt)
               enddo
             enddo

             do k=1,3
               gdav(jc+k) = ZERO
               do l=1,3
                 gdav(jc+k) = gdav(jc+k) + b0(k,l)*fcov(l)
               enddo
             enddo

             jc = jc+3
           enddo
         enddo

       endif

!      lattice parameters

       if(minstr == 1) then

         do i=1,3
         do j=1,3
           strimb(i,j) = strext(i,j)                                    &
     &                 + press*vcell*adinv(i,j)-stress(i,j)
         enddo
         enddo

         do i=1,3
         do j=1,3
           jc = jc + 1
           xdav(jc) = ZERO
           do k=1,3
           do l=1,3
             xdav(jc) = xdav(jc) + b0(i,k)*adot(k,l)*b0(j,l)
           enddo
           enddo
           gdav(jc) = ZERO
           do k=1,3
           do l=1,3
             gdav(jc) = gdav(jc) - a0(i,k)*strimb(k,l)*a0(j,l)
           enddo
           enddo
         enddo
         enddo

       elseif(minstr == 2) then

         jc = jc + 1
         xdav(jc) = adot(3,3)/ad33
         gdav(jc) = stress(3,3)*ad33

       endif
    
       ndav = jc
       do i=1,ndav
         wa(i) = dxmax
         wa(ndav+i) = pgtol
       enddo


!       write(6,*) 'istmin', istmin

       if (istmin == 0) cntrol = ZERO

       call lbfgs(ndav,mkeep,xdav,gdav,entpy,                            &
     &                wa(1),wa(ndav+1),cntrol,wa(2*ndav+1),iwa)

       iconv = -1
       if(nint(cntrol) == 0) iconv=1
       if(nint(cntrol) == 1) iconv=0

!      converts back

!      atomic coordinates

       jc = 0
       if(minrat == 1) then

         do nt=1,ntype
           do j=1,natom(nt)

             do k=1,3
               rat(k,j,nt) = ZERO 
             enddo

             do k=1,3
               jc = jc + 1
               rat(1,j,nt) = rat(1,j,nt) + b0(k,1)*xdav(jc)
               rat(2,j,nt) = rat(2,j,nt) + b0(k,2)*xdav(jc)
               rat(3,j,nt) = rat(3,j,nt) + b0(k,3)*xdav(jc)
             enddo

           enddo
         enddo

       endif

!      lattice parameters

       if(minstr == 1) then

         do i=1,3
         do j=1,3
           jc = jc + 1
           temp(i,j) = xdav(jc)
         enddo
         enddo
         do i=1,3
         do j=1,3
           adotnew(i,j) = ZERO
           do k=1,3
           do l=1,3
             adotnew(i,j) = adotnew(i,j) + a0(k,i)*temp(k,l)*a0(l,j)
           enddo
           enddo
         enddo
         enddo

         adotnew(1,2) = (adotnew(1,2) + adotnew(2,1)) / DOIS
         adotnew(2,1) = adotnew(1,2)
         adotnew(2,3) = (adotnew(2,3) + adotnew(3,2)) / DOIS
         adotnew(3,2) = adotnew(2,3)
         adotnew(3,1) = (adotnew(3,1) + adotnew(1,3)) / DOIS
         adotnew(1,3) = adotnew(3,1)

         call metric_test_g(adotnew,lmetric)
         
         if(lmetric) then
           do i=1,3
           do j=1,3
             adot(i,j) = adotnew(i,j)
           enddo
           enddo
         else
           iconv = -1
         endif

       elseif(minstr == 2) then

         do i=1,3
         do j=1,3
           adotnew(i,j) = adot(i,j)
         enddo
         enddo
         jc = jc + 1
         adotnew(3,3) = xdav(jc)*ad33

         call metric_test_g(adotnew,lmetric)
         
         if(lmetric) then
           adot(3,3) = adotnew(3,3)
         else
           iconv = -1
         endif

       endif

       return
       end subroutine move_lbfgs_call
