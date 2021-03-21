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

!>     adds a Keating model repulsive term to the LDA energy force and stress.
!>     the parameters are fitted so that the equilibrium geometry, etc. is similar
!>     to the experimental values

       subroutine vff_add_keating(iprint, iowrite,                       &
     & ntype, natom, nameat, rat, adot,                                  &
     & force, energy, stress,                                            &
     & mxdtyp, mxdatm)

!      Modified, documentation, name, details, December 2019. JLM
!      copyright inesc-mn. Jose Luis Martins

!      version 4.94

       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)


!      input:

       integer, intent(in)                :: iprint                      !<  if iprint == 0 does not print details
       integer, intent(in)                :: iowrite                     !<  tape number

       integer, intent(in)                ::  mxdtyp                     !<  array dimension of types of atoms
       integer, intent(in)                ::  mxdatm                     !<  array dimension of types of atoms
       integer, intent(in)                ::  ntype                      !<  number of types of atoms
       integer, intent(in)                ::  natom(mxdtyp)              !<  number of atoms of type i
       character(len=2)                   ::  nameat(mxdtyp)             !<  chemical symbol for the type i
       real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)       !<  k-th component (in lattice coordinates) of the position of the n-th atom of type i
       real(REAL64), intent(in)           ::  adot(3,3)                  !<  metric in direct space

!      input and output:

       real(REAL64), intent(inout)        ::  force(3,mxdatm,mxdtyp)     !<  k-th component (in contravariant lattice coordinates)
       real(REAL64), intent(inout)        ::  energy                     !<  energy
       real(REAL64), intent(inout)        ::  stress(3,3)                !<  stress tensor (in contravariant lattice coordinates)

!      local variables

       real(REAL64), allocatable          ::  force_keat(:,:,:)          !  k-th component (in contravariant lattice coordinates) of the keating force of the n-th atom of type i

       real(REAL64)                       ::  energy_keat                !  keating energy
       real(REAL64)                       ::  stress_keat(3,3)           !  Keating stress tensor (in contravariant lattice coordinates)

!      counters 

       integer nt1,j1,j2
       

       allocate(force_keat(3,mxdatm,mxdtyp))
            
       call vff_keating(iprint, iowrite,                                 &
     &   ntype, natom, nameat, rat, adot,                                &
     &   force_keat, energy_keat, stress_keat,                           &        
     &   mxdtyp, mxdatm)
            
!      call born_mayer(ntype,natom,rat,adot,                              &
!     & force_bm,energy_bm,stress_bm,                                     &
!     & abm,bbm,xrange,                                                   &
!     & mxdtyp,mxdatm)
      
       energy = energy + energy_keat
       
       do nt1=1,ntype
       do j1=1,natom(nt1)
         force(1,j1,nt1) = force(1,j1,nt1) + force_keat(1,j1,nt1)
         force(2,j1,nt1) = force(2,j1,nt1) + force_keat(2,j1,nt1)
         force(3,j1,nt1) = force(3,j1,nt1) + force_keat(3,j1,nt1)
       enddo
       enddo

       do j1=1,3
       do j2=1,3
         stress(j2,j1) = stress(j2,j1) + stress_keat(j2,j1)
       enddo
       enddo
       
       deallocate(force_keat)    
       
       end subroutine vff_add_keating
