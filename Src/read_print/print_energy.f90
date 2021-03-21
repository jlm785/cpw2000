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

!>     prints several energy terms

       subroutine print_energy(ipr,entype,energy,force,stress,           &
     & adot,ntype,natom,nameat,                                          &
     & mxdtyp,mxdatm)

!      written 15 january 1999. jlm
!      modified for f90, 21 October 2015. JLM
!      Modified, documentation, June 2019. JLM
!      copyright INESC-MN/Jose Luis Martins

!      version 4.94 of cpw
!      version 1.5 of md


       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  mxdtyp                     !<  array dimension of types of atoms
       integer, intent(in)                ::  mxdatm                     !<  array dimension of number of atoms of a given type

       integer, intent(in)                ::  ipr                        !<  should be equal to one if information is to be printed.
       character(len=*),intent(in)        ::  entype                     !<  type of energy contribution
       real(REAL64), intent(in)           ::  adot(3,3)                  !<  metric in direct space
       integer, intent(in)                ::  ntype                      !<  number of types of atoms
       integer, intent(in)                ::  natom(mxdtyp)              !<  number of atoms of type i
       character(len=2), intent(in)       ::  nameat(mxdtyp)             !<  chemical symbol for the type i

       real(REAL64), intent(in)            ::  energy                    !<  energy
       real(REAL64), intent(in)            ::  force(3,mxdatm,mxdtyp)    !<  k-th component (in contravariant lattice coordinates)  of the force of the n-th atom of type i
       real(REAL64), intent(in)            ::  stress(3,3)               !<  stress tensor (in contravariant lattice coordinates)
       

!      local variables

       real(REAL64)       ::  strcar(3,3)
       real(REAL64)       ::  avec(3,3),bvec(3,3),bdot(3,3),vcell
       real(REAL64)       ::  press
       real(REAL64)       ::  car(3)

!      parameters

       real(REAL64), parameter :: ZERO = 0.0_REAL64
       real(REAL64), parameter :: AUTOGPA = 29421.58_REAL64

!      counters

       integer       ::  i,j,i1,i2,j1,j2
       integer       ::  nt,ntt,ja,jmax


       if(ipr /= 1) return

!      writes energy

       write(6,*)
       write(6,*) '  Energy, force, stress from ',entype
       write(6,*)
       write(6,'(5x,f15.8,5x," energy ",a50)') energy,entype

!      writes stress tensor

       call adot_to_avec_sym(adot,avec,bvec)
       call adot_to_bdot(adot,vcell,bdot)

       do i1=1,3
       do i2=1,3
         strcar(i1,i2) = ZERO
         do j1=1,3
         do j2=1,3
           strcar(i1,i2) = strcar(i1,i2) +                               &
     &                 avec(i1,j1) * stress(j1,j2) * avec(i2,j2)
         enddo
         enddo
         strcar(i1,i2) = strcar(i1,i2) * autogpa / vcell
       enddo
       enddo

       write(6,*)
       write(6,'(8x,"Contravariant stress tensor  (a.u.)",15x,           &
     &         "Cartesian stress (GPa)")')
       do i=1,3
         write (6,'(5x,3(2x,f10.6),5x,3(2x,e14.6),"   stress ",i1,a20)') &
     &         (stress(i,j),j=1,3),(strcar(i,j),j=1,3),i,entype
      enddo
      write(6,*)

!      pressure

       press = zero
       do i=1,3
       do j=1,3
         press = press + adot(i,j)*stress(j,i)
       enddo
       enddo
       press = press/(3*vcell)
       write(6,'(5x,f15.8,5x,f15.8,5x," pressure (au and GPa) ",a20)')   &
     &         press,press*autogpa,entype

!      force

       write(6,*)
       write(6,'(10x,"Force (Lattice coord.)",12x,                       &
     &          "Force (Cartesian coord. a.u)",7x,"no. type ")')
        write(6,*)
       ntt = 0
       do nt=1,ntype
         jmax = natom(nt)
         do ja=1,jmax
           ntt = ntt + 1
           car(1) = avec(1,1)*force(1,ja,nt) +                           &
     &              avec(1,2)*force(2,ja,nt) +                           &
     &              avec(1,3)*force(3,ja,nt)
           car(2) = avec(2,1)*force(1,ja,nt) +                           &
     &              avec(2,2)*force(2,ja,nt) +                           &
     &              avec(2,3)*force(3,ja,nt)
           car(3) = avec(3,1)*force(1,ja,nt) +                           &
     &              avec(3,2)*force(2,ja,nt) +                           &
     &              avec(3,3)*force(3,ja,nt)
           write(6,'(2x,3(2x,f9.5),3x,3(1x,e12.5),1x,i3,3x,a2,3x,        &
     &     " force ",a20)') (force(i,ja,nt),i=1,3),                      &
     &                      (car(i),i=1,3),ntt,nameat(nt),entype
         enddo
       enddo

       return
       end subroutine print_energy
