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

!      calculates and prints the oscillator strengths

       subroutine oscillator_strength(neig,ei,dh0drk,adot,ztot,mxdbnd)

!      Written July 2, 2014. JLM
!      Modified, 4 March 2020, documentation. JLM
!      copyright  Jose Luis Martins/INESC-MN

!      version 4.95


       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  mxdbnd                     !<  array dimension for number of bands

       integer, intent(in)                ::  neig                       !<  number of wavefunctions
       real(REAL64), intent(in)           ::  ei(mxdbnd)                 !<  eigenvalues (Hartree) for rk0
       complex(REAL64), intent(in)        ::  dh0drk(mxdbnd,mxdbnd,3)    !<  d <Psi|H|Psi> d k
       real(REAL64), intent(in)           ::  adot(3,3)                  !<  metric in real space
       real(REAL64), intent(in)           ::  ztot                       !<  total charge density (electrons/cell)

!      local allocatable arrays

       real(REAL64), allocatable          ::  efmei(:)                   !  E_f - E_i
       complex(REAL64), allocatable       ::  ematif(:)                  !  |<i|d H / d k|f>|**2
       complex(REAL64), allocatable       ::  ematcar(:,:)               !  |<i|d H / d k_x|f>|**2 cartesian

!      local variables

       integer           ::  ncond,nval
       real(REAL64)      ::  avec(3,3)                  !  primitive lattice vectors that generate adot in canonical orientation
       real(REAL64)      ::  bvec(3,3)                  !  reciprocal lattice vectors

!      constants

       real(REAL64), parameter  :: ZERO = 0.0_REAL64
       real(REAL64), parameter :: PI = 3.14159265358979323846_REAL64
       real(REAL64), parameter :: HARTREE = 27.21138386_REAL64
       complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)
!      counters

       integer    ::  i, j, k, ij, n, m


       call adot_to_avec_sym(adot,avec,bvec)
       
       write(6,*)
       write(6,*) '  Dipole matrix elements and Oscillator strengths '
       write(6,*)
       write(6,*) '  The cartesian matrix elements are defined'
       write(6,*) '  with respect to the following orientation of'
       write(6,*) '  the crystal primitive lattice vectors'
       write(6,*)
       write(6,'("  a1 = (",f9.4,")     a2 = (",f9.4,")     a3 = (",     &
     &                       f9.4,")")') (avec(1,k),k=1,3)
       write(6,'("       (",f9.4,")          (",f9.4,")          (",     &
     &                       f9.4,")")') (avec(2,k),k=1,3)
       write(6,'("       (",f9.4,")          (",f9.4,")          (",     &
     &                       f9.4,")")') (avec(3,k),k=1,3)

       nval = nint(0.5*ztot)
       ncond = neig - nval

       allocate(efmei(nval*ncond))
       allocate(ematif(nval*ncond))
       allocate(ematcar(nval*ncond,3))

       do i = 1,nval
       do j = 1,ncond
         ij = (i-1)*ncond + j
         efmei(ij) = ei(nval+j) - ei(i)
         ematif(ij) = C_ZERO
         do n = 1,3
         do m = 1,3
           ematif(ij) = ematif(ij) + dh0drk(nval+j,i,n)*     &
     &                          adot(n,m)*dh0drk(i,nval+j,m)
         enddo
         enddo
         ematif(ij) = ematif(ij) / (4*PI*PI)
         do k = 1,3
           ematcar(ij,k) = C_ZERO
         enddo
         do n = 1,3
         do k = 1,3
           ematcar(ij,k) = ematcar(ij,k) + avec(k,n)*dh0drk(nval+j,i,n)
         enddo
         enddo
         do k = 1,3
           ematcar(ij,k) = ematcar(ij,k) / (2*PI)
         enddo
       enddo
       enddo

       write(6,*)
       write(6,'(4x,"i",3x,"E_i(eV)",4x,"f",3x,"E_f(eV)",5x,             &
     &    "|<i|dH/dk|f>|**2",8x,"F_if",9x,"|M_x|**2",4x,                 &
     &    "|M_y|**2",4x,"|M_z|**2")')
       write(6,*)

       do i = 1,nval
       do j = 1,ncond
         ij = (i-1)*ncond + j
         write(6,'(i5,f10.3,i5,f10.3,6x,f12.5,4x,f12.5,4x,3f12.5)')      &
     &       i,ei(i)*HARTREE,nval+j,ei(nval+j)*HARTREE,real(ematif(ij)), &
     &              (2.0/3.0)*real(ematif(ij))/efmei(ij),                &
     &         (abs(ematcar(ij,k))**2,k=1,3)
       enddo
       write(6,*)
       enddo
       write(6,*)

       deallocate(efmei)
       deallocate(ematif)
       deallocate(ematcar)

       return

       end subroutine oscillator_strength
