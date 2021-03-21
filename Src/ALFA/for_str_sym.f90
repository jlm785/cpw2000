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

!>     Symmetrizes the Hellman-Feynman forces and stress.
!>     Adapted from Sverre Froyen plane wave program.

       subroutine for_str_sym(funsym,sunsym,fsym,ssym,                   &
     & ntype,natom,rat,                                                  &
     & ntrans,mtrx,tnp,                                                  &
     & mxdtyp,mxdatm)

!      Written January 18 1988. jlm
!      Modified March 12, 1999. jlm
!      Modified December 2016, f90. JLM
!      Modified, documentation, January 2020. JLM
!      copyright inesc-mn/Jose Luis Martins

!      version 4.94

       implicit none

       integer, parameter :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  mxdtyp                     !<  array dimension of types of atoms
       integer, intent(in)                ::  mxdatm                     !<  array dimension of number of atoms of a given type

       real(REAL64), intent(in)           ::  funsym(3,mxdatm,mxdtyp)    !<  unsymmetrized force  (in covariant lattice coordinates)
       real(REAL64), intent(in)           ::  sunsym(3,3)                !<  unsymmetrized stress tensor in covariant lattice coordinates (hartree/bohr) 

       integer, intent(in)                ::  ntype                      !<  number of types of atoms
       integer, intent(in)                ::  natom(mxdtyp)              !<  number of atoms of type i
       real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)       !<  lattice coordinates of atom j of type i

       integer, intent(in)                ::  ntrans                     !<  number of symmetry operations in the factor group
       integer, intent(in)                ::  mtrx(3,3,48)               !<  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the factor group
       real(REAL64), intent(in)           ::  tnp(3,48)                  !<  2*pi* i-th component (in lattice coordinates) of the fractional translation vector associated with the k-th symmetry operation of the factor group

!      output

       real(REAL64), intent(out)          :: fsym(3,mxdatm,mxdtyp)       !<  symmetrized force  (in lattice coordinates)
       real(REAL64), intent(out)          ::  ssym(3,3)                  !<  symmetrized stress tensor in covariant lattice coordinates (hartree/bohr)

!      local variables

       real(REAL64)           ::  cro(3), xdif
       integer                ::  kk, idist, idif

!      constants

       real(REAL64), parameter  :: ZERO = 0.0_REAL64
       real(REAL64), parameter :: PI = 3.14159265358979323846_REAL64
       real(REAL64), parameter :: EPS = 1.0d-6

!      counters

       integer       ::  i, j, k, l, m, nt


!      initializes

       do i = 1,ntype
         do j = 1,natom(i)
         do k = 1,3
           fsym(k,j,i) = ZERO
         enddo
         enddo
       enddo
!         
       do i = 1,3
       do j = 1,3
         ssym(j,i) = ZERO
       enddo
       enddo

!      symmetrizes forces

       do i = 1,ntrans
         do j = 1,ntype
           do k = 1,natom(j)

!                     -1
!            find mtrx    * (rat - tnp)

             do l = 1,3
               cro(l) = ZERO
               do m = 1,3
                 cro(l) = cro(l) + mtrx(m,l,i)*                      &
     &                 (rat(m,k,j) - tnp(m,i)/(2*PI))
               enddo
             enddo

             kk = 0
             do l = 1,natom(j)
               idist = 0
               do m = 1,3
                 xdif = abs(cro(m) - rat(m,l,j))
                 idif = nint(xdif)
                 if (abs(xdif - real(idif,REAL64)) > EPS)                &
     &                 idist = idist + 1
               enddo
               if(idist == 0) then
                 kk = l
               endif
             enddo

             if(kk == 0) then
               write(6,'("   STOPPED in for_str_sym:    sym ",           &
     &            "operation = ",i3," type of atom = ",i5,               &
     &            "  atom = ",i5)') i, j, k

               stop

             endif

!            rotate force and add

             do l = 1,3
             do m = 1,3
               fsym(l,k,j) = fsym(l,k,j)                                 &
     &                           + mtrx(l,m,i) * funsym(m,kk,j) 
             enddo
             enddo

           enddo
         enddo
       enddo

!      symmetrize stress

       do i = 1,3
       do j = 1,3
         do nt = 1,ntrans
           do k = 1,3
           do l = 1,3
             ssym(j,i) = ssym(j,i)                                       &
     &             + sunsym(k,l) * (mtrx(j,k,nt) * mtrx(i,l,nt))
           enddo
           enddo
         enddo
         ssym(j,i) = ssym(j,i) / ntrans
       enddo
       enddo

!      renormalizes

       do i = 1,ntype
         do j = 1,natom(i)
         do k = 1,3
           fsym(k,j,i) = fsym(k,j,i) / ntrans
         enddo
         enddo
       enddo
!         

!         
       return
       end subroutine for_str_sym
