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

!>     symmetrizes the atomic positions

       subroutine sym_rat(ntype,natom,rat,ntrans,mtrx,tnp,               &
     & mxdtyp,mxdatm)

!      written February 22, 2000. jlm
!      modified March 12, 1999. jlm (time seems non-linear around here...)
!      modified February 6 2012. jlm   symmetrizes in presence of noticeable noise
!      converted to f90, 25 Septenber 2013. jlm
!      Name changed from rat_sym to sym_rat. December 2016. JLM
!      Modified, documentation, December 2019. JLM
!      copyright INESC-MN/Jose Luis Martins

!      version 4.94


       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input:

       integer, intent(in)             ::  mxdtyp                       !<  array dimension of types of atoms
       integer, intent(in)             ::  mxdatm                       !<  array dimension of number of atoms of a given type
       integer, intent(in)             ::  ntype                        !<  number of types of atoms
       integer, intent(in)             ::  natom(mxdtyp)                !<  number of atoms of type i
       integer, intent(in)             ::  ntrans                       !<  number of symmetry operations in the factor group
       integer, intent(in)             ::  mtrx(3,3,48)                 !<  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the factor group
       real(REAL64), intent(in)        ::  tnp(3,48)                    !<  2*pi* i-th component (in lattice coordinates) of the fractional translation vector associated with the k-th symmetry operation of the factor group
       
!      input and output

       real(REAL64), intent(inout)     ::  rat(3,mxdatm,mxdtyp)         !<  lattice coordinates of atom j of type i

!      constants

       real(REAL64), parameter :: PI = 3.14159265358979323846_REAL64
       real(REAL64), parameter :: TWOPI = 2.0_REAL64*PI       
       real(REAL64), parameter :: ZERO = 0.0_REAL64, UM = 1.0_REAL64       

!      scratch

       real(REAL64), allocatable  ::  ratsym(:,:,:)

!      local variables

       integer :: kk,idist,idif(3),icro(3)
       real(REAL64) :: cro(3),xdif,xmaxd

!      counters

       integer ::   i,j,k,l,m


!      initializes ratsym

       allocate(ratsym(3,mxdatm,mxdtyp))

       do i=1,ntype
       do j=1,natom(i)
       do k=1,3
         ratsym(k,j,i) = ZERO
       enddo
       enddo
       enddo

!      finds criterium for same atom ("nearest neighbor distance")

       xmaxd = UM
       do j=1,ntype
         if(natom(j) > 1) then
           do k=1,natom(j)-1
             do l=k+1,natom(j)
               do m=1,3
                 cro(m) = rat(m,k,j)-rat(m,l,j)
                 cro(m) = cro(m) - dble(nint(cro(m)))
                 cro(m) = abs(cro(m))
               enddo
               xmaxd = min(xmaxd,max(cro(1),cro(2),cro(3)))
             enddo
           enddo
         endif
       enddo
       xmaxd = 0.49*xmaxd

!      symmetrizes positions

       do i=1,ntrans
         do j=1,ntype
           do k=1,natom(j)

!                     -1
!            find mtrx    * (rat - tnp)

             do l=1,3
               cro(l) = ZERO
               do m=1,3
                 cro(l) = cro(l)+mtrx(m,l,i)*(rat(m,k,j)-tnp(m,i)/TWOPI)
               enddo
             enddo

             kk = 0
             do l=1,natom(j)
               idist = 0
               do m=1,3
                 xdif = cro(m)-rat(m,l,j)
                 idif(m) = nint(xdif)
                 if (abs(xdif-UM*idif(m)) > xmaxd) idist = idist + 1
               enddo
               if(idist == 0) then
                 kk = l
                 icro(1) = idif(1)
                 icro(2) = idif(2)
                 icro(3) = idif(3)
               endif
             enddo
             if(kk == 0) then
               write(6,'("  *** stopped in rat_sym.  sym operation = ", &
     &         i3," type of atom =",i5,"  atom = ",i5)') i,j,k

               stop

             endif

!            add to position

             do m=1,3
               ratsym(m,kk,j) = ratsym(m,kk,j) + cro(m) - UM*icro(m)
             enddo

           enddo
         enddo
       enddo


!      renormalizes

       do i=1,ntype
       do j=1,natom(i)
       do k=1,3
         rat(k,j,i) = ratsym(k,j,i) / (UM*ntrans)
       enddo
       enddo
       enddo

       deallocate(ratsym)

       return

       end subroutine sym_rat
