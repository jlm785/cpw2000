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

!>     calculates the bonds and bond angles in a tetrahedraly bonded
!>     crystal.

       subroutine vff_bond(iprint, iowrite,                               &
     &                     adot, natotal, rpr, nbond, ibond, itaubd,      &
     &                     nangl, iangl, itauan, mxat)

!      Modified January 2013. J.L.Martins, C.S.Loia.
!      Modified, documentation, details, printing, December 2019. JLM
!      Modified, bug of translation symmetry, 28 June 2020. JLM
!      copyright   INESC-MN, J.L.Martins, C.S.Loia.

!      version 4.982

       implicit none

       integer, parameter        :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)       :: iprint                               !<  if iprint == 0 does not print details
       integer, intent(in)       :: iowrite                              !<  tape number

       integer, intent(in)       :: mxat                                 !<  maximum number of atoms
       real(REAL64), intent(in)  :: adot(3,3)                            !<  real space metric (a_i . a_j)
       integer, intent(in)       :: natotal                              !<  total number of atoms
       real(REAL64), intent(in)  :: rpr(3,mxat)                          !<  atom positions in lattice coordinates

!      output

       integer, intent(out)      :: nbond                                !<  number of bonds (should be 2*natotal)
       integer, intent(out)      :: ibond(2,2*mxat)                      !<  atoms at each end of the bond
       integer, intent(out)      :: itaubd(3,2*mxat)                     !<  lattice translation associated with periodic image for a given atom in bond
       integer, intent(out)      :: nangl                                !<  number of angles (should be 6*natotal)
       integer, intent(out)      :: iangl(3,6*mxat)                      !<  atoms that define the angle
       integer, intent(out)      :: itauan(3,2,6*mxat)                   !<  lattice translation associated with periodic image for a given atom in angle

!      internal variables

       real(REAL64)   :: adist(5)
       integer        :: jd(5),ntau(3,5)
       integer        :: jqd(4,mxat),nqtau(3,4,mxat),iflag(4,mxat)
       real(REAL64)   :: d(3)
       real(REAL64)   :: xsum

!      counters, etc...

       integer  :: i,j,k,l,m,n
       integer  :: mbi,mbj,nct,mbjs,js
       integer  :: ntetra,nc,jc,kc,n1,n2,n3

!      find four nearest neighbours

!      loop over atoms in cell

       do i=1,natotal

         do m=1,4
           adist(m) = 1000000.0*adot(1,1)
         enddo

!        loop over possible neighbours (atoms in cell and neighboring cells)

         do j=1,natotal
        
           do n1=-2,2
           do n2=-2,2
           do n3=-2,2
             if(n1 /=0 .or. n2 /=0 .or. n3/=0 .or. i /= j) then
               xsum = 0.0
               ntau(1,5) = n1 - nint(rpr(1,j)) + nint(rpr(1,i))
               ntau(2,5) = n2 - nint(rpr(2,j)) + nint(rpr(2,i))
               ntau(3,5) = n3 - nint(rpr(3,j)) + nint(rpr(3,i))
               d(1) = rpr(1,j) - rpr(1,i) + ntau(1,5)
               d(2) = rpr(2,j) - rpr(2,i) + ntau(2,5)
               d(3) = rpr(3,j) - rpr(3,i) + ntau(3,5)
               do k=1,3
               do l=1,3
                 xsum = xsum + d(k)*adot(k,l)*d(l)
               enddo
               enddo
               adist(5) = xsum
               jd(5)=j
              
               call vff_order(adist,jd,ntau)

             endif
           enddo
           enddo
           enddo

         enddo                !       loop over possible neighbours
        
         do m=1,4
           jqd(m,i) = jd(m)
           do l=1,3
             nqtau(l,m,i) = ntau(l,m)
           enddo
         enddo
        
       enddo                  ! loop over atoms in cell

!      check if topology is tetrahedral and select bonds

       ntetra=0
       nc=0
       do i=1,natotal
       do m=1,4
         iflag(m,i)=0
       enddo
       enddo

!      loop over atoms in cell

       do i=1,natotal

!        loop over bonds

         do mbi=1,4

           if(iflag(mbi,i) == 0) then

             j = jqd(mbi,i)
             nct = 0
             do mbj=1,4
               if(iflag(mbj,j) == 0) then
                 if(jqd(mbj,j) == i) then
                   n1 = nqtau(1,mbi,i)+nqtau(1,mbj,j)
                   n2 = nqtau(2,mbi,i)+nqtau(2,mbj,j)
                   n3 = nqtau(3,mbi,i)+nqtau(3,mbj,j)
                   if(n1 == 0 .and. n2 == 0 .and. n3 == 0) then
                     nct = nct+1
                     mbjs = mbj
                     js = j
                   endif
                 endif
               endif
             enddo
             
             if(nct == 1) then
               iflag(mbi,i) = js
               iflag(mbjs,js) = i
               nc = nc+1
               ibond(1,nc) = i
               ibond(2,nc) = js
               do k=1,3
                 itaubd(k,nc) = nqtau(k,mbi,i)
               enddo
             else
               ntetra = ntetra+1
             endif

           endif

         enddo         !       loop over bonds
        
       enddo           !       loop over atoms in cell
      
       nbond=nc

!      checks if topology is really tetrahedral

       if(ntetra /= 0 .or. nbond /= 2*natotal) then
      
         write(iowrite,*)
         write(iowrite,*)  ' The topology is not tetrahedral   '
         write(iowrite,*)
        
         write(iowrite,'(3i8)')  ntetra,nbond,natotal
         do i=1,natotal
         do m=1,4
           if(iflag(m,i) == 0) then
             write(iowrite,'(3i8)')  i,m,jqd(m,i)
             write(iowrite,'(3i8)')  (nqtau(l,m,i),l=1,3)
           endif
         enddo
         enddo

         stop

       endif

!      fills bond angle information

       nc=0
       do i=1,natotal
         do jc=1,3
           do kc=jc+1,4
             nc = nc+1
             iangl(1,nc) = i
             iangl(2,nc) = jqd(jc,i)
             iangl(3,nc) = jqd(kc,i)
             do k=1,3
               itauan(k,1,nc) = nqtau(k,jc,i)
               itauan(k,2,nc) = nqtau(k,kc,i)
             enddo
           enddo
         enddo
       enddo
       nangl=nc
      
!      prints bond and angle information

       if(iprint /= 0) then
      
         write(iowrite,*)
         do n=1,nbond
           write(iowrite,'("  bond ",i3,"   between atoms",2i3,          &
     &             "    with tau ",3i3)') n,(ibond(k,n),k=1,2),          &
     &             (itaubd(l,n),l=1,3)
         enddo
         write(iowrite,*)
         write(iowrite,*)

         do n=1,nangl
           write(iowrite,'("  angle ",i3,"   between atoms ",i3,i5,i3,    &
     &         "    with tau",3i3,2x,3i3)') n,(iangl(k,n),k=1,3),         &
     &              ((itauan(l,k,n),l=1,3),k=1,2)
         enddo
         write(iowrite,*)
         write(iowrite,*)

       endif

       return

       end subroutine vff_bond
