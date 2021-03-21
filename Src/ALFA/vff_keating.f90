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

!>     calculates the keating valence force field energy per unit cell

       subroutine vff_keating(iprint, iowrite,                           &
     & ntype, natom, nameat, rat, adot,                                  &
     & force, energy, stress,                                            &
     & mxdtyp, mxdatm)

!      modified January 2013, J.L.Martins,C.S.Loia
!      Modified, documentation, details, December 2019. JLM           WARNING INCOMPATIBLE WITH EARLIER VERSIONS
!      Modified, vff_check, November 10 2020. JLM
!      copyright inesc-mn. Jose Luis Martins, C.S.Loia

!      version 4.99
!      version 1.6 of MD

       implicit none

       integer, parameter        :: REAL64 = selected_real_kind(12)

!      input variables

       integer, intent(in)         :: iprint                             !<  if iprint == 0 does not print details
       integer, intent(in)         :: iowrite                            !<  tape number

       integer, intent(in)         :: mxdtyp                             !<  maximum number of type of atoms
       integer, intent(in)         :: mxdatm                             !<  maximum number of atoms of a given type
       integer, intent(in)         :: ntype                              !<  number of types of atoms
       integer, intent(in)         :: natom(mxdtyp)                      !<  number of atoms of given type
       character(len=2),intent(in) :: nameat(mxdtyp)                     !<  chemical symbol of type of atom
       real(REAL64), intent(in)    :: rat(3,mxdatm,mxdtyp)               !<  k-th component (in contravariant lattice coordinates) of the position of the n-th atom of type i
       real(REAL64), intent(in)    :: adot(3,3)                          !<  real space metric

!      output variables

       real(REAL64), intent(out)  :: energy                              !<  energy in eV per cell
       real(REAL64), intent(out)  :: force(3,mxdatm,mxdtyp)              !<  k-th component (in contravariant lattice coordinates)  of the keating vff force of the n-th atom of type i
       real(REAL64), intent(out)  :: stress(3,3)                         !<  keating vff stress tensor (in contravariant lattice coordinates)

!      main local variables

       real(REAL64), allocatable, save   :: alfa(:)            ! keating coefficient for bond
       real(REAL64), allocatable, save   :: dist(:)            ! equilibrium distance for bond
       real(REAL64), allocatable, save   :: beta(:,:)          ! keating coefficient for angle

       integer, save                     :: natotal            ! total number of atoms
       real(REAL64), allocatable, save   :: rpr(:,:)           ! atom positions in lattice coordinates
       real(REAL64), allocatable, save   :: frc(:,:)           ! force in lattice coordinates
       integer, allocatable, save        :: ityp(:)            ! type of atom
       integer, save                     :: nbond              ! number of bonds (should be 2*natotal)
       integer, allocatable, save        :: ibond(:,:)         ! atoms at each end of the bond
       integer, allocatable, save        :: itaubd(:,:)        ! lattice translation associated with periodic image for a given atom in bond
       integer, save                     :: nangl              ! number of angles (should be 6*natotal)
       integer, allocatable, save        :: iangl(:,:)         ! atoms that define the angle
       integer, allocatable, save        :: itauan(:,:,:)      ! lattice translation associated with periodic image for a given atom in angle

!      other local variables

       real(REAL64) :: d(3),c(3)
       real(REAL64) :: dd,dij2,aij,cijk,bet

!      counters, etc...

       integer  :: n,k,l
       integer  :: ia,ja,ka,itia,itja,itka,itija,itika,itjka

!      allow execution on first call
  
       integer,save    ::  ifirst = 0

!      parameters

       real(REAL64), parameter  :: ZERO = 0.0_REAL64, UM = 1.0_REAL64

       if(ifirst == 0) then                  
         natotal = 0
         do n=1,ntype
           natotal = natotal + natom(n)
         enddo

         allocate(alfa((mxdtyp*(mxdtyp+1))/2))
         allocate(dist((mxdtyp*(mxdtyp+1))/2))
         allocate(beta(mxdtyp,(mxdtyp*(mxdtyp+1))/2))

         call vff_constants(iprint, iowrite,                             &
     &                 ntype, nameat, alfa, dist, beta, mxdtyp)

         allocate(rpr(3,natotal))
         allocate(frc(3,natotal))
         allocate(ityp(natotal))

         k=0
         do n=1,ntype
           do l=1,natom(n)
             k = k+1
             ityp(k) = n
             rpr(1,k) = rat(1,l,n)
             rpr(2,k) = rat(2,l,n)
             rpr(3,k) = rat(3,l,n)
           enddo
         enddo

         allocate(ibond(2,2*natotal))
         allocate(itaubd(3,2*natotal))
         allocate(iangl(3,6*natotal))
         allocate(itauan(3,2,6*natotal))

         call vff_bond(iprint, iowrite,                                  &
     &                adot, natotal, rpr, nbond, ibond, itaubd,          &
     &                nangl, iangl, itauan, natotal)

         call vff_check(iowrite, nbond, ibond, nangl, iangl,           &
     &                  dist, nameat, ityp, natotal, mxdtyp)

         ifirst = 1

       else
         k=0
         do n=1,ntype
           do l=1,natom(n)
             k = k+1
             rpr(1,k) = rat(1,l,n)
             rpr(2,k) = rat(2,l,n)
             rpr(3,k) = rat(3,l,n)
           enddo
         enddo
  
       endif

!      initializes all the outputs

       energy = ZERO

       do n=1,3
       do k=1,3
         stress(n,k) = ZERO
       enddo
       enddo

       do n=1,natotal
       do k=1,3
         frc(k,n) = ZERO
       enddo
       enddo

!      loop over bonds

       do n=1,nbond

         ia = ibond(1,n)
         ja = ibond(2,n)
         itia = ityp(ia)
         itja = ityp(ja)
         if(itia > itja) then
           itija = (itia*(itia-1))/2 + itja
         else
           itija = (itja*(itja-1))/2 + itia
         endif
         dd = ZERO
         do k=1,3
           d(k)=rpr(k,ja) - rpr(k,ia) + UM*itaubd(k,n)
         enddo
         do k=1,3
         do l=1,3
           dd = dd + d(k)*adot(k,l)*d(l)
         enddo
         enddo
         dij2 = dist(itija)*dist(itija)
         aij = 3*alfa(itija) /( 8*dij2)
         energy = energy + aij*(dd-dij2)*(dd-dij2)
         frc(1,ia) = frc(1,ia) + 2*aij*(dd-dij2)*2*d(1)
         frc(2,ia) = frc(2,ia) + 2*aij*(dd-dij2)*2*d(2)
         frc(3,ia) = frc(3,ia) + 2*aij*(dd-dij2)*2*d(3)
         frc(1,ja) = frc(1,ja) - 2*aij*(dd-dij2)*2*d(1)
         frc(2,ja) = frc(2,ja) - 2*aij*(dd-dij2)*2*d(2)
         frc(3,ja) = frc(3,ja) - 2*aij*(dd-dij2)*2*d(3)
         do k=1,3
         do l=1,3
           stress(k,l) = stress(k,l) - 4*aij*(dd-dij2)*d(k)*d(l)
         enddo
         enddo

       enddo                          !     loop over bonds

!      loop over angles

       do n=1,nangl

         ia = iangl(1,n)
         ja = iangl(2,n)
         ka = iangl(3,n)
         itia = ityp(ia)
         itja = ityp(ja)
         itka = ityp(ka)
         if(itja > itka) then
           itjka = (itja*(itja-1))/2 + itka
         else
           itjka = (itka*(itka-1))/2 + itja
         endif
         if(itia > itja) then
           itija = (itia*(itia-1))/2 + itja
         else
           itija = (itja*(itja-1))/2 + itia
         endif
         if(itia > itka) then
           itika = (itia*(itia-1))/2 + itka
         else
           itika = (itka*(itka-1))/2 + itia
         endif
         dd = ZERO
         do k=1,3
           d(k) = rpr(k,ja) - rpr(k,ia) + UM*itauan(k,1,n)
           c(k) = rpr(k,ka) - rpr(k,ia) + UM*itauan(k,2,n)
         enddo
         do k=1,3
         do l=1,3
           dd = dd + d(k)*adot(k,l)*c(l)
         enddo
         enddo
         cijk = dist(itija)*dist(itika)
         bet = 3*beta(itia,itjka) / (8*cijk)
         cijk =-cijk/(3*UM)
         energy = energy + bet*(dd-cijk)*(dd-cijk)
         frc(1,ia) = frc(1,ia) + 2*bet*(dd-cijk)*(d(1)+c(1))
         frc(2,ia) = frc(2,ia) + 2*bet*(dd-cijk)*(d(2)+c(2))
         frc(3,ia) = frc(3,ia) + 2*bet*(dd-cijk)*(d(3)+c(3))
         frc(1,ja) = frc(1,ja) - 2*bet*(dd-cijk)*c(1)
         frc(2,ja) = frc(2,ja) - 2*bet*(dd-cijk)*c(2)
         frc(3,ja) = frc(3,ja) - 2*bet*(dd-cijk)*c(3)
         frc(1,ka) = frc(1,ka) - 2*bet*(dd-cijk)*d(1)
         frc(2,ka) = frc(2,ka) - 2*bet*(dd-cijk)*d(2)
         frc(3,ka) = frc(3,ka) - 2*bet*(dd-cijk)*d(3)        
         do k=1,3
         do l=1,3
!          keeps stress symmetric
           stress(k,l) = stress(k,l)- 2*bet*(dd-cijk)*                   &
     &                                (d(k)*c(l)+d(l)*c(k))
         enddo
         enddo

       enddo                           !     loop over angles

!      changes storage mode of forces

       k=0
       do n=1,ntype
         do l=1,natom(n)
           k = k+1
           force(1,l,n) = frc(1,k) 
           force(2,l,n) = frc(2,k) 
           force(3,l,n) = frc(3,k) 
         enddo
       enddo

       return

       end subroutine vff_keating
