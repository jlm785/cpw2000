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

!>     table of keating vff constants.
!>     that correct LDA results for C, Si, Ge, Sn
!>     Distances in angstroms and Force constants in N/m from table.
!>     On output they are converted to atomic (Hartree/bohr) units.

       subroutine vff_constants(iprint, iowrite,                         &
     &                     ntype, nameat, alfa, dist, beta, mxdtyp)

!      If the force constants are unavailable they are given a HUGE value
!      so it is obvious in the output...

!      written January 2013, J.L.Martins
!      modified May 2014, C.L. Reis,
!      to correct LDA results for C, Si, Ge, Sn
!      the beta constant of Sn is found to be negative.
!      geometric averages are replaced by arithmetic averages
!      when beta is negative.
!      Modified, documentation, details, printing, December 2019. JLM
!      copyright INESC-MN/Jose Luis Martins/C.L. Reis

!      version 4.94

       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input variables

       integer, intent(in)         :: iprint                             !<  if iprint == 0 does not print details
       integer, intent(in)         :: iowrite                            !<  tape number

       integer, intent(in)         :: mxdtyp                             !<  maximum number of type of atoms
       integer, intent(in)         :: ntype                              !<  number of types of atoms
       character(len=2),intent(in) :: nameat(mxdtyp)

!      output variables

       real(REAL64), intent(out)   :: alfa((mxdtyp*(mxdtyp+1))/2)        !<  keating coefficient for bond
       real(REAL64), intent(out)   :: dist((mxdtyp*(mxdtyp+1))/2)        !<  equilibrium distance for bond
       real(REAL64), intent(out)   :: beta(mxdtyp,(mxdtyp*(mxdtyp+1))/2) !< keating coefficient for angle

!      counters

       integer :: n,m,k,nn2,mm2,nm2,mk2,kk2

!      physical constants

       real(REAL64), parameter  :: BOHR = 0.52917721E-10
       real(REAL64), parameter  :: EV = 27.211385
       real(REAL64), parameter  :: HARTREE = EV * 1.6021765E-19
       real(REAL64), parameter  :: ANG = 1E-10


!      same species alpha,dist,beta (group IV)

       do n=1,ntype
         nn2 = (n*(n+1))/2
         if(nameat(n) == 'C ' .or. nameat(n) == ' C') then
           dist(nn2) = 1.54455630764955*1.1563

           alfa(nn2)   = 13.6248740466374851
           beta(n,nn2) = 6.1553571220124841

         elseif(nameat(n) == 'Si') then        
           dist(nn2) = 2.35169198397664*1.07506

           alfa(nn2)   = 6.3503462524617307
           beta(n,nn2) = 1.2601737048539403

         elseif(nameat(n) == 'Ge') then          
           dist(nn2)  = 2.44994256603599*0.97698

           alfa(nn2)  = 2.9182378149662425
           beta(n,nn2)= 0.4584868373542476

         elseif(nameat(n) == 'Sn') then
           dist(nn2) = 1.00786*2.80722134636726

           alfa(nn2)   = 14.7836507262404933
           beta(n,nn2) = -3.8734705709774992     

         elseif(nameat(n) == 'Pb') then
           dist(nn2) = 2.99_REAL64

           alfa(nn2) = 0.0_REAL64
           beta(n,nn2) = 0.25_REAL64*alfa(nn2)

         else
           dist(nn2) = -1.0

           alfa(nn2) = 1000000.0
           beta(n,nn2) = alfa(nn2)
         endif
       enddo

       if(ntype > 1) then

!        pairs of species alpha,dist,beta

         do n=1,ntype-1
         do m=n+1,ntype
           nm2 = (m*(m-1))/2+n
           nn2 = (n*(n+1))/2
           mm2 = (m*(m+1))/2

!          group IV

           if( (nameat(n) == 'C ' .or. nameat(n) == ' C' .or.             &
     &          nameat(n) == 'Si' .or. nameat(n) == 'Ge' .or.             &
     &          nameat(n) == 'Sn' .or. nameat(n) == 'Pb') .and.           &
     &         (nameat(m) == 'C ' .or. nameat(m) == ' C' .or.             &
     &          nameat(m) == 'Si' .or. nameat(m) == 'Ge' .or.             &
     &          nameat(m) == 'Sn' .or. nameat(m) == 'Pb')  ) then

             alfa(nm2)   = sqrt(alfa(nn2)*alfa(mm2))
             beta(n,mm2) = sqrt(beta(n,nn2)*beta(m,mm2))
             dist(nm2)   = sqrt(dist(nn2)*dist(mm2))

             if(((nameat(n) == 'C ' .or. nameat(n) == ' C') .and.         &
     &            nameat(m) == 'Sn') .or. (nameat(n) == 'Sn' .and.        &
     &            (nameat(m) == 'C ' .or. nameat(m) == ' C'))) then

               beta(n,mm2) = (beta(n,nn2)+beta(m,mm2))/2

             elseif((nameat(n) == 'Si' .and. nameat(m) == 'Sn') .or.      &
     &            (nameat(n) == 'Sn' .and. nameat(m) == 'Si')) then

               beta(n,mm2) = (beta(n,nn2)+beta(m,mm2))/2

             elseif((nameat(n) == 'Ge' .and. nameat(m) == 'Sn') .or.      &
     &            (nameat(n) == 'Sn' .and. nameat(m) == 'Ge')) then

               beta(n,mm2) = (beta(n,nn2)+beta(m,mm2))/2

             endif

           else
             alfa(nm2) = 1000000.0
             dist(nm2) = -1.0
             beta(n,mm2) = 1.0*alfa(nm2)

             stop 'Keating Correction not available for these elements'

           endif

!          beta(m,nn2) = beta(n,mm2)
!          beta(n,nm2) = sqrt(beta(n,mm2)*beta(n,nn2))
!          beta(m,nm2) = sqrt(beta(m,nn2)*beta(m,mm2))

           beta(m,nn2) = beta(n,mm2)
           if((beta(n,mm2)<0.0D0) .or.(beta(n,nn2)<0.0D0)) then
             beta(n,nm2) = (beta(n,mm2)+beta(n,nn2))/2
           else
             beta(n,nm2) = sqrt(beta(n,mm2)*beta(n,nn2))
          endif

           if((beta(m,nn2)<0.0D0) .or. (beta(m,mm2)<0.0D0)) then
             beta(m,nm2) =  (beta(m,nn2)+beta(m,mm2))/2
           else
             beta(m,nm2) = sqrt(beta(m,nn2)*beta(m,mm2))
           endif


         enddo
         enddo

         if(ntype > 2) then

!          triplets of species beta

           do n=1,ntype
             do m=1,ntype-1
               do k=m+1,ntype
                 if(m /= n .and. k /= n) then
                   mk2 = (k*(k-1))/2 + m
                   mm2 = (m*(m+1))/2
                   kk2 = (k*(k+1))/2
                   if((beta(n,mm2)<0.0D0) .or. (beta(n,kk2)<0.0D0)) then
                     beta(n,mk2) = (beta(n,mm2)+beta(n,kk2))/2                  
                   else 
                     beta(n,mk2) = sqrt(beta(n,mm2)*beta(n,kk2))                  
                   endif
                 endif
               enddo
             enddo
           enddo

         endif

       endif

       if(iprint /= 0) then

         write(iowrite,*)
         write(iowrite,*)
         write(iowrite,*) '    alpha  constants (N/m) '
         write(iowrite,*)
         write(iowrite,'(6x,10(10x,a2))') (nameat(m),m=1,ntype)
         do n=1,ntype
           write(iowrite,'(2x,a2,4x,10f12.2)') nameat(n),                &
     &                     (alfa((n*(n-1))/2+k),k=1,n)
         enddo

         write(iowrite,*)
         write(iowrite,*) '    distances (Angstrom)'
         write(iowrite,*)
         write(iowrite,'(6x,10(10x,a2))') (nameat(m),m=1,ntype)
         do n=1,ntype
           write(iowrite,'(2x,a2,4x,10f12.3)') nameat(n),                &
     &                     (dist((n*(n-1))/2+k),k=1,n)
         enddo     

         write(iowrite,*)
         write(iowrite,*) '    beta constants (N/m)'
         write(iowrite,*)
         do n=1,ntype
           write(iowrite,*)
           write(iowrite,'(" corner atom:  ",a2)') nameat(n)
           write(iowrite,*)
           write(iowrite,'(6x,10(10x,a2))') (nameat(m),m=1,ntype)
           do m=1,ntype
           write(iowrite,'(2x,a2,4x,10f12.2)') nameat(m),                &
     &                     (beta(n,(m*(m-1))/2+k),k=1,m)
           enddo
         enddo

       endif

!      convert to atomic units

       do n=1,(ntype*(ntype+1))/2
         alfa(n) = alfa(n) * BOHR*BOHR / HARTREE
         dist(n) = dist(n) * ANG / BOHR
         do m=1,ntype
           beta(m,n) = beta(m,n) * BOHR*BOHR / HARTREE
         enddo
       enddo

       return

       end subroutine vff_constants

