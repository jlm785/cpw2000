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

!>     writes a file for later use in band a structure filtered by atomic
!>     orbital file

       subroutine out_band_info_write(filename,io,                       &
     &      title,subtitle,nstyle,                                       &
     &      neig,nrk,rk,rk_fld,xk,e_of_k,eref,nocc,                      &
     &      nbaslcao,infolcao,basxpsi,pkn,                               &
     &      nvert,xcvert,nlines,ljump,nkstep,label,xklab,ntype,nameat)

!      Adapted by Carlos Loia Reis from out_band_xmgrace.
!      Date unknown.
!      Modified, documentation 29 May 2020. JLM
!      Modiified to write information to QtBandViewer June 2021. CLR
!      copyright  Carlos Loia Reis/Jose Luis Martins/INESC-MN

       implicit none
       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       character(len=*), intent(in)       ::  filename                   !<  file to be written
       integer, intent(in)                ::  io                         !<  tape number 

       character(len=50), intent(in)      ::  title                      !<  title for plots
       character(len=140), intent(in)     ::  subtitle                   !<  subtitle for plots
       integer, intent(in)                ::  nstyle                     !<  choice of plot style

       integer, intent(in)                ::  neig                       !<  number of bands
       integer, intent(in)                ::  nrk                        !<  number of k-vectors
       real(REAL64), intent(in)           ::  xk(nrk)                    !<  x coordinate of k-point in plot
       real(REAL64), intent(in)           ::  e_of_k(neig,nrk)           !<  band energies of k-point in plot

       real(REAL64), intent(in)           ::  eref                       !<  reference energy for plot
       integer, intent(in)                ::  nocc                       !<  number of occupied states (different color)

       integer, intent(in)                ::  nbaslcao                   !<  number of atomic orbitals 
       integer, intent(in)                ::  infolcao(5,nbaslcao)       !<  information about the original atomic orbital.  (type of atom, atom of that type, n,l,m)
       real(REAL64), intent(in)           ::  basxpsi(nbaslcao,neig,nrk) !<  |<bas|psi>|^2 for each k 
       real(REAL64), intent(in)           ::  pkn(nrk,neig)

       integer, intent(in)                ::  nvert                      !<  number of vertical lines in plot
       real(REAL64), intent(in)           ::  xcvert(nvert)              !<  x coordinate of vertical line
       integer, intent(in)                ::  nlines                     !<  number of lines in reciprocal space
       logical, intent(in)                ::  ljump(nlines)              !<  indicates if the new line contains a jump from the preceeding
       integer, intent(in)                ::  nkstep(nlines)             !<  number of steps in line
       character(len=6), intent(in)       ::  label(nvert+nlines)        !<  label of symmetry k-points
       real(REAL64), intent(in)           ::  xklab(nvert+nlines)        !<  x coordinate of label

       integer, intent(in)                ::  ntype                      !<  number of types of atoms
       character(len=2)                   ::  nameat(ntype)              !<  chemical symbol for the type i
       
       real(REAL64), intent(in)           :: rk(3,nrk)
       real(REAL64), intent(in)           :: rk_fld(3,nrk)

!      counters

       integer      ::   n, n10, k, j, jj

!      extra stuff

!      constants
  
       real(REAL64), parameter  :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
       real(REAL64), parameter  ::  EV = 27.21138505_REAL64

       integer      :: i, irk, iband, iorb

       integer, parameter                     :: REAL32 = selected_real_kind(6)

       real(REAL32), allocatable              ::  pkn_out(:,:)
       real(REAL32), allocatable              ::  basxpsi_out(:,:,:)         
       real(REAL32), allocatable              ::  e_of_k_out(:,:)            
       real(REAL32), allocatable              ::  xk_out(:)  
       
       real(REAL64) ymin, ymax, ymtmp                

       real(REAL64), allocatable :: xk_start(:), xk_end(:)
       
       real(REAL32) :: scl
       
!      begin

       allocate(xk_out(nrk))
       allocate(e_of_k_out(neig, nrk))
       allocate(pkn_out(nrk,neig))      
       allocate (basxpsi_out(nbaslcao,neig,nrk))

       allocate (xk_start(nlines))
       allocate (xk_end(nlines))
             
!      finds the energy range for the bands

       ymin = e_of_k(1,1) - eref
       do irk=1,nrk
         ymtmp = e_of_k(neig,irk) - eref
         do n=1,neig
           if(ymin > e_of_k(n,irk)-eref) ymin = e_of_k(n,irk)-eref
           if(ymtmp < e_of_k(n,irk)-eref) ymtmp = e_of_k(n,irk)-eref
         enddo
         if(irk == 1) then
           ymax = ymtmp
         else
           if(ymax > ymtmp) ymax = ymtmp
         endif
       enddo
       ymin = real(nint(ymin*27.212)) - UM
       ymax = real(nint(ymax*27.212)) + UM

!      find xk_start and xk_end
       j=1
       do i=1, nlines
         if (ljump(i)) then
           xk_start(i) = xklab(j)
           xk_end(i)   = xklab(j+1)
           j=j+2
         else
           xk_start(i) = xklab(j-1)
           xk_end(i)   = xklab(j)
           j=j+1
         endif      
       enddo
       
!      transfer large arrays to REAL32

       do irk=1, nrk       
         do iband = 1, neig         
           xk_out(irk) = xk(irk)           
           e_of_k_out(iband, irk) = (e_of_k(iband, irk) -eref)*EV
           pkn_out(irk, iband) = pkn(irk,iband)
           
           if ( pkn_out(irk, iband) > 1.0_REAL32) then
            pkn_out(irk, iband) = 1.0_REAL32
           endif
           
           if ( pkn_out(irk, iband) < 1.0e-8) then
            pkn_out(irk, iband) = 0.0_REAL32
           endif
           do iorb=1, nbaslcao
             basxpsi_out(iorb,iband,irk) = 0.9999*basxpsi(iorb,iband,irk)  ! dont like this a bit , but python does not know how to sum !!!
             if (basxpsi_out(iorb,iband,irk) > 1.0_REAL32 ) then 
              basxpsi_out(iorb,iband,irk) = 1.0_REAL32
              write(*,*) 'out of domain in basxpsi', basxpsi_out(iorb,iband,irk)
             endif
             if (basxpsi_out(iorb,iband,irk) < 1.0e-8 ) then 
              basxpsi_out(iorb,iband,irk) = 0.0
             endif
             
           enddo             ! loop over orbitals
         enddo               ! loop over eigenvalues 
                  
!         stop
         
       enddo                 ! loop over k-points
       
!       basxpsi_out(:,:,:) = 0.0
       
                     

!      begin binary write

       open(unit=504, file=trim(filename)//".bv", form='unformatted')
                     
       write(504) title                     
       write(504) subtitle                  
       write(504) eref                      
       write(504) xcvert(1),xcvert(nvert)                      
       write(504) ymin, ymax                
       write(504) nvert + nlines                   
       write(504) xklab(:)                  
       write(504) label(:)                  
       write(504) nlines  
       write(504) xk_start
       write(504) xk_end                      
       write(504) nrk, neig        
       write(504) xk_out(:)                        
       write(504) e_of_k_out         
       write(504) pkn_out(:,:)        
  
       write(504) nbaslcao
       write(504) infolcao
       write(504) basxpsi_out
       write(504) ntype
       write(504) nameat
       
       write(504) rk
       write(504) rk_fld
       
       close(504)
       
!      end binary write


       open(unit=io,file=filename,form='formatted')

       write(io,'(a50)') title
       write(io,'(a140)') subtitle
       write(io,'(i10,5x,"plot style")') nstyle

       write(io,'(2i10,5x,"number of states and k-points")') neig,nrk

       n10 = neig/10
       do n = 1,nrk
         write(io,'(f14.6,5x,"plot coor. for k-point",i5," followed",    &
     &            " by band energies in Hartree")') xk(n),n
         if(n10 > 0) then
           do j = 1,n10
             write(io,'(5x,10f14.6)') (e_of_k(jj,n),jj=10*(j-1)+1,10*j)
           enddo
         endif
         if(neig > 10*n10) then
           write(io,'(5x,10f14.6)') (e_of_k(jj,n),jj=10*n10+1,neig)
         endif
       enddo

       write(io,'(5x,f14.6,i8,5x,"eref,nocc")') eref,nocc

       write(io,'(2i10,5x,"nvert,nlines")') nvert,nlines

       n10 = nvert/10
       if(n10 > 0) then
         do j = 1,n10
           write(io,'(5x,10f14.6)') (xcvert(jj),jj=10*(j-1)+1,10*j)
         enddo
       endif
       if(nvert > 10*n10) then
         write(io,'(5x,10f14.6)') (xcvert(jj),jj=10*n10+1,nvert)
       endif
 
       n10 = nlines/10
       if(n10 > 0) then
         do j = 1,n10
           write(io,'(5x,10(l5,i5))') (ljump(jj),nkstep(jj),             &
     &                                     jj=10*(j-1)+1,10*j)
         enddo
       endif
       if(nlines > 10*n10) then
         write(io,'(5x,1010(l5,i5))') (ljump(jj),nkstep(jj),             &
     &                                     jj=10*n10+1,nlines)
       endif
  
       n10 = (nvert+nlines)/10
       if(n10 > 0) then
         do j = 1,n10
           write(io,'(5x,10(2x,a6,f14.6))') (label(jj),xklab(jj),        &
     &                                     jj=10*(j-1)+1,10*j)
         enddo
       endif
       if(nvert+nlines > 10*n10) then
         write(io,'(5x,1010(2x,a6,f14.6))') (label(jj),xklab(jj),        &
     &                                     jj=10*n10+1,nvert+nlines)
       endif

       write(io,'(5x,i8,5x,"number of atomic orbitals")') nbaslcao

       do n = 1,nbaslcao
         write(io,'(5x,5i8,5x,"type, atom, n, l ,m")')                   &
     &                                (infolcao(j,n),j=1,5)
       enddo

       n10 = nbaslcao/10
       do k = 1,nrk
         do n = 1,neig
           if(n10 > 0) then
             do j = 1,n10
               write(io,'(5x,10f14.6)') (basxpsi(jj,n,k),                &
     &                                      jj=10*(j-1)+1,10*j)
             enddo
           endif
           if(nbaslcao > 10*n10) then
               write(io,'(5x,10f14.6)') (basxpsi(jj,n,k),                &
     &                                      jj=10*n10+1,nbaslcao)
           endif
         enddo
       enddo
       
       !! extra unfolding info
       do k = 1,nrk
         do n = 1,neig
            write(io,'(10f14.6)') pkn(k,n)
          enddo
       enddo
       !!
       
       !! extra atomic information: ntype natom(ntype)
       
         write(*,*) 'ntype is', ntype
       
         write(io,'(i5)') ntype
         
         do n = 1,ntype
          write(io,'(a2)') nameat(n)
         enddo
             
       !!
       
    
       close(unit=io)

       return
       end subroutine out_band_info_write
