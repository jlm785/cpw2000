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

!>     identifies the space group of the crystal

       subroutine sym_space_group(ntrans,mtrx,tnp,tol,                   &
     & nrot,irotdir,irotrec,                                             &
     & ntype,natom,rat,                                                  &
     & mxdtyp,mxdatm)

!      Written 12 November 2003. JLM
!      Modified 19 March 2004.
!      Modified 23 May 2014, tnp noise bug, f90.  JLM
!      Modified, rotrat, 8 June 2014. JLM
!      Modified, documentation, December 2019. JLM
!      Copyright INESC-MN/Jose Luis Martins


!      version 4.94


       implicit none
       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  mxdtyp                     !<  array dimension of types of atoms
       integer, intent(in)                ::  mxdatm                     !<  array dimension of number of atoms of a given type

       real(REAL64), intent(in)           ::  tol                        !<  tolerance for same position

       integer, intent(in)                ::  nrot                       !<  number of operations in lattice point group
       integer, intent(in)                ::  irotdir(3,3,48)            !<  n-th rotation matrix in lattice point group in lattice coordinates x_i = irot(i,j)y_j
       integer, intent(in)                ::  irotrec(3,3,48)            !<  n-th rotation matrix in lattice point group in reciprocal lattice coordinates x_i = irot(i,j)y_j

       integer, intent(in)                ::  ntype                      !<  number of types of atoms
       integer, intent(in)                ::  natom(mxdtyp)              !<  number of atoms of type i
       real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)       !<  k-th component (in lattice coordinates) of the position of the n-th atom of type i

!      output

       integer, intent(out)               ::  ntrans                     !<  number of symmetry operations in the factor group
       integer, intent(out)               ::  mtrx(3,3,48)               !<  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the factor group
       real(REAL64), intent(out)          ::  tnp(3,48)                  !<  2*pi* i-th component (in lattice coordinates) of the fractional translation vector associated with the k-th symmetry operation of the factor group

!      allocatable array

       real(REAL64), allocatable          ::  rotrat(:,:,:)

!      local variables

       real(REAL64)       ::  frac(3), diff(3), sdif
       integer            ::  ib(48)
       integer            ::  itymin, natmin, nfsucc
       integer            ::  ntotal, ncount, jsucc

       logical            ::  lf                         !  an approximation was found
       logical            ::  lsq                        !  it is the square of a rati onal
       integer            ::  nnum, nden                 !  x ~ nnum/ndem
       

!      constants

       real(REAL64), parameter  :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
       real(REAL64), parameter  :: PI = 3.14159265358979323846_REAL64
       real(REAL64), parameter  :: EPS = 0.00001_REAL64                  !  tolerance
       integer, parameter       :: NTRY = 124                            !  tries denominators up to ntry

!      counters

       integer    ::  i, j, k, n, nf, j1, j2


       allocate(rotrat(3,mxdatm,mxdtyp))

!      identifies type of atoms with smallest number

       itymin = 1
       natmin = natom(1)
       do i=1,ntype
         if(natom(i) < natmin) then
           itymin = i
           natmin = natom(i)
         endif
       enddo

!      loop over lattice rotations

       ntrans = 0
       do n=1,nrot

!        rotated atoms

         do i=1,ntype
         do j=1,natom(i)
         do k=1,3
           rotrat(k,j,i) = irotdir(k,1,n)*rat(1,j,i) +                   &
     &                     irotdir(k,2,n)*rat(2,j,i) +                   &
     &                     irotdir(k,3,n)*rat(3,j,i)
         enddo
         enddo
         enddo

!        loop over -fractional translation candidates

         nfsucc = 0
         do nf=1,natom(itymin)

           do k=1,3
             frac(k) = rotrat(k,1,itymin) - rat(k,nf,itymin)
             frac(k) = frac(k) - UM*nint(frac(k))
           enddo

!          tests operation

           ntotal = 0
           ncount = 0
           do i=1,ntype

             ntotal = ntotal + natom(i)
             do j1=1,natom(i)
               jsucc = 0
               do j2=1,natom(i)

                 sdif = ZERO
                 do k=1,3
                   diff(k) = rotrat(k,j2,i) - rat(k,j1,i) - frac(k)
                   diff(k) = diff(k) - UM*nint(diff(k))
                   sdif = sdif + abs(diff(k))
                 enddo
                 if(sdif < tol) then
                   jsucc = 1

                   exit

                 endif
               enddo

               if(jsucc == 1) ncount = ncount + 1

             enddo

           enddo

           if(ncount == ntotal) then
             nfsucc = 1

             exit

           endif

         enddo

         if(nfsucc == 1) then
           ntrans = ntrans+1
           ib(ntrans) = n

!          should be close to a rational number

           do k=1,3

             call near_rational(frac(k),lf,lsq,nnum,nden,NTRY,EPS)

             if(lf .and. .not. lsq) then
               frac(k) = UM*nnum/(UM*nden)
             endif

           enddo

           do k=1,3
             tnp(k,ntrans) = -frac(k)
           enddo
         endif

       enddo

!      stores the result

       do n=1,ntrans
         do j=1,3
         do i=1,3
           mtrx(i,j,n) = irotrec(i,j,ib(n))
         enddo
         enddo
         do j=1,3
           tnp(j,n) = 2*PI*tnp(j,n)
         enddo
       enddo

       deallocate(rotrat)

       return

       end subroutine sym_space_group
