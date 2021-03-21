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

!>     computes k-points for integration
!>     over the brillouin zone. periodic mesh.

       subroutine int_pnt(nbandi,nx,ny,nz,sx,sy,sz,ipr,                  &
     & adot,                                                             &
     & ntrans,mtrx,                                                      &
     & nrk,rk,wgk,nband,indk, kmap,                                      &
     & mxdnrk,mxdbnd)

!      adapted from sverre froyen plane wave program
!      written june 1 1987. jlm
!      modified july 31 1987. jlm
!      version 4.0. 14 october 93. jlm
!      modified 22 march 1999. jlm
!      modified 29 september 2004. jlm  bug found by Richard Hennig and Cyrus Umrigar
!      modified, f90, November 17, 2013. jlm
!      modified, xk, name, kmap, May 26, 2014. JLM
!      modified to be faster O(N^2) -> O(NlogN). 31 July 2014. JLM
!      modified iback. 31 July 2014. JLM
!      modified mxdpnt  11 September 2015.  JLM
!      Modified, documentation, January 2020. JLM

!      copyright  J.L.Martins, INESC-MN.

!      nbandi is the number of eigenvalues to be computed
!      nx, ny, and nz are the number of points in the three
!      directions dermined by the lattice wave vectors. sx, sy, and
!      sz shifts the grid of integration points from the origin.

!      version 4.94

       implicit none

       integer, parameter  :: REAL64 = selected_real_kind(12)

!      input:

       integer, intent(in)             ::  mxdnrk                       !<  size of k-points
       integer, intent(in)             ::  mxdbnd                       !<  size of number of bands

       integer, intent(in)             ::  nbandi                       !<  number of bands for every k-point
       integer, intent(in)             ::  nx,ny,nz                     !<  size of the integration mesh in k-space (nx*ny*nz)
       real(REAL64), intent(in)        ::  sx,sy,sz                     !<  offset of the integration mesh (0.5 for Monkhorst-Pack, 0.0 for DOS)
       integer, intent(in)             ::  ipr                          !<  controls printing (0,1,2). Higher value for more details
       real(REAL64), intent(in)        ::  adot(3,3)                    !<  metric in direct space
       integer, intent(in)             ::  ntrans                       !<  number of symmetry operations in the factor group
       integer, intent(in)             ::  mtrx(3,3,48)                 !<  rotation matrix (in reciprocal lattice coordinates)for the k-th symmetry operation of the factor group

!      output

       integer, intent(out)            ::  nrk                          !<  number of k-points for integration in the irreducible wedge of the brillouin zone
       real(REAL64), intent(out)       ::  rk(3,mxdnrk)                 !<  component in lattice coordinates of the k-point in the mesh
       real(REAL64), intent(out)       ::  wgk(mxdnrk)                  !<  weight in the integration of k-point
       integer, intent(out)            ::  nband(mxdnrk)                !<  number of bands for each k-points
       integer, intent(out)            ::  indk(6,mxdnrk)               !<  index of the six k-points neighbouring k-point i
       integer, intent(out)            ::  kmap(3,nx,ny,nz)             !<  kmap(1,...) corresponding k-point, kmap(2,...) = 1 additional inversion, kmap(3,...) symmetry operation

!      local allocatable arrays

       real(REAL64), allocatable       ::  xk(:,:)
       integer, allocatable            ::  iback(:,:)
       integer, allocatable            ::  kmark(:)
       real(REAL64), allocatable       ::  dkg(:)
       integer, allocatable            ::  indx(:), indxinv(:)

!      local arrays

       real(REAL64)  ::  rktran(3)
       real(REAL64)  ::  bvec(3,3), bdot(3,3), vcell, avec(3,3)



!      local variables

       integer       ::  kdel, ip, jp, kp, im, jm, km
       integer       ::  irk
       real(REAL64)  ::  dx,dy,dz,dw,sumt
       real(REAL64)  ::  diff
       logical       ::  lsame, lsame2
       real(REAL64)  ::  xt(3),xtmod

       integer       ::  mxdpnt

!      counters

       integer   ::  i, j, k, l, n
       integer   ::  i1, i2, i3, im1, im2, im3
       integer   ::  kmax, kmin, kmc

!      constants

       real(REAL64), parameter :: ZERO = 0.0_REAL64 , UM = 1.0_REAL64
       real(REAL64), parameter :: EPS = 0.000001_REAL64


       call adot_to_avec_sym(adot,avec,bvec)
       call adot_to_bdot(adot,vcell,bdot)

       mxdpnt = nx*ny*nz

       allocate(xk(3,mxdpnt))
       allocate(iback(3,mxdpnt))
       allocate(kmark(mxdpnt))
       allocate(dkg(mxdpnt))

!      dimension test

       if (nbandi > mxdbnd) then
         write(6,'("  *** STOPPED in int_pnt.   Change mxdbnd from ",    &
     &        i8,"    to at least  ",i8)') mxdbnd,nbandi

         stop

       endif

!      set up uniform array of k points
!      kmark is used to mark reducible k points and also to
!      map reducible to irreducible k points

       dx = UM / nx
       dy = UM / ny
       dz = UM / nz
       n = 0
       do i=1,nx
       do j=1,ny
       do k=1,nz
         n = n+1
         xk(1,n) = ((i-1) + sx) * dx
         xk(2,n) = ((j-1) + sy) * dy
         xk(3,n) = ((k-1) + sz) * dz
         kmark(n) = n
         iback(1,n) = i
         iback(2,n) = j
         iback(3,n) = k
       enddo
       enddo
       enddo

!      finds distance to closest G-point

       im1 = nint((abs(adot(1,2))+abs(adot(1,3)))/adot(1,1)) + 1
       im2 = nint((abs(adot(2,1))+abs(adot(2,3)))/adot(2,2)) + 1
       im3 = nint((abs(adot(3,1))+abs(adot(3,2)))/adot(3,3)) + 1
       do n = 1,mxdpnt
         dkg(n) = ZERO
         do i=1,3
         do j=1,3
           dkg(n) = dkg(n) + xk(i,n)*bdot(i,j)*xk(j,n)
         enddo
         enddo
         do i1 = -im1,im1
         do i2 = -im2,im2
         do i3 = -im3,im3
           xt(1) = xk(1,n) - i1*UM
           xt(2) = xk(2,n) - i2*UM
           xt(3) = xk(3,n) - i3*UM
           xtmod = ZERO
           do i=1,3
           do j=1,3
             xtmod = xtmod + xt(i)*bdot(i,j)*xt(j)
           enddo
           enddo
           if(xtmod < dkg(n)) dkg(n) = xtmod
         enddo
         enddo
         enddo
       enddo

!      creates a list of increasing lengths

       allocate(indx(mxdpnt))
       allocate(indxinv(mxdpnt))

       call sort(mxdpnt,dkg,indx)

       do i = 1,mxdpnt
         indxinv(indx(i)) = i
       enddo

!      reduce to irreducible zone

       n = mxdpnt
       nrk = 0
       dw = UM / n
       do i=1,n
         if (kmark(i) == i) then

!          new irreducible point
!          mark with negative kmark

           nrk = nrk + 1

           if (nrk > mxdnrk) then
             write(6,'("  *** STOPPED in int_pnt.    Increase mxdnrk",   &
     &           " from ",i8,".  Suggested value: ",i8)')                &
     &                      mxdnrk, min((mxdnrk*n)/i,n/2)

             stop

           endif

           rk(1,nrk) = xk(1,i)
           rk(2,nrk) = xk(2,i)
           rk(3,nrk) = xk(3,i)
           kmark(i) = -nrk
           kmap(1,iback(1,i),iback(2,i),iback(3,i)) = nrk
           kmap(2,iback(1,i),iback(2,i),iback(3,i)) = 0
           kmap(3,iback(1,i),iback(2,i),iback(3,i)) = 1
           wgk(nrk) = dw

           if (i /= n) then

!            first identifies those which have the same shortest distance to a G-vector

             kmin = indxinv(i)
             do k = indxinv(i),1,-1

               if(abs(dkg(i) - dkg(indx(k))) > EPS) exit

               kmin = k
             enddo

             kmax = indxinv(i)
             do k = indxinv(i),n

               if(abs(dkg(i) - dkg(indx(k))) > EPS) exit

               kmax = k
             enddo

!            operate on irreducible rk with the symmetry operations

             do j=1,ntrans

               rktran(1) = mtrx(1,1,j)*rk(1,nrk) +                       &
     &                     mtrx(1,2,j)*rk(2,nrk) +                       &
     &                     mtrx(1,3,j)*rk(3,nrk)
               rktran(2) = mtrx(2,1,j)*rk(1,nrk) +                       &
     &                     mtrx(2,2,j)*rk(2,nrk) +                       &
     &                     mtrx(2,3,j)*rk(3,nrk)
               rktran(3) = mtrx(3,1,j)*rk(1,nrk) +                       &
     &                     mtrx(3,2,j)*rk(2,nrk) +                       &
     &                     mtrx(3,3,j)*rk(3,nrk)

!              translate to interval 0-1.

               do k=1,3
                 kdel = int(rktran(k))
                 if (rktran(k) < 0) kdel = kdel - 1
                 rktran(k) = rktran(k) - kdel
               enddo

!              remove (mark) k points related to irreducible rk by symmetry

               do kmc = kmin,kmax

                 k = indx(kmc)

                 if (kmark(k) == k) then

!                  both the transformed rk ...

                   lsame = .TRUE.
                   do l=1,3
                     diff = rktran(l) - xk(l,k)
                     diff = abs(diff - nint(diff))
                     if (diff > EPS) then
                       lsame = .FALSE.

                       exit

                     endif
                   enddo

                   if(lsame) then
                     wgk(nrk) = wgk(nrk) + dw
                     kmark(k) = nrk
                     kmap(1,iback(1,k),iback(2,k),iback(3,k)) = nrk
                     kmap(2,iback(1,k),iback(2,k),iback(3,k)) = 0
                     kmap(3,iback(1,k),iback(2,k),iback(3,k)) = j
                   else

!                    ... and its inverse (time reversal)

                     lsame2 = .TRUE.
                     do l=1,3
                       diff = rktran(l) + xk(l,k)
                       diff = abs(diff - nint(diff))
                       if (diff > EPS) then
                         lsame2 = .FALSE.

                         exit

                       endif
                     enddo

                     if(lsame2) then
                       wgk(nrk) = wgk(nrk) + dw
                       kmark(k) = nrk
                       kmap(1,iback(1,k),iback(2,k),iback(3,k)) = nrk
                       kmap(2,iback(1,k),iback(2,k),iback(3,k)) = 1
                       kmap(3,iback(1,k),iback(2,k),iback(3,k)) = j
                     endif
                   endif

                 endif

               enddo                                                     ! loop over remaining k-points of same distace

             enddo                                                       ! loop over symmetry operations

           endif

         endif

       enddo                                                             ! loop over nx*ny*nz k-points

!      set up array indk
!      indk(1-6,i) gives the index for irreducible points that
!      are neighbors of point i

       n = 0
       do i=1,nx
       do j=1,ny
       do k=1,nz

         n = n + 1
         if (kmark(n) <= 0) then
           irk = -kmark(n)
           ip = i + 1
           jp = j + 1
           kp = k + 1
           if (i == nx) ip = 1
           if (j == ny) jp = 1
           if (k == nz) kp = 1
           im = i - 1
           jm = j - 1
           km = k - 1
           if (i == 1) im = nx
           if (j == 1) jm = ny
           if (k == 1) km = nz
           ip = ((ip-1)*ny + j -1)*nz + k
           indk(1,irk) = iabs(kmark(ip))
           im = ((im-1)*ny + j -1)*nz + k
           indk(2,irk) = iabs(kmark(im))
           jp = ((i -1)*ny + jp-1)*nz + k
           indk(3,irk) = iabs(kmark(jp))
           jm = ((i -1)*ny + jm-1)*nz + k
           indk(4,irk) = iabs(kmark(jm))
           kp = ((i -1)*ny + j -1)*nz + kp
           indk(5,irk) = iabs(kmark(kp))
           km = ((i -1)*ny + j -1)*nz + km
           indk(6,irk) = iabs(kmark(km))
         endif

       enddo
       enddo
       enddo

       sumt = ZERO
       do i=1,nrk
         nband(i) = nbandi
         sumt = sumt + wgk(i)
       enddo

!      This should never happen, paranoid check...

       if(abs(sumt-UM) > EPS) then
         write(6,'("  *** STOPPED in int_pnt.    Sum of weights",        &
     &             " different from 1")')

         stop

       endif

!      prints the result

       call print_int_pnt(ipr, nbandi, nx,ny,nz, sx,sy,sz,               &
     & adot,                                                             &
     & nrk, rk, wgk, dkg,                                                &
     & mxdnrk)

       deallocate(xk)
       deallocate(iback)
       deallocate(kmark)
       deallocate(dkg)
       deallocate(indx)
       deallocate(indxinv)

       return

       end subroutine int_pnt
