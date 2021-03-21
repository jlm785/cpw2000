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

!>     Expands the atomic orbitals in the plane-wave basis
!>     Orbitals may have the inversion symmetry applied to them.

       subroutine atomic_orbital_c16(rkpt,mtxd,isort,icmplx,             &
     & nbaslcao,baslcao,infolcao,                                        &
     & ng,kgv,                                                           &
     & norbat,nqwf,delqwf,wvfao,lorb,                                    &
     & ntype,natom,rat,adot,                                             &
     & mxdtyp,mxdatm,mxdlqp,mxddim,mxdorb,mxdgve,mxdlao)


!      written f90/complex  19 June 2012
!      Modified 7 January 2014, style. JLM
!      modified, dimensions vkb, real rylm, March 31, 2014. jlm
!      written January 2008. JLM
!      Modified April 18, 2014. JLM
!      Modified April 21, 2014, icmplx. JLM
!      Modified (kgv in fi), 2014. CLR
!      Final adaptation from two sources, November 2015. JLM
!      Added infolacao, 30 May 2019.  JLM
!      Modified documentation, January 2020. JLM
!      Merge with version 2 (infolcao). JLM

!      copyright INESC-MN/Jose Luis Martins, Carlos Loia Reis

!      version 4.98

       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  mxdtyp                     !<  array dimension of types of atoms
       integer, intent(in)                ::  mxdatm                     !<  array dimension of number of atoms of a given type
       integer, intent(in)                ::  mxddim                     !<  array dimension of plane-waves
       integer, intent(in)                ::  mxdlqp                     !<  array dimension for local potential
       integer, intent(in)                ::  mxdorb                     !<  array dimension of number of local orbitals
       integer, intent(in)                ::  mxdgve                     !<  array dimension of G-space vectors
       integer, intent(in)                ::  mxdlao                     !<  array dimension of orbital per atom type

       real(REAL64), intent(in)           ::  rkpt(3)                    !<  k-point reciprocal lattice coordinates
       integer, intent(in)                ::  mtxd                       !<  wavefunction dimension
       integer, intent(in)                ::  isort(mxddim)              !<  G-vector corresponding to coefficient i of wavefunction

       integer, intent(in)                ::  icmplx                     !<  indicates if the structure factor is complex

       integer, intent(in)                ::  ng                         !<  size of g-space
       integer, intent(in)                ::  kgv(3,mxdgve)              !<  G-vectors in reciprocal lattice coordinates

       integer, intent(in)                ::  norbat(mxdtyp)             !<  number of atomic orbitals for atom k
       integer, intent(in)                ::  nqwf(mxdtyp)               !<  number of points for wavefunction interpolation for atom k
       real(REAL64), intent(in)           ::  delqwf(mxdtyp)             !<  step used in the wavefunction interpolation for atom k
       real(REAL64), intent(in)      ::  wvfao(-2:mxdlqp,mxdlao,mxdtyp)  !<  wavefunction for atom k, ang. mom. l (normalized to vcell)
       integer, intent(in)                ::  lorb(mxdlao,mxdtyp)        !<  angular momentum of orbital n of atom k

       integer, intent(in)                ::  ntype                      !<  number of types of atoms
       integer, intent(in)                ::  natom(mxdtyp)              !<  number of atoms of type i
       real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)       !<  lattice coordinates of atom j of type i
       real(REAL64), intent(in)           ::  adot(3,3)                  !<  metric in real space

!      output

       integer, intent(out)               ::  nbaslcao                   !<  number of atomic orbitals
       complex(REAL64), intent(out)       ::  baslcao(mxddim,mxdorb)     !<  atomic orbitals in plane-wave basis
       integer, intent(out)               ::  infolcao(5,mxdorb)         !<  information about the original atomic orbital.  (type of atom, atom of that type, n,l,m)

!      local allocatable variables

       complex(REAL64), allocatable  ::  st(:)
       integer, allocatable  ::  irow(:)
       integer, allocatable  ::  nq(:)

!      local variables

       real(REAL64)    ::  avec(3,3),bvec(3,3)
       real(REAL64)    ::  vcell, bdot(3,3),volfac
       real(REAL64)    ::  qi,qk(3),qcar(3),fi,xni,xi,vq0

       real(REAL64)    ::  rylm(0:3,-3:3)                 !   r**l sqrt(4 pi / 2l+1) "Y_lm(l,m)+-Y_lm(l,-m)"
       real(REAL64)    ::  drylmdqc(3,0:3,-3:3)           !   unused
       real(REAL64)    ::  d2rylmdqc2(3,3,0:3,-3:3)       !   unused

       complex(REAL64) ::  zylm(0:3,-3:3)                 !   unused
       complex(REAL64) ::  dzylmdqc(3,0:3,-3:3)           !   unused
       complex(REAL64) ::  d2zylmdqc2(3,3,0:3,-3:3)       !   unused

       complex(REAL64) ::  pitol                          !   i**l

       integer         ::  lmax, ifail, icount
       real(REAL64)    ::  diff, dist
       complex(REAL64) ::  olda, oldb

!      constants

       real(REAL64), parameter :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
       complex(REAL64), parameter :: C_I = cmplx(ZERO,UM,REAL64)
       real(REAL64), parameter :: PI = 3.14159265358979323846_REAL64
       real(REAL64), parameter :: EPS = 0.00000001_REAL64

!      counters

       integer         ::  ind,i,j,k,l,m,kk,ni
       integer         ::  nt, n


!      paranoid check, stops complaint about unused ng...

       if(ng > mxdgve) stop

       allocate(st(mxdatm))

!      cell stuff

       call adot_to_avec(adot,avec,bvec)

       call adot_to_bdot(adot,vcell,bdot)

       volfac = UM / sqrt(vcell)

!      checks dimensions and maximum angular momentum

       ind = 0
       lmax = 0
       do k = 1,ntype
       do j = 1,norbat(k)
         ind = ind + (2*lorb(j,k)+1)*natom(k)
         lmax = max(lmax,lorb(j,k))
       enddo
       enddo
       nbaslcao = ind

       if(nbaslcao > mxdorb) then
         write(6,'("  atomic_orbital_c16:    increase mxdorb from ",i8,  &
     &            " to ",i8)')  mxdorb,nbaslcao

         stop

       endif

       if(lmax > mxdlao) then
         write(6,'("  atomic_orbital_c16:    increase mxdlao from ",i8,  &
     &            " to ",i8)')  mxdlao,lmax

         stop

       endif

!      starts loop over g-vectors

       do i=1,mtxd
         qk(1) = rkpt(1) + UM*(kgv(1,isort(i)))
         qk(2) = rkpt(2) + UM*(kgv(2,isort(i)))
         qk(3) = rkpt(3) + UM*(kgv(3,isort(i)))
         qcar(1) = bvec(1,1)*qk(1)+bvec(1,2)*qk(2)+bvec(1,3)*qk(3)
         qcar(2) = bvec(2,1)*qk(1)+bvec(2,2)*qk(2)+bvec(2,3)*qk(3)
         qcar(3) = bvec(3,1)*qk(1)+bvec(3,2)*qk(2)+bvec(3,3)*qk(3)
         qi = qcar(1)*qcar(1) + qcar(2)*qcar(2) + qcar(3)*qcar(3)
         qi = sqrt(qi)

!        angular functions

         call ylm_rc(qcar,rylm,drylmdqc,d2rylmdqc2,.FALSE.,              &
     &            zylm,dzylmdqc,d2zylmdqc2,lmax,0,ifail)

!        starts loop over second index

         ind = 0

         do k=1,ntype

!          structure phase factor

           do kk=1,natom(k)
             fi = qk(1)*rat(1,kk,k) + qk(2)*rat(2,kk,k) +                &
     &            qk(3)*rat(3,kk,k)
             fi = -2*PI*fi
             st(kk) = cmplx(cos(fi),sin(fi),REAL64)
           enddo

           do j = 1,norbat(k)
             l = lorb(j,k)

             xni = qi/delqwf(k)
             ni = nint(xni)
             if(ni < 0) ni = 0

             vq0 = ZERO
             if(ni < nqwf(k)-1) then
               xi  = xni - UM*ni
               vq0 = wvfao(ni,j,k) * (UM+xi) * (UM-xi)                   &
     &              + (UM/2) * (wvfao(ni+1,j,k)*(UM+xi)                  &
     &                         -wvfao(ni-1,j,k)*(UM-xi)) * xi
               vq0 = volfac * vq0
             endif

!            adds phase so that orbital is real for rkpt = 0

             pitol = C_I**l

             do m=-l,l
               do kk=1,natom(k)
                 ind = ind + 1

                 baslcao(i,ind) = pitol*st(kk)*rylm(l,m)*vq0

               enddo
             enddo


           enddo    !  j=1,norbat(k)

         enddo      !  k=1,ntype

       enddo        !  i=1,mtxd

       deallocate(st)

!      performs changes of basis for special case

       if(icmplx == 0) then

!        constructs symetrical and antisymetrical functions
!        for the case of real hamiltonians

         allocate(irow(mxdatm))

         icount = 0
         do nt=1,ntype
           do i=1,natom(nt)
             irow(i) = 0
           enddo
           do i=1,natom(nt)
             do j=i,natom(nt)
               dist = ZERO
               do k=1,3
                 diff = rat(k,j,nt)+rat(k,i,nt)
                 diff = diff - nint(diff)
                 dist = dist + diff*diff
               enddo
               if(dist < EPS) then
                 irow(i) = j
                 irow(j) = i
               endif
             enddo
             if(irow(i) == 0) then
               write(6,'("  Stopped in atomic_orbital_c16:   Failed ",   &
     &                   "to find inversion")')

               stop

             endif
           enddo

           do n=1,norbat(nt)
             l = lorb(n,nt)
             do m=-l,l
               do i=1,natom(nt)
                 if(irow(i) > i) then
                   if(mod(l,2) == 0) then
                     do j=1,mtxd
                       olda = baslcao(j,icount+i)
                       oldb = baslcao(j,icount+irow(i))
                       baslcao(j,icount+i) = (olda+oldb)/sqrt(2*UM)
                       baslcao(j,icount+irow(i)) = C_I *                 &
     &                                       (olda-oldb)/sqrt(2*UM)
                    enddo
                   else
                     do j=1,mtxd
                       olda = baslcao(j,icount+i)
                       oldb = baslcao(j,icount+irow(i))
                       baslcao(j,icount+i) = (olda-oldb)/sqrt(2*UM)
                       baslcao(j,icount+irow(i)) = C_I *                 &
     &                                       (olda+oldb)/sqrt(2*UM)
                     enddo
                   endif
                 elseif(irow(i) == i) then
                   if(mod(l,2) == 1) then
                     do j=1,mtxd
                       baslcao(j,icount+i) = C_I * baslcao(j,icount+i)
                     enddo
                   endif
                 endif
               enddo
               icount = icount + natom(nt)
             enddo
           enddo
         enddo

         deallocate(irow)

       endif

!      fills infolcao, loops must be consistent with previous loops!

       lmax = 0
       do k = 1,ntype
       do j = 1,norbat(k)
         if(lmax < lorb(j,k)) lmax = lorb(j,k)
       enddo
       enddo

       allocate(nq(0:lmax))

       ind = 0

       do k = 1,ntype

         do l = 0,lmax
           nq(l) = 0
         enddo

         do j = 1,norbat(k)

           l = lorb(j,k)
           nq(l) = nq(l) + 1

           do m = -l,l
             do kk = 1,natom(k)
               ind = ind + 1

               infolcao(1,ind) = k
               infolcao(2,ind) = kk
               infolcao(3,ind) = nq(l)
               infolcao(4,ind) = l
               infolcao(5,ind) = m

             enddo
           enddo

         enddo    !  j=1,norbat(k)

       enddo      !  k=1,ntype

       deallocate(nq)

       return

       end subroutine atomic_orbital_c16
