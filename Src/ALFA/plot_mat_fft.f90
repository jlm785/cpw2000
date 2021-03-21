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

!>  the matrix for the contour plot is generated
!>  by a super-quadratic interpolation in the fft mesh

subroutine plot_mat_fft(rhopl,nx,ny,chd,id,n1,n2,n3,c0,dx,dy)

! Modified, f90, 27 May 2014. JLM
! Modified, quad_3D32pt, 18 october 2020. JLM
! Documentation, name, 4 february 2021. JLM
! copyright  Jose Luis Martins/INESC-MN

  implicit none

! version 4.99


  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  nx, ny                          !<  Dimensions of grid in plane
  integer, intent(in)                ::  id, n1, n2, n3                  !<  Dimensions of fft grid

  complex(REAL64), intent(in)        ::  chd(id,n2,n3)                   !<  charge density on FFT grid

  real(REAL64), intent(in)           ::  c0(3)                           !<  corner (origin) of plane
  real(REAL64), intent(in)           ::  dx(3),dy(3)                     !<  step vectors that define the plane (lattice coordinates)

! output

  complex(REAL64), intent(out)       ::  rhopl(nx,ny)                    !<  charge density interpolated on the planar grid

! other variables

  real(REAL64)      ::  coor(3)
  real(REAL64)      ::  bl,bu
  integer           ::  nn(3), kl(3)
  real(REAL64)      ::  x(3), fp(-1:2,-1:2,-1:2)
  real(REAL64)      ::  fl(-1:2,-1:2,-1:2)                               !  weight of function on the 4x4x4 grid

! counters

  integer      ::  i, j, k
  integer      ::  i1, i2, i3
  integer      ::  j1, j2, j3


  nn(1) = n1
  nn(2) = n2
  nn(3) = n3
  do i=1,nx
  do j=1,ny
!   This is probably a BUG                              BUG
    do k=1,3
      coor(k) = c0(k) + dx(k)*(i-1) + dy(k)*(j-1)
      coor(k) = nn(k)*coor(k)
      kl(k) = floor(coor(k))
      x(k) = coor(k) - kl(k)
      kl(k) = mod(kl(k),nn(k))
      if(kl(k) < 0) kl(k) = kl(k) + nn(k)
      kl(k) = kl(k) + 1

! loop over 4x4x4 interpolation grid

      do i1 = -1,2
      do i2 = -1,2
      do i3 = -1,2

        j1 = kl(1) + i1
        j2 = kl(2) + i2
        j3 = kl(3) + i3
        if(j1 <= 0) j1 = j1 + nn(1)
        if(j1 > nn(1)) j1 = j1 - nn(1)
        if(j2 <= 0) j2 = j2 + nn(2)
        if(j2 > nn(2)) j2 = j2 - nn(2)
        if(j3 <= 0) j3 = j3 + nn(3)
        if(j3 > nn(3)) j3 = j3 - nn(3)

        fp(i1,i2,i3) = chd(j1,j2,j3)

      enddo
      enddo
      enddo

      call quad_3D32pt(x,fl)

      rhopl(i,j) = 0.0d0

      do i1 = -1,2
      do i2 = -1,2
      do i3 = -1,2
        rhopl(i,j) = rhopl(i,j) + fp(i1,i2,i3) * fl(i1,i2,i3)
      enddo
      enddo
      enddo
    enddo
  enddo
  enddo

  bl = real(rhopl(1,1))
  bu = real(rhopl(1,1))
  do i=1,nx
  do j=1,ny
    if (real(rhopl(i,j)) < bl) bl = real(rhopl(i,j))
    if (real(rhopl(i,j)) > bu) bu = real(rhopl(i,j))
  enddo
  enddo

  write(6,*)
  write(6,'("  min and max of function in plane:  ",2f12.5)') bl,bu
  write(6,*)

  return

end subroutine plot_mat_fft
