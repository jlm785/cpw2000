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

!>  complementary error function erfc=1-erf
!>
!>  \author       Jose Luis Martins
!>  \version      4.94
!>  \date         Before December 1984.  6 June 2024.
!>  \copyright    GNU Public License v2

function erfc(x)

! July 1977 edition.  W. Fullerton, C3, Los Alamos Sci. Lab.
! Modified by J.L.Martins Dec 84
! 16 significant figures
! Modified documentation August 2019.  JLM
! Converted to f90, December 2016.  JLM
! Indentation, 6 June 2024. JLM


  implicit none

  integer, parameter  :: REAL64 = selected_real_kind(12)

! input

  real(REAL64), intent(in)   ::  x                                       !<  argument

! result

  real(REAL64)               ::  erfc                                    !<  erfc(x)

! local variables

  real(REAL64)       ::  b0, b1, b2
  real(REAL64)       ::  y, yt, z, twoz
  real(REAL64)       ::  xmax, xsml, sqeps

  real(REAL64), parameter    ::  SPI = 1.772453850905516_REAL64


  real(REAL64), parameter   ::  erfcs(13) =                              &
      (/   -.049046121234691808_REAL64 ,                                 &
           -.142261205103713640_REAL64 ,                                 &
            .010035582187599796_REAL64 ,                                 &
           -.000576876469976748_REAL64 ,                                 &
            .000027419931252196_REAL64 ,                                 &
           -.000001104317550734_REAL64 ,                                 &
            .000000038488755420_REAL64 ,                                 &
           -.000000001180858253_REAL64 ,                                 &
            .000000000032334215_REAL64 ,                                 &
           -.000000000000799101_REAL64 ,                                 &
            .000000000000017990_REAL64 ,                                 &
           -.000000000000000371_REAL64 ,                                 &
            .000000000000000007_REAL64   /)

  real(REAL64), parameter   ::  erc2cs(23) =                             &
      (/   -.069601346602309501_REAL64 ,                                 &
           -.041101339362620893_REAL64 ,                                 &
            .003914495866689626_REAL64 ,                                 &
           -.000490639565054897_REAL64 ,                                 &
            .000071574790013770_REAL64 ,                                 &
           -.000011530716341312_REAL64 ,                                 &
            .000001994670590201_REAL64 ,                                 &
           -.000000364266647159_REAL64 ,                                 &
            .000000069443726100_REAL64 ,                                 &
           -.000000013712209021_REAL64 ,                                 &
            .000000002788389661_REAL64 ,                                 &
           -.000000000581416472_REAL64 ,                                 &
            .000000000123892049_REAL64 ,                                 &
           -.000000000026906391_REAL64 ,                                 &
            .000000000005942614_REAL64 ,                                 &
           -.000000000001332386_REAL64 ,                                 &
            .000000000000302804_REAL64 ,                                 &
           -.000000000000069666_REAL64 ,                                 &
            .000000000000016208_REAL64 ,                                 &
           -.000000000000003809_REAL64 ,                                 &
            .000000000000000904_REAL64 ,                                 &
           -.000000000000000216_REAL64 ,                                 &
            .000000000000000052_REAL64 /)

  real(REAL64), parameter   ::  erfccs(24) =                             &
     (/   0.071517931020292500_REAL64 ,                                  &
          -.026532434337606719_REAL64 ,                                  &
           .001711153977920853_REAL64 ,                                  &
          -.000163751663458512_REAL64 ,                                  &
           .000019871293500549_REAL64 ,                                  &
          -.000002843712412769_REAL64 ,                                  &
           .000000460616130901_REAL64 ,                                  &
          -.000000082277530261_REAL64 ,                                  &
           .000000015921418724_REAL64 ,                                  &
          -.000000003295071356_REAL64 ,                                  &
           .000000000722343973_REAL64 ,                                  &
          -.000000000166485584_REAL64 ,                                  &
           .000000000040103931_REAL64 ,                                  &
          -.000000000010048164_REAL64 ,                                  &
           .000000000002608272_REAL64 ,                                  &
          -.000000000000699105_REAL64 ,                                  &
           .000000000000192946_REAL64 ,                                  &
          -.000000000000054704_REAL64 ,                                  &
           .000000000000015901_REAL64 ,                                  &
          -.000000000000004729_REAL64 ,                                  &
           .000000000000001432_REAL64 ,                                  &
          -.000000000000000439_REAL64 ,                                  &
           .000000000000000138_REAL64 ,                                  &
          -.000000000000000048_REAL64 /)

  real(REAL64), parameter  :: ZERO = 0.0_REAL64, UM = 1.0_REAL64

! counters

  integer     ::  i

! smallest positive number on computer
! small=6.8d-8618
! relative precision
! rel=3.6d-15
! xsml=-sqrt(-dlog(spi*rel))
! xmax=sqrt(-dlog(spi*small))
! xmax=xmax-0.5*dlog(xmax)/xmax-0.01
! sqeps=sqrt(2.0*rel)

  xsml=-5.7
  xmax=140.0
  sqeps=8.5d-8

! assymptotic values

  if(x < xsml) then

    erfc = 2*UM

  elseif(x > xmax) then

    erfc = ZERO

  else

    y=abs(x)

    if(y < sqeps) then

      erfc = UM - 2*x/SPI

    elseif(y < UM) then
      z = 2*x*x - UM
      twoz = 2*z
      b1 = ZERO
      b0 = ZERO

      do i = 1,13
        b2 = b1
        b1 = b0
        b0 = twoz*b1 - b2 + erfcs(14-i)
      enddo

      erfc=UM - x*(UM + (b0-b2)/2)

    else

      yt = y*y

      if(yt <= 4.0) then

        z = (8*UM/yt - 5*UM) / (3*UM)
        b1 = ZERO
        b0 = ZERO
        twoz = 2*z

        do i = 1,23
          b2 = b1
          b1 = b0
          b0 = twoz*b1 - b2 + erc2cs(24-i)
        enddo

        erfc = (exp(-yt)/y)*(UM + b0-b2) / 2

      else

        z = 8*UM/yt - UM
        b1 = ZERO
        b0 = ZERO
        twoz = 2*z

        do i = 1,24
          b2 = b1
          b1 = b0
          b0 = twoz*b1 -b2 + erfccs(25-i)
        enddo

        erfc = (exp(-yt)/y)*(UM + b0-b2) / 2

      endif

    endif

    if(x < ZERO) erfc = 2*UM - erfc

  endif

  return

end function erfc
