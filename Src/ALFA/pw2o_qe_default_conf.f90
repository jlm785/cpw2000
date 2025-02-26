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

!>  subroutine provides a default electronic configuration for each atom.
!>
!>  \author       Nikolas Garofil
!>  \version      5.11
!>  \date         12 February 2021.
!>  \copyright    GNU Public License v2

subroutine pw2o_qe_default_conf(ized,config)

! Stolen from Nikolas Garofil, Quantum espresso, 12 February 2021. JLM

  implicit none

  integer, parameter  :: REAL64 = selected_real_kind(12)


  integer, intent(in)                      :: ized                       !<  atomic number
  character(len=80), intent(out)           ::  config                    !<  default atomic configuration

  config = 'unknown'
  if (ized==1)   config='1s1.0'
  if (ized==2)   config='1s2.0'
  if (ized==3)   config='[He] 2s1.0'
  if (ized==4)   config='[He] 2s2.0'
  if (ized==5)   config='[He] 2s2.0 2p1.0'
  if (ized==6)   config='[He] 2s2.0 2p2.0'
  if (ized==7)   config='[He] 2s2.0 2p3.0'
  if (ized==8)   config='[He] 2s2.0 2p4.0'
  if (ized==9)   config='[He] 2s2.0 2p5.0'
  if (ized==10)  config='[He] 2s2.0 2p6.0'
  if (ized==11)  config='[Ne] 3s1.0'
  if (ized==12)  config='[Ne] 3s2.0'
  if (ized==13)  config='[Ne] 3s2.0 3p1.0'
  if (ized==14)  config='[Ne] 3s2.0 3p2.0'
  if (ized==15)  config='[Ne] 3s2.0 3p3.0'
  if (ized==16)  config='[Ne] 3s2.0 3p4.0'
  if (ized==17)  config='[Ne] 3s2.0 3p5.0'
  if (ized==18)  config='[Ne] 3s2.0 3p6.0'
  if (ized==19)  config='[Ar] 4s1.0'
  if (ized==20)  config='[Ar] 4s2.0'
  if (ized==21)  config='[Ar] 4s2.0 3d1.0'
  if (ized==22)  config='[Ar] 4s2.0 3d2.0'
  if (ized==23)  config='[Ar] 4s2.0 3d3.0'
  if (ized==24)  config='[Ar] 4s1.0 3d5.0'
  if (ized==25)  config='[Ar] 4s2.0 3d5.0'
  if (ized==26)  config='[Ar] 4s2.0 3d6.0'
  if (ized==27)  config='[Ar] 4s2.0 3d7.0'
  if (ized==28)  config='[Ar] 4s2.0 3d8.0'
  if (ized==29)  config='[Ar] 4s1.0 3d10.0'
  if (ized==30)  config='[Ar] 4s2.0 3d10.0'
  if (ized==31)  config='[Ar] 4s2.0 3d10.0 4p1.0'
  if (ized==32)  config='[Ar] 4s2.0 3d10.0 4p2.0'
  if (ized==33)  config='[Ar] 4s2.0 3d10.0 4p3.0'
  if (ized==34)  config='[Ar] 4s2.0 3d10.0 4p4.0'
  if (ized==35)  config='[Ar] 4s2.0 3d10.0 4p5.0'
  if (ized==36)  config='[Ar] 4s2.0 3d10.0 4p6.0'
  if (ized==37)  config='[Kr] 5s1.0'
  if (ized==38)  config='[Kr] 5s2.0'
  if (ized==39)  config='[Kr] 5s2.0 4d1.0'
  if (ized==40)  config='[Kr] 5s2.0 4d2.0'
  if (ized==41)  config='[Kr] 5s1.0 4d4.0'
  if (ized==42)  config='[Kr] 5s1.0 4d5.0'
  if (ized==43)  config='[Kr] 5s2.0 4d5.0'
  if (ized==44)  config='[Kr] 5s1.0 4d7.0'
  if (ized==45)  config='[Kr] 5s1.0 4d8.0'
  if (ized==46)  config='[Kr] 5s0.0 4d10.0'
  if (ized==47)  config='[Kr] 5s1.0 4d10.0'
  if (ized==48)  config='[Kr] 5s2.0 4d10.0'
  if (ized==49)  config='[Kr] 5s2.0 4d10.0 5p1.0'
  if (ized==50)  config='[Kr] 5s2.0 4d10.0 5p2.0'
  if (ized==51)  config='[Kr] 5s2.0 4d10.0 5p3.0'
  if (ized==52)  config='[Kr] 5s2.0 4d10.0 5p4.0'
  if (ized==53)  config='[Kr] 5s2.0 4d10.0 5p5.0'
  if (ized==54)  config='[Kr] 5s2.0 4d10.0 5p6.0'
  if (ized==55)  config='[Xe] 6s1.0'
  if (ized==56)  config='[Xe] 6s2.0'
  if (ized==57)  config='[Xe] 6s2.0 5d1.0 4f0.0'
  if (ized==58)  config='[Xe] 6s2.0 5d1.0 4f1.0'
  if (ized==59)  config='[Xe] 6s2.0 5d0.0 4f3.0'
  if (ized==60)  config='[Xe] 6s2.0 5d0.0 4f4.0'
  if (ized==61)  config='[Xe] 6s2.0 5d0.0 4f5.0'
  if (ized==62)  config='[Xe] 6s2.0 5d0.0 4f6.0'
  if (ized==63)  config='[Xe] 6s2.0 5d0.0 4f7.0'
  if (ized==64)  config='[Xe] 6s2.0 5d1.0 4f7.0'
  if (ized==65)  config='[Xe] 6s2.0 5d0.0 4f9.0'
  if (ized==66)  config='[Xe] 6s2.0 5d0.0 4f10.0'
  if (ized==67)  config='[Xe] 6s2.0 5d0.0 4f11.0'
  if (ized==68)  config='[Xe] 6s2.0 5d0.0 4f12.0'
  if (ized==69)  config='[Xe] 6s2.0 5d0.0 4f13.0'
  if (ized==70)  config='[Xe] 6s2.0 5d0.0 4f14.0'
  if (ized==71)  config='[Xe] 6s2.0 5d1.0 4f14.0'
  if (ized==72)  config='[Xe] 6s2.0 5d2.0 4f14.0'
  if (ized==73)  config='[Xe] 6s2.0 5d3.0 4f14.0'
  if (ized==74)  config='[Xe] 6s2.0 5d4.0 4f14.0'
  if (ized==75)  config='[Xe] 6s2.0 5d5.0 4f14.0'
  if (ized==76)  config='[Xe] 6s2.0 5d6.0 4f14.0'
  if (ized==77)  config='[Xe] 6s2.0 5d7.0 4f14.0'
  if (ized==78)  config='[Xe] 6s1.0 5d9.0 4f14.0'
  if (ized==79)  config='[Xe] 6s1.0 5d10.0 4f14.0'
  if (ized==80)  config='[Xe] 6s2.0 5d10.0 4f14.0'
  if (ized==81)  config='[Xe] 6s2.0 5d10.0 4f14.0 6p1.0'
  if (ized==82)  config='[Xe] 6s2.0 5d10.0 4f14.0 6p2.0'
  if (ized==83)  config='[Xe] 6s2.0 5d10.0 4f14.0 6p3.0'
  if (ized==84)  config='[Xe] 6s2.0 5d10.0 4f14.0 6p4.0'
  if (ized==85)  config='[Xe] 6s2.0 5d10.0 4f14.0 6p5.0'
  if (ized==86)  config='[Xe] 6s2.0 5d10.0 4f14.0 6p6.0'
  if (ized==87)  config='[Rn] 7s1.0'
  if (ized==88)  config='[Rn] 7s2.0'
  if (ized==89)  config='[Rn] 7s2.0 6d1.0 5f0.0'
  if (ized==90)  config='[Rn] 7s2.0 6d2.0 5f0.0'
  if (ized==91)  config='[Rn] 7s2.0 6d1.0 5f2.0'
  if (ized==92)  config='[Rn] 7s2.0 6d1.0 5f3.0'
  if (ized==93)  config='[Rn] 7s2.0 6d1.0 5f4.0'
  if (ized==94)  config='[Rn] 7s2.0 6d0.0 5f6.0'
  if (ized==95)  config='[Rn] 7s2.0 6d0.0 5f7.0'
  if (ized==96)  config='[Rn] 7s2.0 6d1.0 5f7.0'
  if (ized==97)  config='[Rn] 7s2.0 6d0.0 5f9.0'
  if (ized==98)  config='[Rn] 7s2.0 6d0.0 5f10.0'
  if (ized==99)  config='[Rn] 7s2.0 6d0.0 5f11.0'
  if (ized==100) config='[Rn] 7s2.0 6d0.0 5f12.0'
  if (ized==101) config='[Rn] 7s2.0 6d0.0 5f13.0'
  if (ized==102) config='[Rn] 7s2.0 6d0.0 5f14.0'
  if (ized==103) config='[Rn] 7s2.0 7p1.0 5f14.0'

  return

end subroutine pw2o_qe_default_conf
