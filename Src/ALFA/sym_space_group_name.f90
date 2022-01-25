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

!>  Prints the possible space groups
!>
!>  \author       Jose Luis Martins
!>  \version      5.01
!>  \date         23 April 2021
!>  \copyright    GNU Public License v2

subroutine sym_space_group_name(ibravais, ntrans, mtrx, tnp)

! information taken from http://img.chem.ucl.ac.uk/sgp/large/sgp.htm

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  ibravais                         !<  bravais lattice (1=sc,2=fcc,3=bcc,4=st,5=ct,6=so,7=bco,8=fco,9=bco,10=hex,11=romb,12=sm,13=cm,14=tri)


  integer, intent(in)               ::  ntrans                           !<  number of symmetry operations in the factor group
  integer, intent(in)               ::  mtrx(3,3,48)                     !<  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the factor group
  real(REAL64), intent(in)          ::  tnp(3,48)                        !<  2*pi* i-th component (in lattice coordinates) of the fractional translation vector associated with the k-th symmetry operation of the factor group

! local variables

  logical        ::  linversion
  logical        ::  lfrac

  integer        ::  inv
  integer        ::  mtmp(3,3)

  integer        ::  io

! parameters

!   real(REAL64), parameter    ::  UM = 1.0_REAL64
!   real(REAL64), parameter    ::  ZERO = 0.0_REAL64
  real(REAL64), parameter    ::  EPS = 1.0E-6_REAL64

! counters

  integer        ::  i, j

! checks if inversion is present

  io = 6

  linversion = .FALSE.
  if(ntrans > 1) then
    linversion = .TRUE.
    inv = ntrans/2+1
    do i = 1,3
      do j = 1,3
        mtmp(i,j) = mtrx(i,j,inv)
      enddo
      mtmp(i,i) = mtmp(i,i) + 1
    enddo
    do i = 1,3
    do j = 1,3
      if(mtmp(i,j) /= 0) linversion = .FALSE.
    enddo
    enddo
  endif

! checks if fractional translations are present ( beware: origin may be badly picked)

  lfrac = .FALSE.

  do j = 1,ntrans
    do i = 1,3
      if(abs(tnp(i,j)) > EPS) lfrac = .TRUE.
    enddo
  enddo

  write(io,*)
  if(ibravais == 1) then

    write(io,*) '  The Bravais lattice is  simple cubic'
    write(io,*)
    if(linversion) then

      if(ntrans == 24) then
        if(lfrac) then
          write(io,*) '  The space group could be:'
          write(io,*) '  Pn-3 (201), Pa-3 (205)'
          write(io,*) '  Pm-3 (200) if origin is badly chosen'
        else
          write(io,*) '  The space group could be:'
          write(io,*) '  Pm-3 (200)'
        endif
      elseif(ntrans == 48) then
        if(lfrac) then
          write(io,*) '  The space group could be:'
          write(io,*) '  Pn-3n (222), Pm-3n (223), Pn-3m (224)'
          write(io,*) '  Pm-3m (221) if origin is badly chosen'
        else
          write(io,*) '  The space group could be:'
          write(io,*) '  Pm-3m (221)'
        endif
      endif

    else

      if(ntrans == 12) then
        if(lfrac) then
          write(io,*) '  The space group could be:'
          write(io,*) '  P2_13 (198)'
          write(io,*) '  P23 (195) if origin is badly chosen'
        else
          write(io,*) '  The space group could be:'
          write(io,*) '  P23 (195)'
        endif
      elseif(ntrans == 24) then
        if(lfrac) then
          write(io,*) '  The space group could be:'
          write(io,*) '  P4_232 (208), P4_232 (208), P4_332 (212)'
          write(io,*) '  P4_132 (213), P-43n (218)'
          write(io,*) '  P432 (207), P-43m (215) if origin is badly chosen'
        else
          write(io,*) '  The space group could be:'
          write(io,*) '  P432 (207), P-43m (215)'
        endif
      endif

    endif

  elseif(ibravais == 2) then

    write(io,*) '  The Bravais lattice is fcc  face centered cubic'
    write(io,*)
    if(linversion) then

      if(ntrans == 24) then
        if(lfrac) then
          write(io,*) '  The space group could be:'
          write(io,*) '  Fd-3 (203)'
          write(io,*) '  Fm-3 (202) if origin is badly chosen'
        else
          write(io,*) '  The space group could be:'
          write(io,*) '  Fm-3 (202)'
        endif
      elseif(ntrans == 48) then
        if(lfrac) then
          write(io,*) '  The space group could be:'
          write(io,*) '  Fm-3c (226), Fd-3m (227), Fd-3c (228)'
          write(io,*) '  Fm-3m (225) if origin is badly chosen'
        else
          write(io,*) '  The space group could be:'
          write(io,*) '  Fm-3m (225)'
        endif
      endif

    else

      if(ntrans == 12) then
        if(lfrac) then
          write(io,*) '  The space group could be:'
          write(io,*) '  F23 (196) if origin is badly chosen'
        else
          write(io,*) '  The space group could be:'
          write(io,*) '  F23 (196)'
        endif
      elseif(ntrans == 24) then
        if(lfrac) then
          write(io,*) '  The space group could be:'
          write(io,*) '  F4_132 (210), F-43c (219)'
          write(io,*) '  F432 (209), F-432 (216) if origin is badly chosen'
        else
          write(io,*) '  The space group could be:'
          write(io,*) '  F432 (209), F-432 (216)'
        endif
      endif

    endif

  elseif(ibravais == 3) then

    write(io,*) '  The Bravais lattice is bcc  body centered cubic'
    write(io,*)
    if(linversion) then

      if(ntrans == 24) then
        if(lfrac) then
          write(io,*) '  The space group could be:'
          write(io,*) '  Ia-3 (206)'
          write(io,*) '  Im-3 (204) if origin is badly chosen'
        else
          write(io,*) '  The space group could be:'
          write(io,*) '  Im-3 (204)'
        endif
      elseif(ntrans == 48) then
        if(lfrac) then
          write(io,*) '  The space group could be:'
          write(io,*) '  Ia-3d (230)'
          write(io,*) '  Im-3m (229) if origin is badly chosen'
        else
          write(io,*) '  The space group could be:'
          write(io,*) '  Im-3m (229)'
        endif
      endif

    else

      if(ntrans == 12) then
        if(lfrac) then
          write(io,*) '  The space group could be:'
          write(io,*) '  I2_13 (198)'
          write(io,*) '  I23 (197) if origin is badly chosen'
        else
          write(io,*) '  The space group could be:'
          write(io,*) '  I23 (197)'
        endif
      elseif(ntrans == 24) then
        if(lfrac) then
          write(io,*) '  The space group could be:'
          write(io,*) '  I4_132 (214), I-43d (220)'
          write(io,*) '  I432 (211), I-43m (217) if origin is badly chosen'
        else
          write(io,*) '  The space group could be:'
          write(io,*) '  I432 (211), I-43m (217)'
        endif
      endif

    endif

  elseif(ibravais == 4) then

    write(io,*) '  The Bravais lattice is simple tetragonal'
    write(io,*)
    if(linversion) then

      if(ntrans == 8) then
        if(lfrac) then
          write(io,*) '  The space group could be:'
          write(io,*) '  P4_2/m (84), P4/n (85), P4_2/n (86)'
          write(io,*) '  P4/m (83) if origin is badly chosen'
        else
          write(io,*) '  The space group could be:'
        endif
      elseif(ntrans == 16) then
        if(lfrac) then
          write(io,*) '  The space group could be:'
          write(io,*) '  P4/mcc (124), P4/nbm (125), P4/nnc (126)'
          write(io,*) '  P4/mbm (127), P4/mnc (128), P4/nmm (129)'
          write(io,*) '  P4/ncc (130), P4_2/mmc (131), P4_2/mcm (132)'
          write(io,*) '  P4_2/nbc (133), P4_2/nnm (134), P4_2/mbc (135)'
          write(io,*) '  P4_2/mnm (136), P4_2/nmc (137), P4_2/ncm (138)'
          write(io,*) '  P4/mmm (123) if origin is badly chosen'
        else
          write(io,*) '  The space group could be:'
          write(io,*) '  P4/mmm (123)'
        endif
      endif

    else

      if(ntrans == 4) then
        if(lfrac) then
          write(io,*) '  The space group could be:'
          write(io,*) '  P4_1 (76), P4_2 (77), P4_3 (78)'
          write(io,*) '  P4 (75), P-4 (81) if origin is badly chosen'
        else
          write(io,*) '  The space group could be:'
          write(io,*) '  P4 (75), P-4 (81)'
        endif
      elseif(ntrans == 8) then
        if(lfrac) then
          write(io,*) '  The space group could be:'
          write(io,*) '  P42_12 (90), P4_122 (91), P4_2_12 (92), P4_222 (93)'
          write(io,*) '  P42_22_12 (94), P4_322 (95), P4_32_12 (96)'
          write(io,*) '  P4_322 (95), P4bm (100), P4-2cm (101)'
          write(io,*) '  P4_2nm (102), P4cc (103), P4nc (104)'
          write(io,*) '  P4_2mc(105), P4_2bc(106), P-42c (112)'
          write(io,*) '  P-42_1m (113), P-42_1c (114), P-4c2 (116)'
          write(io,*) '  P-4b2 (117), P-4n2 (118)'
          write(io,*) '  P422 (89), P4mm (99), P-42m (111) if origin is badly chosen'
          write(io,*) '  P-4m2 (115) if origin is badly chosen'
       else
          write(io,*) '  The space group could be:'
          write(io,*) '  P422 (89), P4mm (99), P-42m (111)'
          write(io,*) '  P-4m2 (115)'
        endif
      endif

    endif

  elseif(ibravais == 5) then

    write(io,*) '  The Bravais lattice is centered tetragonal'
    write(io,*)
    if(linversion) then

      if(ntrans == 8) then
        if(lfrac) then
          write(io,*) '  The space group could be:'
          write(io,*) '  I4_1/a (88)'
          write(io,*) '  I4/m (87) if origin is badly chosen'
        else
          write(io,*) '  The space group could be:'
          write(io,*) '  I4/m (87)'
        endif
      elseif(ntrans == 16) then
        if(lfrac) then
          write(io,*) '  The space group could be:'
          write(io,*) '  I4/mcm (140), I4_1/amd, I4_1/acd'
          write(io,*) '  I4/mmm (139) if origin is badly chosen'
        else
          write(io,*) '  The space group could be:'
          write(io,*) '  I4/mmm (139)'
        endif
      endif

    else

      if(ntrans == 4) then
        if(lfrac) then
          write(io,*) '  The space group could be:'
          write(io,*) '  I4_1 (80)'
          write(io,*) '  I4 (79), I-4 (82) if origin is badly chosen'
        else
          write(io,*) '  The space group could be:'
          write(io,*) '  I4 (79), I-4 (82)'
        endif
      elseif(ntrans == 8) then
        if(lfrac) then
          write(io,*) '  The space group could be:'
          write(io,*) '  I4_122 (98), I4cm (108), I4_1md(109)'
          write(io,*) '  I4_1cd (110), I-4c2 (120), I-42m (121)'
          write(io,*) '  I-42d (122)'
          write(io,*) '  I422 (97), I4mm (107), I-4m2 (119) if origin is badly chosen'
        else
          write(io,*) '  The space group could be:'
          write(io,*) '  I422 (97), I4mm (107), I-4m2 (119)'
        endif
      endif

    endif

  elseif(ibravais == 6) then

    write(io,*) '  The Bravais lattice is simple orthorhombic'
    write(io,*)
    if(linversion) then

      if(ntrans == 8) then
        if(lfrac) then
          write(io,*) '  The space group could be:'
          write(io,*) '  Pnnn (48), Pccm (49), Pban (50)'
          write(io,*) '  Pmma (51), Pnna (52), Pmna (53)'
          write(io,*) '  Pcca (54), Pbam (55), Pccn (56)'
          write(io,*) '  Pbcm (57), Pnnm (58), Pmmn (59)'
          write(io,*) '  Pbcn (60), Pbca (61), Pnma (62)'
          write(io,*) '  Pmmm (47) if origin is badly chosen'
        else
          write(io,*) '  The space group could be:'
          write(io,*) '  Pmmm (47)'
        endif
      endif

    else

      if(ntrans == 4) then
        if(lfrac) then
          write(io,*) '  The space group could be:'
          write(io,*) '  P222_1 (17), P2_12_12 (18), P2_12_12_1 (19)'
          write(io,*) '  Pmc2_1 (26), Pcc2 (27), Pma2 (28)'
          write(io,*) '  Pca2_1 (29), Pnc2 (30), Pmn2_1 (31)'
          write(io,*) '  Pba2 (32), Pna2_1 (33), Pnn2 (34)'
          write(io,*) '  P222 (16), Pmm2 (25) if origin is badly chosen'
        else
          write(io,*) '  The space group could be:'
          write(io,*) '  P222 (16), Pmm2 (25)'
        endif
      endif

    endif

  elseif(ibravais == 7) then

    write(io,*) '  The Bravais lattice is base centered orthorhombic'
    write(io,*)
    if(linversion) then

      if(ntrans == 8) then
        if(lfrac) then
          write(io,*) '  The space group could be:'
          write(io,*) '  Cmcm (63), Cmca (64), Cccm (66)'
          write(io,*) '  Cmma (67), Ccca (68)'
          write(io,*) '  Cmmm (65) if origin is badly chosen'
        else
          write(io,*) '  The space group could be:'
          write(io,*) '  Cmmm (65)'
        endif
      endif

    else

      if(ntrans == 4) then
        if(lfrac) then
          write(io,*) '  The space group could be:'
          write(io,*) '  C222_1 (20), Cmc2_1 (36), Ccc2 (37)'
          write(io,*) '  Abm2 (39), Ama2 (40), Aba2 (41)'
          write(io,*) '  C222 (21), Cmm2 (35), Amm2 (38) if origin is badly chosen'
        else
          write(io,*) '  The space group could be:'
          write(io,*) '  C222 (21), Cmm2 (35), Amm2 (38)'
        endif
      endif

    endif

  elseif(ibravais == 8) then

    write(io,*) '  The Bravais lattice is face centered orthorhombic'
    write(io,*)
    if(linversion) then

      if(ntrans == 8) then
        if(lfrac) then
          write(io,*) '  The space group could be:'
          write(io,*) '  Fddd (70)'
          write(io,*) '  Fmmm (69) if origin is badly chosen'
        else
          write(io,*) '  The space group could be:'
          write(io,*) '  Fmmm (69)'
        endif
      endif

    else

      if(ntrans == 4) then
        if(lfrac) then
          write(io,*) '  The space group could be:'
          write(io,*) '  Fdd2 (43)'
          write(io,*) '  F222 (22), Fmm2 (42) if origin is badly chosen'
        else
          write(io,*) '  The space group could be:'
          write(io,*) '  F222 (22), Fmm2 (42)'
        endif
      endif

    endif

  elseif(ibravais == 9) then

    write(io,*) '  The Bravais lattice is body centered orthorhombic'
    write(io,*)
    if(linversion) then

      if(ntrans == 8) then
        if(lfrac) then
          write(io,*) '  The space group could be:'
          write(io,*) '  Ibam (72), Ibca (73), Imma (74)'
          write(io,*) '  Immm (71) if origin is badly chosen'
        else
          write(io,*) '  The space group could be:'
          write(io,*) '  Immm (71)'
        endif
      endif

    else

      if(ntrans == 4) then
        if(lfrac) then
          write(io,*) '  The space group could be:'
          write(io,*) '  I2_12_12_1 (24), Iba2(45), Ima2 (46)'
          write(io,*) '  I222 (23), Imm2 (44) if origin is badly chosen'
        else
          write(io,*) '  The space group could be:'
          write(io,*) '  I222 (23), Imm2 (44)'
        endif
      endif

    endif

  elseif(ibravais == 10) then

    write(io,*) '  The Bravais lattice is hexagonal'
    write(io,*) '  including non-rhombohedral trigonal'
    write(io,*)
    if(linversion) then

      if(ntrans == 6) then
        if(lfrac) then
          write(io,*) '  The space group could be:'
          write(io,*) '  P-3 (147) if origin is badly chosen'
        else
          write(io,*) '  The space group could be:'
          write(io,*) '  P-3 (147)'
        endif
      elseif(ntrans == 12) then
        if(lfrac) then
          write(io,*) '  The space group could be:'
          write(io,*) '  P-31c (163), P-3c1 (165), P6_3/m (176)'
          write(io,*) '  P-31m (162), P-3m1 (164), P6/m (175) if origin is badly chosen'
        else
          write(io,*) '  The space group could be:'
          write(io,*) '  P-31m (162), P-3m1 (164), P6/m (175)'
        endif
      elseif(ntrans == 24) then
        if(lfrac) then
          write(io,*) '  The space group could be:'
          write(io,*) '  P6/mmc (192), P6_3/mcm (193), P6_3/mmc (194)'
          write(io,*) '  P6/mmm (191) if origin is badly chosen'
        else
          write(io,*) '  The space group could be:'
          write(io,*) '  P6/mmm (191)'
        endif
      endif

    else

      if(ntrans == 3) then
        if(lfrac) then
          write(io,*) '  The space group could be:'
          write(io,*) '  P3_1 (144), P3_2 (145)'
          write(io,*) '  P3 (143) if origin is badly chosen'
        else
          write(io,*) '  The space group could be:'
          write(io,*) '  P3 (143)'
        endif
      elseif(ntrans == 6) then
        if(lfrac) then
          write(io,*) '  The space group could be:'
          write(io,*) '  P3_112 (151), P3_121 (152), P3_212 (153)'
          write(io,*) '  P3_221 (154), P3c1 (158), P31c (159)'
          write(io,*) '  P6_1 (169), P6_5 (170), P6_2 (171)'
          write(io,*) '  P6_4 (172), P6_3 (173)'
          write(io,*) '  P312 (149), P321 (150), P3m1 (156) if origin is badly chosen'
          write(io,*) '  P31m (157), P6 (168), P-6 (174) if origin is badly chosen'
        else
          write(io,*) '  The space group could be:'
          write(io,*) '  P312 (149), P321 (150), P3m1 (156)'
          write(io,*) '  P31m (157), P6 (168), P-6 (174)'
        endif
      elseif(ntrans == 12) then
        if(lfrac) then
          write(io,*) '  The space group could be:'
          write(io,*) '  P6_122 (178), P6_522 (179), P6_222 (180)'
          write(io,*) '  P6_422 (181), P6_322 (182), P6cc (184)'
          write(io,*) '  P6_3cm (185), P6_3mc (186), P-6c2 (188)'
          write(io,*) '  P-62c (190)'
          write(io,*) '  P622 (177), P6mm (183), P-6m2 (187) if origin is badly chosen'
          write(io,*) '  P-62m (189) if origin is badly chosen'
        else
          write(io,*) '  The space group could be:'
          write(io,*) '  P622 (177), P6mm (183), P-6m2 (187)'
          write(io,*) '  P-62m (189)'
        endif
      endif

    endif

  elseif(ibravais == 11) then

    write(io,*) '  The Bravais lattice is rhombohedral'
    write(io,*)
    if(linversion) then

      if(ntrans == 6) then
        if(lfrac) then
          write(io,*) '  The space group could be:'
          write(io,*) '  R-3 (148) if origin is badly chosen'
        else
          write(io,*) '  The space group could be:'
          write(io,*) '  R-3 (148)'
        endif
      elseif(ntrans == 12) then
        if(lfrac) then
          write(io,*) '  The space group could be:'
          write(io,*) '  R-3c (167)'
          write(io,*) '  R-3m (166) if origin is badly chosen'
        else
          write(io,*) '  The space group could be:'
          write(io,*) '  R-3m (166)'
        endif
      endif

    else

      if(ntrans == 3) then
        if(lfrac) then
          write(io,*) '  The space group could be:'
          write(io,*) '  R3 (146) if origin is badly chosen'
        else
          write(io,*) '  The space group could be:'
          write(io,*) '  R3 (146)'
        endif
      elseif(ntrans == 6) then
        if(lfrac) then
          write(io,*) '  The space group could be:'
          write(io,*) '  R3c (161)'
          write(io,*) '  R32 (155), R3m (160) if origin is badly chosen'
        else
          write(io,*) '  The space group could be:'
          write(io,*) '  R32 (155), R3m (160)'
        endif
      endif

    endif

  elseif(ibravais == 12) then

    write(io,*) '  The Bravais lattice is simple monoclinic'
    write(io,*)
    if(linversion) then

      if(ntrans == 4) then
        if(lfrac) then
          write(io,*) '  The space group could be:'
          write(io,*) '  P2_1/m (11), P2/c (13), P2_1/c'
          write(io,*) '  P2/m (10) if origin is badly chosen'
        else
          write(io,*) '  The space group could be:'
          write(io,*) '  P2/m (10)'
        endif
      endif

    else

      if(ntrans == 2) then
        if(lfrac) then
          write(io,*) '  The space group could be:'
          write(io,*) '  P2_1 (4), Pc (7)'
          write(io,*) '  P2 (3), Pm (6) if origin is badly chosen'
        else
          write(io,*) '  The space group could be:'
          write(io,*) '  P2 (3), Pm (6)'
        endif
      endif

    endif

  elseif(ibravais == 13) then

    write(io,*) '  The Bravais lattice is centered monoclinic'
    write(io,*)
    if(linversion) then

      if(ntrans == 4) then
        if(lfrac) then
          write(io,*) '  The space group could be:'
          write(io,*) '  C2/c (15)'
          write(io,*) '  C2/m (12) if origin is badly chosen'
        else
          write(io,*) '  The space group could be:'
          write(io,*) '  C2/m (12)'
        endif
      endif

    else

      if(ntrans == 2) then
        if(lfrac) then
          write(io,*) '  The space group could be:'
          write(io,*) '  Cc (9)'
          write(io,*) '  C2 (5), Cm (8) if origin is badly chosen'
        else
          write(io,*) '  The space group could be:'
          write(io,*) '  C2 (5), Cm (8)'
        endif
      endif

    endif

  elseif(ibravais == 14) then

    write(io,*) '  The Bravais lattice is triclinic'
    write(io,*)
    if(linversion) then

          write(io,*) '  The space group is:'
          write(io,*) '  P-1 (2)'

    else

          write(io,*) '  The space is:'
          write(io,*) '  P1 (1)'

    endif

  endif

  write(io,*)
  write(io,*) '  As not all 230 possible space groups were tested,'
  write(io,*) '  check this space group information carefuly.'
  write(io,*)

  return
end subroutine sym_space_group_name
