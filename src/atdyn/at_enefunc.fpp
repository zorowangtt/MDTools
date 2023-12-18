!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_enefunc_mod
!> @brief   define potential energy functions
!! @authors Yuji Sugita (YS), Takaharu Mori (TM), Chigusa Kobayashi (CK)
!!          Kiyoshi Yagi (KY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_enefunc_mod

  use at_enefunc_gromacs_mod
  use at_enefunc_amber_mod
  use at_enefunc_charmm_mod
  use at_enefunc_go_mod
  use at_energy_mod
  use at_restraints_str_mod
  use at_enefunc_str_mod
  use at_energy_str_mod
  use at_boundary_str_mod
  use molecules_str_mod
  use select_mod
  use select_atoms_mod
  use select_atoms_str_mod
  use fitting_mod
  use fitting_str_mod
  use fileio_grotop_mod
  use fileio_prmtop_mod
  use fileio_par_mod
  use fileio_gpr_mod
  use fileio_eef1_mod
  use messages_mod
  use constants_mod
  use mpi_parallel_mod

  implicit none
  private

  ! subroutines
  public  :: define_enefunc
  public  :: setup_fitting_atdyn
  private :: setup_enefunc_dispcorr
  private :: setup_enefunc_fit_refcoord
  private :: remove_fixatm
  private :: remove_fixatm_bond
  private :: remove_fixatm_angl
  private :: remove_fixatm_urey
  private :: remove_fixatm_dihe
  private :: remove_fixatm_rb_dihe
  private :: remove_fixatm_impr
  private :: remove_fixatm_cmap

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    define_enefunc
  !> @brief        a driver subroutine for defining potential energy functions
  !! @authors      TM, CK
  !! @param[in]    ene_info   : ENERGY section control parameters
  !! @param[in]    boundary   : boundary conditions information
  !! @param[in]    par        : CHARMM PAR information
  !! @param[in]    gpr        : GO model parameter information
  !! @param[in]    prmtop     : AMBER parameter topology information
  !! @param[in]    grotop     : GROMACS parameter topology information
  !! @param[in]    molecule   : molecule information
  !! @param[in]    restraints : restraints information
  !! @param[out]   enefunc    : structure of enefunc
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine define_enefunc(ene_info, boundary, par, gpr, prmtop, grotop, eef1,&
                            molecule, restraints, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info 
    type(s_boundary),        intent(in)    :: boundary
    type(s_par),             intent(in)    :: par
    type(s_gpr),             intent(in)    :: gpr
    type(s_prmtop),          intent(in)    :: prmtop
    type(s_grotop),          intent(in)    :: grotop
    type(s_eef1),            intent(in)    :: eef1
    type(s_molecule),        intent(in)    :: molecule
    type(s_restraints),      intent(in)    :: restraints
    type(s_enefunc),         intent(inout) :: enefunc


    call init_enefunc(enefunc)


    enefunc%forcefield      = ene_info%forcefield
    enefunc%output_style    = ene_info%output_style

    enefunc%switchdist      = ene_info%switchdist
    enefunc%cutoffdist      = ene_info%cutoffdist
    enefunc%pairlistdist    = ene_info%pairlistdist
    enefunc%dielec_const    = ene_info%dielec_const
    enefunc%force_switch    = ene_info%vdw_force_switch
    enefunc%vdw_shift       = ene_info%vdw_shift
    enefunc%dispersion_corr = ene_info%dispersion_corr
    enefunc%nonb_limiter    = ene_info%nonb_limiter
    enefunc%contact_check   = ene_info%contact_check
    enefunc%minimum_contact = ene_info%minimum_contact


    if (par%num_bonds > 0) then

      call define_enefunc_charmm &
                         (ene_info, boundary, par, eef1, molecule, &
                          restraints, enefunc)

    else if (gpr%num_bonds > 0) then

      call define_enefunc_go &
                         (ene_info, gpr, molecule, restraints, enefunc)

    else if (prmtop%num_atoms > 0) then

      call define_enefunc_amber &
                          (ene_info, boundary, prmtop, molecule, &
                           restraints, enefunc)

    else if (grotop%num_atomtypes > 0) then

      call define_enefunc_gromacs &
                          (ene_info, boundary, grotop, molecule, &
                           restraints, enefunc)

    end if

    ! remove functions for fixed atoms
    !
    if (boundary%num_fixatm > 0) &
      call remove_fixatm(boundary, enefunc)

    ! dispersion correction
    !
    call setup_enefunc_dispcorr(ene_info, molecule, enefunc)

    ! spherical potential
    !
    if (boundary%sph_pot) then
      enefunc%spot_use = .true.

    else
      enefunc%spot_use = .false.

    end if

    return

  end subroutine define_enefunc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_fitting_atdyn
  !> @brief        setup fitting in atdyn
  !! @authors      CK
  !! @param[in]    is_rpath   : flag if rpath or not
  !! @param[in]    fit_info   : fitting information
  !! @param[in]    sel_info   : selection information
  !! @param[inout] molecule : molecule information
  !! @param[inout] enefunc  : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  subroutine setup_fitting_atdyn(is_rpath, fit_info, sel_info, &
                                 molecule, enefunc)

    ! formal arguments
    logical,                 intent(in)      :: is_rpath
    type(s_fit_info),        intent(in)      :: fit_info
    type(s_sel_info),        intent(in)      :: sel_info
    type(s_molecule),        intent(in)      :: molecule
    type(s_enefunc), target, intent(inout)   :: enefunc

    ! local variables
    integer                                :: i
    integer                                :: fitting_atom_idx
    integer                                :: fitting_move
    logical                                :: fit_check, fit_check_rpath
    type(s_selatoms)                       :: selatoms

    enefunc%fitting_method  = fit_info%fitting_method

    if (enefunc%pressure_rmsd) then
      if (enefunc%fitting_method .ne. FittingMethodTR .and.      &
          enefunc%fitting_method .ne. FittingMethodTR_ROT .and.  &
          enefunc%fitting_method .ne. FittingMethodTR_ZROT) then
          call error_msg('Setup_Fitting_Atdyn> pressure_rmsd option is '//&
                      'allowed only in Translation fitting')
      endif
    endif

    if (enefunc%fitting_method .eq. FittingMethodNO) then
      if (main_rank) then
        fit_check = .false.
        do i = 1,enefunc%num_restraintfuncs
          if (enefunc%restraint_kind(i) .eq. RestraintsFuncRMSD .or. &
            enefunc%restraint_kind(i) .eq. RestraintsFuncRMSDCOM) then
            fit_check = .true.
          endif
        end do
        if (fit_check .and. .not. fit_info%force_no_fitting) then
          call error_msg('Setup_Fitting_Atdyn> No fit is not allowed '//&
                      'in RMSD restraint')
        endif
        if (fit_check .and. fit_info%force_no_fitting) then
          write(MsgOut,'(A)') "Setup_Fitting_Atdyn> RMSD restraint without FITTING"
          write(MsgOut,*) 
        endif
      endif
      do i = 1,enefunc%num_restraintfuncs
        if (enefunc%restraint_kind(i) .eq. RestraintsFuncPC .or. &
            enefunc%restraint_kind(i) .eq. RestraintsFuncPCCOM) then
           call error_msg('Setup_Fitting_Atdyn> No fit is not allowed '//&
                       'in PC/PCCOM ')
        endif
      end do
      return
    endif

    fit_check       = .false.
    do i = 1,enefunc%num_restraintfuncs
      if (enefunc%restraint_kind(i) .eq. RestraintsFuncRMSD .or. &
          enefunc%restraint_kind(i) .eq. RestraintsFuncRMSDCOM .or. &
          enefunc%restraint_kind(i) .eq. RestraintsFuncPC .or. &
          enefunc%restraint_kind(i) .eq. RestraintsFuncPCCOM) then
          fit_check = .true.
          exit
      endif
    end do

    fit_check_rpath = .false.
    do i = 1,enefunc%num_restraintfuncs
      if ( enefunc%restraint_kind(i) .eq. RestraintsFuncPOSI .and. &
           is_rpath) then
        fit_check_rpath = .true.
        exit
      endif
    end do

    if (.not. fit_check .and. .not. fit_check_rpath) then
       if (main_rank) then
         write(MsgOut,'(A)') "Setup_Fitting_Atdyn> NO fitting is applied, skip"
         write(MsgOut,'(A)') "  fitting method  =  NO"
         write(MsgOut,*) 
       endif
       enefunc%fitting_method=FittingMethodNO
       return
    endif

    if (enefunc%fitting_method .ne. FittingMethodTR_ROT .and. &
        enefunc%fitting_method .ne. FittingMethodXYTR_ZROT) &
      call error_msg('Setup_Fitting_Atdyn> NO/TR+ROT/XYTR+ZROT is allowed')

    fitting_atom_idx          = fit_info%fitting_atom

    call select_atom(molecule, sel_info%groups(fitting_atom_idx), selatoms)
    enefunc%num_fitting = size(selatoms%idx)
    enefunc%mass_weight   = fit_info%mass_weight

    if (is_rpath) then
      enefunc%fitting_file  = FittingFileFIT
      enefunc%fitting_move  = FittingMoveSYS
    else
      enefunc%fitting_file  = FittingFileREF
      enefunc%fitting_move  = FittingMoveREF
    endif

    do i = 1,enefunc%num_restraintfuncs
      if (enefunc%restraint_kind(i) .eq. RestraintsFuncRMSDCOM .and. &
        (.not. enefunc%mass_weight) ) &
        call error_msg('Setup_Fitting_Atdyn> RESTRAINTS and FITTING '//&
                       'is inconsistent. Please check [FITTING] section.')
      if (enefunc%restraint_kind(i) .eq. RestraintsFuncPCCOM .and. &
        (.not. enefunc%mass_weight) ) &
        call error_msg('Setup_Fitting_Atdyn> RESTRAINTS and FITTING '//&
                       'is inconsistent.  Please check [FITTING] section.')
      if (enefunc%restraint_kind(i) .eq. RestraintsFuncRMSD .and. &
        (enefunc%mass_weight) ) &
        call error_msg('Setup_Fitting_Atdyn> RESTRAINTS and FITTING'//&
                       'is inconsistent.  Please check [FITTING] section.')
      if (enefunc%restraint_kind(i) .eq. RestraintsFuncPC .and. &
        (enefunc%mass_weight) ) &
        call error_msg('Setup_Fitting_Atdyn> RESTRAINTS and FITTING'//&
                       'is inconsistent.  Please check [FITTING] section.')
      if (enefunc%restraint_kind(i) .eq. RestraintsFuncRMSD .or.    &
          enefunc%restraint_kind(i) .eq. RestraintsFuncRMSDCOM .or. &
          enefunc%restraint_kind(i) .eq. RestraintsFuncPC .or.      &
          enefunc%restraint_kind(i) .eq. RestraintsFuncPCCOM) then
        if (fitting_atom_idx .ne. enefunc%restraint_grouplist(1,i)) then
          call error_msg('Setup_Fitting_Atdyn> Fitting group should be'//&
                       'idenical to the group in restraint')
        endif
      endif
    end do

    call setup_enefunc_fit_refcoord(enefunc%fitting_file, molecule, enefunc)
    enefunc%fitting_atom(1:enefunc%num_fitting) =     &
        selatoms%idx(1:enefunc%num_fitting)

    return

  end subroutine setup_fitting_atdyn

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_dispcorr
  !> @brief        define dispersion correction term 
  !! @authors      CK
  !! @param[in]    ene_info : ENERGY section control parameters information
  !! @param[inout] molecule : molecule information
  !! @param[inout] enefunc  : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  subroutine setup_enefunc_dispcorr(ene_info, molecule, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)      :: ene_info
    type(s_molecule),        intent(in)      :: molecule
    type(s_enefunc), target, intent(inout)   :: enefunc

    ! local variables
    integer                  :: i, j, iatmcls, ntypes
    integer                  :: i1, i2, i3, i4
    integer                  :: natom2, nexpair
    real(wp)                 :: lj6_tot, lj6_diff, lj6_ex
    real(wp)                 :: factor, rpair
    real(wp)                 :: diff_cs, diff_cs2, diff_cs3, diff_cs4
    real(wp)                 :: cutoff , cutoff2, cutoff3, cutoff4
    real(wp)                 :: cutoff5, cutoff6, cutoff7, cutoff8
    real(wp)                 :: cutoff14
    real(wp)                 :: inv_cutoff3, inv_cutoff6, inv_cutoff12
    real(wp)                 :: switchdist , switchdist2, switchdist3
    real(wp)                 :: switchdist4, switchdist5
    real(wp)                 :: switchdist6, switchdist7, switchdist8
    real(wp)                 :: shift_a, shift_b, shift_c
    real(wp)                 :: vswitch, eswitch, vlong

    integer,         pointer :: bondlist(:,:),anglelist(:,:), ureylist(:,:)
    integer,         pointer :: dihelist(:,:),rb_dihelist(:,:)
    integer,         pointer :: imprlist(:,:)

    integer,     allocatable :: atype(:)


    if (ene_info%dispersion_corr == Disp_corr_NONE) &
      return

    bondlist    => enefunc%bond_list
    anglelist   => enefunc%angl_list
    ureylist    => enefunc%urey_list
    dihelist    => enefunc%dihe_list
    rb_dihelist => enefunc%rb_dihe_list
    imprlist    => enefunc%impr_list

    ntypes = enefunc%num_atom_cls
    allocate(atype(1:ntypes))

    atype(1:ntypes) = 0

    do i = 1, molecule%num_atoms
      iatmcls = molecule%atom_cls_no(i)
      atype(iatmcls) = atype(iatmcls)+1
    end do

    lj6_tot = 0.0_wp
    do i = 1, ntypes
      do j = 1, ntypes
        lj6_tot = lj6_tot + enefunc%nonb_lj6(i,j)*atype(i)*atype(j)
      end do
    end do

    deallocate(atype)

    cutoff      = enefunc%cutoffdist
    cutoff2     = cutoff*cutoff
    cutoff3     = cutoff2*cutoff
    inv_cutoff3 = 1.0_wp/cutoff3

    eswitch = 0.0_wp
    vswitch = 0.0_wp
    vlong   = inv_cutoff3/3.0_wp

    if (enefunc%forcefield == ForcefieldAMBER) then

      factor = 2.0_wp*PI*lj6_tot
      enefunc%dispersion_energy = -factor*vlong
      enefunc%dispersion_virial = -2.0_wp*factor*vlong

    else if (enefunc%forcefield == ForcefieldGROAMBER .or.             &
             enefunc%forcefield == ForcefieldGROMARTINI) then

      ! remove exclusion
      !
      lj6_ex = 0.0_wp
      nexpair = 0
      do i = 1, molecule%num_atoms
        iatmcls = molecule%atom_cls_no(i)
        lj6_ex  = lj6_ex + enefunc%nonb_lj6(iatmcls,iatmcls)
      end do

      ! bonds
      do i = 1, enefunc%num_bonds
        i1     = molecule%atom_cls_no(bondlist(1,i))
        i2     = molecule%atom_cls_no(bondlist(2,i))
        lj6_ex = lj6_ex + enefunc%nb14_lj6(i1,i2)
      end do

      ! angles
      do i = 1, enefunc%num_angles
        i1     = molecule%atom_cls_no(anglelist(1,i))
        i3     = molecule%atom_cls_no(anglelist(3,i))
        lj6_ex = lj6_ex + enefunc%nb14_lj6(i1,i3)
      end do

      ! Urey angles
      do i = 1, enefunc%num_ureys
        i1     = molecule%atom_cls_no(ureylist(1,i))
        i3     = molecule%atom_cls_no(ureylist(3,i))
        lj6_ex = lj6_ex + enefunc%nb14_lj6(i1,i3)
      end do

      ! dihedral
      do i = 1, enefunc%num_dihedrals
        i1     = molecule%atom_cls_no(dihelist(1,i))
        i4     = molecule%atom_cls_no(dihelist(4,i))
        lj6_ex = lj6_ex + enefunc%nb14_lj6(i1,i4)
      end do

      ! RB dihedral
      do i = 1, enefunc%num_rb_dihedrals
        i1     = molecule%atom_cls_no(rb_dihelist(1,i))
        i4     = molecule%atom_cls_no(rb_dihelist(4,i))
        lj6_ex = lj6_ex + enefunc%nb14_lj6(i1,i4)
      end do

      ! improper
      do i = 1, enefunc%num_impropers
        i1     = molecule%atom_cls_no(imprlist(1,i))
        i4     = molecule%atom_cls_no(imprlist(4,i))
        lj6_ex = lj6_ex + enefunc%nb14_lj6(i1,i4)
      end do

      nexpair = nexpair + molecule%num_atoms       &
                        + enefunc%num_bonds        &
                        + enefunc%num_angles       &
                        + enefunc%num_ureys        &
                        + enefunc%num_dihedrals    &
                        + enefunc%num_rb_dihedrals &
                        + enefunc%num_impropers

      lj6_diff = (lj6_tot - lj6_ex)

      natom2 = molecule%num_atoms*molecule%num_atoms
      rpair  = real(natom2/(natom2-nexpair),wp)
      factor       = 2.0_wp*PI*rpair*lj6_diff

      switchdist   = enefunc%switchdist
      diff_cs      = (cutoff - switchdist)

      if (diff_cs > EPS) then

        if (enefunc%vdw_shift) then

          cutoff4      = cutoff3*cutoff
          cutoff5      = cutoff4*cutoff
          cutoff6      = cutoff5*cutoff
          cutoff7      = cutoff6*cutoff
          cutoff8      = cutoff7*cutoff
          cutoff14     = cutoff7*cutoff7
          inv_cutoff6  = inv_cutoff3*inv_cutoff3
          inv_cutoff12 = inv_cutoff6*inv_cutoff6
  
          diff_cs2     = diff_cs*diff_cs
          diff_cs3     = diff_cs2*diff_cs
          diff_cs4     = diff_cs3*diff_cs
  
          switchdist2  = switchdist*switchdist
          switchdist3  = switchdist2*switchdist
          switchdist4  = switchdist3*switchdist
          switchdist5  = switchdist4*switchdist
          switchdist6  = switchdist5*switchdist
          switchdist7  = switchdist6*switchdist
          switchdist8  = switchdist7*switchdist
  
          ! LJ6
          ! 
          shift_a = -(10.0_wp*cutoff - 7.0_wp*switchdist)/(cutoff8*diff_cs2)
          shift_b =  ( 9.0_wp*cutoff - 7.0_wp*switchdist)/(cutoff8*diff_cs3)
  
          shift_c = inv_cutoff6 - 2.0_wp * shift_a * diff_cs3  &
                    - 1.5_wp * shift_b * diff_cs4
  
          eswitch = -2.0_wp * shift_a * ((1.0_wp/6.0_wp)*cutoff6               &
                                        -(3.0_wp/5.0_wp)*cutoff5*switchdist    &
                                        +(3.0_wp/4.0_wp)*cutoff4*switchdist2   &
                                        -(1.0_wp/3.0_wp)*cutoff3*switchdist3   &
                                        +(1.0_wp/6.0e1_wp)*switchdist6)        &
                    -1.5_wp * shift_b * ((1.0_wp/7.0_wp)*cutoff7               &
                                        -(2.0_wp/3.0_wp)*cutoff6*switchdist    &
                                        +(6.0_wp/5.0_wp)*cutoff5*switchdist2   &
                                        -                cutoff4*switchdist3   &
                                        +(1.0_wp/3.0_wp)*cutoff3*switchdist4   &
                                        -(1.0_wp/1.05e2_wp)*switchdist7)       &
                    -(1.0_wp/3.0_wp) * shift_c * (cutoff3-switchdist3)
  
          ! LJ12
          !
          shift_a = -(16.0_wp*cutoff - 13.0_wp*switchdist)/(cutoff14*diff_cs2)
          shift_b =  (15.0_wp*cutoff - 13.0_wp*switchdist)/(cutoff14*diff_cs3)
          shift_c = inv_cutoff12 - 2.0_wp * shift_a * diff_cs3  &
                    - 1.5_wp * shift_b * diff_cs4
  
  
          vswitch = shift_a * ( (1.0_wp/6.0_wp)*cutoff6                        &
                               -(2.0_wp/5.0_wp)*cutoff5*switchdist             &
                               +(1.0_wp/4.0_wp)*cutoff4*switchdist2            &
                               -(1.0_wp/6.0e1_wp)*switchdist6)                 &
                   +shift_b * ( (1.0_wp/7.0_wp)*cutoff7                        &
                               -(1.0_wp/2.0_wp)*cutoff6*switchdist             &
                               +(3.0_wp/5.0_wp)*cutoff5*switchdist2            &
                               -(1.0_wp/4.0_wp)*cutoff4*switchdist3            &
                               +(1.0_wp/1.4e2_wp)*switchdist7)
          enefunc%dispersion_energy = factor*(eswitch-vlong)
          enefunc%dispersion_virial = -2.0_wp*factor*(-vswitch+vlong)

        else

          eswitch = enefunc%eswitch
          vswitch = enefunc%vswitch
          enefunc%dispersion_energy = factor*(eswitch-vlong)
          enefunc%dispersion_virial = -factor*(vswitch+vlong)

        end if

      else

        enefunc%dispersion_energy = factor*(eswitch-vlong)
        enefunc%dispersion_virial = -2.0_wp*factor*(-vswitch+vlong)

      end if
    end if

  end subroutine setup_enefunc_dispcorr

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_fit_refcoord
  !> @brief        setup enefunc fit reference coords
  !! @authors      KT, CK
  !! @param[in]    fitting_file : fitfile or reffile
  !! @param[in]    molecule     : molecule information
  !! @param[inout] enefunc      : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_fit_refcoord(fitting_file, molecule, enefunc)

    ! formal arguments
    integer,                 intent(in)    :: fitting_file
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: j


    if (fitting_file .eq. FittingFileFIT) then
      if (size(molecule%atom_fitcoord(1,:)) /= molecule%num_atoms) &
      call error_msg('Setup_Enefunc_Fit_Refcoord> bad fitfile in [INPUT]')

    else if (fitting_file .eq. FittingFileREF) then
      if (size(enefunc%restraint_refcoord(1,:)) /= molecule%num_atoms) &
      call error_msg('Setup_Enefunc_Fit_Refcoord> bad reffile in [INPUT]')
     
    endif

    enefunc%do_fitting = .true.

    call alloc_enefunc(enefunc, EneFuncFitc, molecule%num_atoms,  &
                       enefunc%num_fitting)

    if (fitting_file .eq. FittingFileFIT) then

      do j = 1, molecule%num_atoms
        enefunc%fit_refcoord(1,j) = molecule%atom_fitcoord(1,j)
        enefunc%fit_refcoord(2,j) = molecule%atom_fitcoord(2,j)
        enefunc%fit_refcoord(3,j) = molecule%atom_fitcoord(3,j)
      end do

    else

      do j = 1, molecule%num_atoms
        enefunc%fit_refcoord(1,j) = enefunc%restraint_refcoord(1,j)
        enefunc%fit_refcoord(2,j) = enefunc%restraint_refcoord(2,j)
        enefunc%fit_refcoord(3,j) = enefunc%restraint_refcoord(3,j)
      end do

    endif

    do j = 1, molecule%num_atoms
      if (enefunc%mass_weight) then
        enefunc%fit_mass(j) = molecule%mass(j)
      else
        enefunc%fit_mass(j) = 1.0_wp
      endif
    end do

    return

  end subroutine setup_enefunc_fit_refcoord

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    remove_fixatm
  !> @brief        remove potential functions for fixed atoms
  !! @authors      KY
  !! @param[in]    boundary : boundary conditions information
  !! @param[inout] enefunc  : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine remove_fixatm(boundary, enefunc)

    ! formal arguments
    type(s_boundary),        intent(in)    :: boundary
    type(s_enefunc),         intent(inout) :: enefunc

    ! bond
    !
    if (enefunc%num_bonds > 0) &
      call remove_fixatm_bond(boundary, enefunc)

    ! angle
    !
    if (enefunc%num_angles > 0) &
      call remove_fixatm_angl(boundary, enefunc)

    ! Urey-Bradley
    !
    if (enefunc%num_ureys > 0) &
      call remove_fixatm_urey(boundary, enefunc)

    ! dihedral
    !
    if (enefunc%num_dihedrals > 0) &
      call remove_fixatm_dihe(boundary, enefunc)

    ! RB dihedral
    !
    if (enefunc%num_rb_dihedrals > 0) &
      call remove_fixatm_rb_dihe(boundary, enefunc)

    ! improper
    !
    if (enefunc%num_impropers > 0) &
      call remove_fixatm_impr(boundary, enefunc)

    ! cmap
    !
    if (enefunc%num_cmaps > 0) &
      call remove_fixatm_cmap(boundary, enefunc)

    if (main_rank) then

      write(MsgOut,'(A)') &
           'Remove_Fixatm> Modified Interactions'
      write(MsgOut,'(A20,I10,A20,I10)')                         &
           '  bond_ene        = ', enefunc%num_bonds,           &
           '  angle_ene       = ', enefunc%num_angles
      if (enefunc%num_ureys > 0)   &
        write(MsgOut,'(A20,I10)')  &
           '  urey_ene        = ', enefunc%num_ureys
      write(MsgOut,'(A20,I10,A20,I10)')                         &
           '  torsion_ene     = ', enefunc%num_dihedrals,       &
           '  improper_ene    = ', enefunc%num_impropers
      if (enefunc%num_rb_dihedrals > 0) &
        write(MsgOut,'(A20,I10)')       &
           '  rb_dihed_ene    = ', enefunc%num_rb_dihedrals
      !write(MsgOut,'(A20,I10)')                                 &
      !     '  cmap_ene        = ', enefunc%num_cmaps
      write(MsgOut,'(A)') ' '
    end if

  end subroutine remove_fixatm

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    remove_fixatm_bond
  !> @brief        remove bond functions for fixed atoms
  !! @authors      KY
  !! @param[in]    boundary : boundary conditions information
  !! @param[inout] enefunc  : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine remove_fixatm_bond(boundary, enefunc)

    ! formal arguments
    type(s_boundary),        intent(in)    :: boundary
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer    :: i, j, k
    integer    :: nbonds, nbonds_new
    integer    :: found
    integer    :: istart, iend
    integer    :: alloc_stat, dealloc_stat

    real(wp),  allocatable   :: temp_fc(:), temp_r0(:)
    integer,   allocatable   :: temp_list(:,:)

    nbonds = enefunc%num_bonds

    ! count the number of bonds including fixed atoms
    !
    found = 0
    do i = 1, nbonds
      if (boundary%fixatm(enefunc%bond_list(1,i)) .and. &
          boundary%fixatm(enefunc%bond_list(2,i)))      &
        found = found + 1
    end do

    if (found == 0) return
    nbonds_new = nbonds - found

    ! copy information from enefunc to temp array
    !
    allocate(temp_fc(nbonds), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    allocate(temp_r0(nbonds), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    allocate(temp_list(2,nbonds), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc

    temp_fc   = enefunc%bond_force_const
    temp_r0   = enefunc%bond_dist_min
    temp_list = enefunc%bond_list

    ! re-allocate EneFuncBond
    !
    call dealloc_enefunc(enefunc, EneFuncBond)

    enefunc%num_bonds = nbonds_new
    call alloc_enefunc(enefunc, EneFuncBond, nbonds_new)

    ! modify enefunc%bond
    !
    found = 0
    do i = 1, nbonds
      if (.not. (boundary%fixatm(temp_list(1,i)) .and. &
                 boundary%fixatm(temp_list(2,i)))) then
        found = found + 1
        enefunc%bond_force_const(found) = temp_fc(i)
        enefunc%bond_dist_min(found)    = temp_r0(i)
        enefunc%bond_list(:,found)      = temp_list(:,i)
      end if

    end do

    ! re-setup loop index
    !
    call get_loop_index(enefunc%num_bonds, istart, iend)
    enefunc%istart_bond = istart
    enefunc%iend_bond   = iend

    ! deallocate temp array
    !
    deallocate(temp_fc, stat=dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc
    deallocate(temp_r0, stat=dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc
    deallocate(temp_list, stat=dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc

  end subroutine remove_fixatm_bond

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    remove_fixatm_angl
  !> @brief        remove angle and UB functions for fixed atoms
  !! @authors      KY
  !! @param[in]    boundary : boundary conditions information
  !! @param[inout] enefunc  : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine remove_fixatm_angl(boundary, enefunc)

    ! formal arguments
    type(s_boundary),        intent(in)    :: boundary
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer    :: i, j, k
    integer    :: nangles, nangles_new
    integer    :: found
    integer    :: istart, iend
    integer    :: alloc_stat, dealloc_stat

    real(wp),  allocatable   :: temp_fc(:), temp_r0(:)
    integer,   allocatable   :: temp_list(:,:)

    nangles = enefunc%num_angles

    ! count the number of angles including fixed atoms
    !
    found = 0
    do i = 1, nangles
      if (boundary%fixatm(enefunc%angl_list(1,i)) .and. &
          boundary%fixatm(enefunc%angl_list(2,i)) .and. &
          boundary%fixatm(enefunc%angl_list(3,i)))      &
        found = found + 1
    end do

    if (found == 0) return

    nangles_new = nangles - found

    ! copy information from enefunc to temp array
    !
    allocate(temp_fc(nangles), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    allocate(temp_r0(nangles), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    allocate(temp_list(3,nangles), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc

    temp_fc   = enefunc%angl_force_const
    temp_r0   = enefunc%angl_theta_min
    temp_list = enefunc%angl_list

    ! re-allocate EneFuncAngl
    !
    call dealloc_enefunc(enefunc, EneFuncAngl)

    enefunc%num_angles = nangles_new
    call alloc_enefunc(enefunc, EneFuncAngl, nangles_new)

    ! modify enefunc%angle
    !
    found = 0
    do i = 1, nangles
      if (.not. (boundary%fixatm(temp_list(1,i)) .and. &
                 boundary%fixatm(temp_list(2,i)) .and. &
                 boundary%fixatm(temp_list(3,i)))) then
        found = found + 1
        enefunc%angl_force_const(found) = temp_fc(i)
        enefunc%angl_theta_min(found)   = temp_r0(i)
        enefunc%angl_list(:,found)      = temp_list(:,i)
      end if

    end do

    ! re-setup loop index
    !
    call get_loop_index(enefunc%num_angles, istart, iend)
    enefunc%istart_angle = istart
    enefunc%iend_angle   = iend

    ! deallocate temp array
    !
    deallocate(temp_fc, stat=dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc
    deallocate(temp_r0, stat=dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc
    deallocate(temp_list, stat=dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc

  end subroutine remove_fixatm_angl

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    remove_fixatm_urey
  !> @brief        remove Urey-Bradley functions for fixed atoms
  !! @authors      KY
  !! @param[in]    boundary : boundary conditions information
  !! @param[inout] enefunc  : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine remove_fixatm_urey(boundary, enefunc)

    ! formal arguments
    type(s_boundary),        intent(in)    :: boundary
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer    :: i, j, k
    integer    :: nureys, nureys_new
    integer    :: found
    integer    :: istart, iend
    integer    :: alloc_stat, dealloc_stat

    real(wp),  allocatable   :: temp_fc(:), temp_r0(:)
    integer,   allocatable   :: temp_list(:,:)

    nureys = enefunc%num_ureys

    ! count the number of ureys including fixed atoms
    !
    found = 0
    do i = 1, nureys
      if (boundary%fixatm(enefunc%urey_list(1,i)) .and. &
          boundary%fixatm(enefunc%urey_list(2,i)))      &
        found = found + 1
    end do

    if (found == 0) return

    nureys_new = nureys - found

    ! copy information from enefunc to temp array
    !
    allocate(temp_fc(nureys), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    allocate(temp_r0(nureys), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    allocate(temp_list(3,nureys), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc

    temp_fc   = enefunc%urey_force_const
    temp_r0   = enefunc%urey_rmin
    temp_list = enefunc%urey_list

    ! re-allocate EneFuncUrey
    !
    call dealloc_enefunc(enefunc, EneFuncUrey)

    enefunc%num_ureys = nureys_new
    call alloc_enefunc(enefunc, EneFuncUrey, nureys_new)

    ! modify enefunc%urey
    !
    found = 0
    do i = 1, nureys
      if (.not. (boundary%fixatm(temp_list(1,i)) .and. &
                 boundary%fixatm(temp_list(2,i)))) then
        found = found + 1
        enefunc%urey_force_const(found) = temp_fc(i)
        enefunc%urey_rmin(found)        = temp_r0(i)
        enefunc%urey_list(:,found)      = temp_list(:,i)
      end if

    end do

    ! re-setup loop index
    !
    call get_loop_index(enefunc%num_ureys, istart, iend)
    enefunc%istart_urey = istart
    enefunc%iend_urey   = iend

    ! deallocate temp array
    !
    deallocate(temp_fc, stat=dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc
    deallocate(temp_r0, stat=dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc
    deallocate(temp_list, stat=dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc

  end subroutine remove_fixatm_urey

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    remove_fixatm_dihe
  !> @brief        remove dihedral functions for fixed atoms
  !! @authors      KY
  !! @param[in]    boundary : boundary conditions information
  !! @param[inout] enefunc  : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine remove_fixatm_dihe(boundary, enefunc)

    ! formal arguments
    type(s_boundary),        intent(in)    :: boundary
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer    :: i, j, k
    integer    :: ndihedrals, ndihedrals_new
    integer    :: found
    integer    :: istart, iend
    integer    :: alloc_stat, dealloc_stat

    real(wp),  allocatable   :: temp_fc(:), temp_phase(:)
    real(wp),  allocatable   :: temp_scee(:), temp_scnb(:)
    integer,   allocatable   :: temp_list(:,:), temp_period(:)

    ndihedrals = enefunc%num_dihedrals

    ! count the number of dihedrals including fixed atoms
    !
    found = 0
    do i = 1, ndihedrals
      if (boundary%fixatm(enefunc%dihe_list(1,i)) .and. &
          boundary%fixatm(enefunc%dihe_list(2,i)) .and. &
          boundary%fixatm(enefunc%dihe_list(3,i)) .and. &
          boundary%fixatm(enefunc%dihe_list(4,i)))      &
        found = found + 1
    end do

    if (found == 0) return

    ndihedrals_new = ndihedrals - found

    ! copy information from enefunc to temp array
    !
    allocate(temp_fc(ndihedrals), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    allocate(temp_phase(ndihedrals), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    allocate(temp_scee(ndihedrals), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    allocate(temp_scnb(ndihedrals), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    allocate(temp_period(ndihedrals), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    allocate(temp_list(4,ndihedrals), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc

    temp_fc     = enefunc%dihe_force_const
    temp_phase  = enefunc%dihe_phase
    temp_scee   = enefunc%dihe_scee
    temp_scnb   = enefunc%dihe_scnb
    temp_period = enefunc%dihe_periodicity
    temp_list   = enefunc%dihe_list

    ! re-allocate EneFuncDihe
    !
    call dealloc_enefunc(enefunc, EneFuncDihe)

    enefunc%num_dihedrals = ndihedrals_new
    call alloc_enefunc(enefunc, EneFuncDihe, ndihedrals_new)

    ! modify enefunc%dihe
    !
    found = 0
    do i = 1, ndihedrals
      if (.not. (boundary%fixatm(temp_list(1,i)) .and. &
                 boundary%fixatm(temp_list(2,i)) .and. &
                 boundary%fixatm(temp_list(3,i)) .and. &
                 boundary%fixatm(temp_list(4,i)))) then
        found = found + 1
        enefunc%dihe_force_const(found) = temp_fc(i)
        enefunc%dihe_phase(found)       = temp_phase(i)
        enefunc%dihe_scee(found)        = temp_scee(i)
        enefunc%dihe_scnb(found)        = temp_scnb(i)
        enefunc%dihe_periodicity(found) = temp_period(i)
        enefunc%dihe_list(:,found)      = temp_list(:,i)
      end if

    end do

    ! re-setup loop index
    !
    call get_loop_index(enefunc%num_dihedrals, istart, iend)
    enefunc%istart_dihedral = istart
    enefunc%iend_dihedral   = iend

    ! deallocate temp array
    !
    deallocate(temp_fc, stat=dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc
    deallocate(temp_phase, stat=dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc
    deallocate(temp_scee, stat=dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc
    deallocate(temp_scnb, stat=dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc
    deallocate(temp_period, stat=dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc
    deallocate(temp_list, stat=dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc

  end subroutine remove_fixatm_dihe

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    remove_fixatm_rb_dihe
  !> @brief        remove RB dihedral functions for fixed atoms
  !! @authors      KY
  !! @param[in]    boundary : boundary conditions information
  !! @param[inout] enefunc  : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine remove_fixatm_rb_dihe(boundary, enefunc)

    ! formal arguments
    type(s_boundary),        intent(in)    :: boundary
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer    :: i, j, k
    integer    :: ndihedrals, ndihedrals_new
    integer    :: found
    integer    :: istart, iend
    integer    :: alloc_stat, dealloc_stat

    real(wp),  allocatable   :: temp_c(:,:)
    integer,   allocatable   :: temp_list(:,:)

    ndihedrals = enefunc%num_rb_dihedrals

    ! count the number of RB dihedrals including fixed atoms
    !
    found = 0
    do i = 1, ndihedrals
      if (boundary%fixatm(enefunc%rb_dihe_list(1,i)) .and. &
          boundary%fixatm(enefunc%rb_dihe_list(2,i)) .and. &
          boundary%fixatm(enefunc%rb_dihe_list(3,i)) .and. &
          boundary%fixatm(enefunc%rb_dihe_list(4,i)))      &
        found = found + 1
    end do

    if (found == 0) return

    ndihedrals_new = ndihedrals - found

    ! copy information from enefunc to temp array
    !
    allocate(temp_c(6,ndihedrals), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    allocate(temp_list(4,ndihedrals), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc

    temp_c    = enefunc%rb_dihe_c
    temp_list = enefunc%rb_dihe_list

    ! re-allocate EneFuncRBDihe
    !
    call dealloc_enefunc(enefunc, EneFuncRBDihe)

    enefunc%num_rb_dihedrals = ndihedrals_new
    call alloc_enefunc(enefunc, EneFuncRBDihe, ndihedrals_new)

    ! modify RB enefunc%dihe
    !
    found = 0
    do i = 1, ndihedrals
      if (.not. (boundary%fixatm(temp_list(1,i)) .and. &
                 boundary%fixatm(temp_list(2,i)) .and. &
                 boundary%fixatm(temp_list(3,i)) .and. &
                 boundary%fixatm(temp_list(4,i)))) then
        found = found + 1
        enefunc%rb_dihe_c(:,found)    = temp_c(:,i)
        enefunc%rb_dihe_list(:,found) = temp_list(:,i)
      end if

    end do

    ! re-setup loop index
    !
    call get_loop_index(enefunc%num_rb_dihedrals, istart, iend)
    enefunc%istart_rb_dihed = istart
    enefunc%iend_rb_dihed   = iend

    ! deallocate temp array
    !
    deallocate(temp_c, stat=dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc
    deallocate(temp_list, stat=dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc

  end subroutine remove_fixatm_rb_dihe

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    remove_fixatm_impr
  !> @brief        remove improper functions for fixed atoms
  !! @authors      KY
  !! @param[in]    boundary : boundary conditions information
  !! @param[inout] enefunc  : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine remove_fixatm_impr(boundary, enefunc)

    ! formal arguments
    type(s_boundary),        intent(in)    :: boundary
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer    :: i, j, k
    integer    :: nimpropers, nimpropers_new
    integer    :: found
    integer    :: istart, iend
    integer    :: alloc_stat, dealloc_stat

    real(wp),  allocatable   :: temp_fc(:), temp_phase(:)
    integer,   allocatable   :: temp_list(:,:), temp_period(:)

    nimpropers = enefunc%num_impropers

    ! count the number of impropers including fixed atoms
    !
    found = 0
    do i = 1, nimpropers
      if (boundary%fixatm(enefunc%impr_list(1,i)) .and. &
          boundary%fixatm(enefunc%impr_list(2,i)) .and. &
          boundary%fixatm(enefunc%impr_list(3,i)) .and. &
          boundary%fixatm(enefunc%impr_list(4,i)))      &
        found = found + 1
    end do

    if (found == 0) return

    nimpropers_new = nimpropers - found

    ! copy information from enefunc to temp array
    !
    allocate(temp_fc(nimpropers), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    allocate(temp_period(nimpropers), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    allocate(temp_phase(nimpropers), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    allocate(temp_list(4,nimpropers), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc

    temp_fc     = enefunc%impr_force_const
    temp_period = enefunc%impr_periodicity
    temp_phase  = enefunc%impr_phase
    temp_list   = enefunc%impr_list

    ! re-allocate EneFuncDihe
    !
    call dealloc_enefunc(enefunc, EneFuncImpr)

    enefunc%num_impropers = nimpropers_new
    call alloc_enefunc(enefunc, EneFuncImpr, nimpropers_new)

    ! modify enefunc%dihe
    !
    found = 0
    do i = 1, nimpropers
      if (.not. (boundary%fixatm(temp_list(1,i)) .and. &
                 boundary%fixatm(temp_list(2,i)) .and. &
                 boundary%fixatm(temp_list(3,i)) .and. &
                 boundary%fixatm(temp_list(4,i)))) then
        found = found + 1
        enefunc%impr_force_const(found) = temp_fc(i)
        enefunc%impr_periodicity(found) = temp_period(i)
        enefunc%impr_phase(found)       = temp_phase(i)
        enefunc%impr_list(:,found)      = temp_list(:,i)
      end if

    end do

    ! re-setup loop index
    !
    call get_loop_index(enefunc%num_impropers, istart, iend)
    enefunc%istart_improper = istart
    enefunc%iend_improper   = iend

    ! deallocate temp array
    !
    deallocate(temp_fc, stat=dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc
    deallocate(temp_period, stat=dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc
    deallocate(temp_phase, stat=dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc
    deallocate(temp_list, stat=dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc

  end subroutine remove_fixatm_impr

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    remove_fixatm_cmap
  !> @brief        remove CMAP functions for fixed atoms
  !! @authors      KY
  !! @param[in]    boundary : boundary conditions information
  !! @param[inout] enefunc  : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine remove_fixatm_cmap(boundary, enefunc)

    ! formal arguments
    type(s_boundary),        intent(in)    :: boundary
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer    :: i, j, k
    integer    :: ncmaps, ncmaps_new
    integer    :: found
    integer    :: istart, iend
    integer    :: alloc_stat, dealloc_stat

    real(wp),  allocatable   :: temp_force(:,:,:), temp_coef(:,:,:,:,:)
    integer,   allocatable   :: temp_list(:,:), temp_res(:), temp_type(:)

    ! TODO

  end subroutine remove_fixatm_cmap

end module at_enefunc_mod
