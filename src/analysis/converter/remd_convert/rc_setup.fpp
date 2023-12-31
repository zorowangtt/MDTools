!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   rc_setup_mod
!> @brief   setup variables and structures in REMD convert
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module rc_setup_mod

  use rc_control_mod
  use rc_option_mod
  use rc_option_str_mod
  use fitting_mod
  use fitting_str_mod
  use input_mod
  use input_str_mod
  use output_mod
  use output_str_mod
  use trajectory_str_mod
  use select_mod
  use molecules_mod
  use molecules_str_mod
  use fileio_grocrd_mod
  use fileio_grotop_mod
  use fileio_ambcrd_mod
  use fileio_prmtop_mod
  use fileio_psf_mod
  use fileio_pdb_mod
  use constants_mod
 
  implicit none
  private

  ! subroutines
  public :: setup

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup
  !> @brief        setup variables and structures in CRD convert
  !! @authors      NT
  !! @param[in]    ctrl_data  : information of control parameters
  !! @param[inout] molecule   : molecule information
  !! @param[inout] trajectory : trajectory information
  !! @param[inout] fitting    : fitting information
  !! @param[inout] option     : option information
  !! @param[inout] input      : input information
  !! @param[inout] output     : output information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup(ctrl_data,  &
                   molecule,   &
                   trajectory, &
                   fitting,    &
                   option,     &
                   input,      &
                   output)

    ! formal arguments
    type(s_ctrl_data),       intent(in)    :: ctrl_data
    type(s_molecule),        intent(inout) :: molecule
    type(s_trajectory),      intent(inout) :: trajectory
    type(s_fitting),         intent(inout) :: fitting
    type(s_option),          intent(inout) :: option
    type(s_input),           intent(inout) :: input
    type(s_output),          intent(inout) :: output

    ! local variables
    type(s_psf)              :: psf
    type(s_pdb)              :: ref, ref_out
    type(s_prmtop)           :: prmtop
    type(s_ambcrd)           :: ambcrd
    type(s_grotop)           :: grotop
    type(s_grocrd)           :: grocrd


    ! input files
    !
    call input_files(ctrl_data%inp_info, &
                     psf=psf,            &
                     ref=ref,            &
                     prmtop=prmtop,      &
                     ambcrd=ambcrd,      &
                     grotop=grotop,      &
                     grocrd=grocrd)

    ! define molecules
    !
    call define_molecules(molecule, pdb=ref,       &
                                    psf=psf,       &
                                    prmtop=prmtop, &
                                    ambcrd=ambcrd, &
                                    grotop=grotop, &
                                    grocrd=grocrd)

    call dealloc_psf_all(psf)
    call dealloc_pdb_all(ref)
    call dealloc_prmtop_all(prmtop)
    call dealloc_ambcrd_all(ambcrd)
    call dealloc_grotop_all(grotop)
    call dealloc_grocrd_all(grocrd)


    call alloc_trajectory(trajectory, molecule%num_atoms)
    trajectory%coord(:,:)   = 0.0_wp
    trajectory%pbc_box(:,:) = 0.0_wp


    ! setup input
    !
    call setup_input(ctrl_data%inp_info, input)


    ! setup selection
    !
    call setup_selection(ctrl_data%sel_info, molecule)


    ! setup fitting
    !
    call setup_fitting(ctrl_data%fit_info, ctrl_data%sel_info, &
                       molecule, fitting)


    ! setup option
    !
    call setup_option(ctrl_data%opt_info, ctrl_data%sel_info, &
                      molecule, option)


    ! setup output 
    !
    call setup_output(ctrl_data%out_info, output)


    ! export reference molecules
    !
    if (output%pdbfile /= '') then

      call export_molecules(molecule, option%trjout_atom, ref_out)
      call output_pdb(output%pdbfile, ref_out)

      call dealloc_pdb_all(ref_out)

    end if

    return

  end subroutine setup

end module rc_setup_mod
