!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_vibration_mod
!> @brief   perform vibrational analysis
!! @authors Kiyoshi Yagi (KY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_vibration_mod

  use at_dynvars_str_mod
  use at_enefunc_str_mod
  use at_boundary_str_mod
  use at_output_str_mod
  use at_pairlist_str_mod
  use at_vibration_str_mod
  use at_energy_mod
  use at_output_mod
  use at_qmmm_mod
  use molecules_str_mod
  use fileio_mod
  use fileio_control_mod
  use select_mod
  use select_atoms_mod
  use select_atoms_str_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef MPI
  use mpi
#endif

  implicit none
  private

  ! structures
  type, public :: s_vib_info

    ! For vibrational analysis
    integer          :: runmode             = RunModeHarm
    integer          :: nreplica            = 1
    character(256)   :: vibatm_select_index = ''
    character(256)   :: output_minfo_atm    = ''
    character(256)   :: minfo_folder        = 'minfo.files'
    real(wp)         :: diff_stepsize       = 0.001_wp
    character(256)   :: gridfile            = ''
    character(256)   :: datafile            = ''
!   obsolent options
    character(256)   :: gridfolderID        = 'a'
    logical          :: dryrun              = .false.
!    real(wp)         :: cutoff              = -1.0

  end type s_vib_info

  ! subroutines
  public  :: show_ctrl_vibration
  public  :: read_ctrl_vibration
  public  :: setup_vibration
  private :: setup_vibatoms
  public  :: run_vib
  private :: harmonic_analysis
  private :: diagonalize
  private :: print_minfo
  private :: print_minfo_grad
  private :: calc_gridpoints

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_vibration
  !> @brief        show usage of VIBRATION section
  !! @authors      KY
  !! @param[in]    show_all : show all usage or not
  !! @param[in]    run_mode : run mode string : "md", "min", "remd", "rpath"
  !!                          "vib"
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine show_ctrl_vibration(show_all, run_mode)

    ! formal arguments
    logical,                 intent(in)    :: show_all
    character(*),            intent(in)    :: run_mode


    if (show_all) then

      select case (run_mode)

      case ('vib')

        write(MsgOut,'(A)') '[SELECTION]'
        write(MsgOut,'(A)') 'group1              = atomno:5-8  # Vib group 1'
        write(MsgOut,'(A)') 'group2              = sid:PROA    # Vib group 2'
        write(MsgOut,'(A)') ' '

        write(MsgOut,'(A)') '[VIBRATION]'
        write(MsgOut,'(A)') 'nreplica            = 1           # number of MPI processes' 
        write(MsgOut,'(A)') 'runmode             = HARM        # HARM, QFF, or GRID'
        write(MsgOut,'(A)') 'vibatm_select_index = 1           # atoms subject to vib analysis'
        write(MsgOut,'(A)') 'output_minfo_atm    = 2           # atoms punched to a minfo file'
        write(MsgOut,'(A)') 'minfo_folder        = minfo.files # a folder where minfo files are created'
        write(MsgOut,'(A)') '# diff_stepsize     = 0.001       # displacement for numerical diff.'
        write(MsgOut,'(A)') '# gridfile          = grid.xyz    # the xyz file containing coordinates of grid points'
        write(MsgOut,'(A)') '# datafile          = grid.dat    # the output file for grid data'
        ! obsolent
        !write(MsgOut,'(A)') '# gridfolderID      = a           # ID of a folder to store the grid point data'
        !write(MsgOut,'(A)') '# dryrun            = false     # if true, only generate input files and exit'
        !write(MsgOut,'(A)') 'cutoff              = -1.0      # cutoff distance (in Angs) for Hessian evaluation.'
        write(MsgOut,'(A)') ' '

      end select

    else

      select case (run_mode)

      case ('vib')

        write(MsgOut,'(A)') '[SELECTION]'
        write(MsgOut,'(A)') 'group1              = atomno:5-8  # Vib group 1'
        write(MsgOut,'(A)') 'group2              = sid:PROA    # Vib group 2'
        write(MsgOut,'(A)') ' '

        write(MsgOut,'(A)') '[VIBRATION]'
        write(MsgOut,'(A)') 'nreplica            = 1           # number of MPI processes' 
        write(MsgOut,'(A)') 'runmode             = HARM        # HARM, QFF, or GRID'
        write(MsgOut,'(A)') 'vibatm_select_index = 1           # fixed atoms in minimization'
        write(MsgOut,'(A)') 'output_minfo_atm    = 2           # atoms punched to a minfo file'
        write(MsgOut,'(A)') 'minfo_folder        = minfo.files # a folder where minfo files are created'
        write(MsgOut,'(A)') ' '

      end select


    end if

    return

  end subroutine show_ctrl_vibration
  
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_vibration
  !> @brief        read VIBRATIOn section in the control file
  !! @authors      KY
  !! @param[in]    handle   : unit number of control file
  !! @param[out]   vib_info : MINIMIZE section in control parameters
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ctrl_vibration(handle, vib_info)

    ! parameters
    character(*),            parameter     :: Section = 'Vibration'

    ! formal arguments
    integer,                 intent(in)    :: handle
    type(s_vib_info),        intent(inout) :: vib_info

    ! local
    logical :: found_error


    ! read parameters from control file
    !
    call begin_ctrlfile_section(handle, Section)

    ! Vibrational analysis
    call read_ctrlfile_type   (handle, Section, 'runmode',              &
                               vib_info%runmode, RunModeTypes)
    call read_ctrlfile_integer(handle, Section, 'nreplica',             &
                               vib_info%nreplica)
    call read_ctrlfile_string (handle, Section, 'vibatm_select_index',  &
                               vib_info%vibatm_select_index)
    call read_ctrlfile_string (handle, Section, 'output_minfo_atm',     &
                               vib_info%output_minfo_atm)
    call read_ctrlfile_real   (handle, Section, 'diff_stepsize',        &
                               vib_info%diff_stepsize)
    call read_ctrlfile_string (handle, Section, 'minfo_folder',          &
                               vib_info%minfo_folder)
!    call read_ctrlfile_real   (handle, Section, 'cutoff',               &
!                               vib_info%cutoff)
    call read_ctrlfile_string (handle, Section, 'gridfile',             &
                               vib_info%gridfile)
    call read_ctrlfile_string (handle, Section, 'datafile',             &
                               vib_info%datafile)
!    call read_ctrlfile_string (handle, Section, 'gridfolderID',         &
!                               vib_info%gridfolderID)
!    call read_ctrlfile_logical(handle, Section, 'dryrun',               &
!                               vib_info%dryrun)

    call end_ctrlfile_section(handle)

    if (trim(vib_info%vibatm_select_index) .eq. '') &
      call error_msg('Read_Ctrl_Vibration> No VIB atoms defined')

    select case (vib_info%runmode)
    case (RunModeQFF)

      ! Generate QFF
      if (vib_info%gridfile == '') vib_info%gridfile = 'makeQFF.xyz'

    case (RunModeGRID)
      ! Generate Grid-PES
      if (vib_info%gridfile == '') vib_info%gridfile = 'makeGrid.xyz'
      if (vib_info%datafile == '') vib_info%datafile = 'makeGrid.dat'

    end select

    ! write parameters to MsgOut
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Read_Ctrl_Vibration> Parameters of VIBRATION'
      write(MsgOut,'(A,A10)')    &
          '  runmode             = ', trim(RunModeTypes(vib_info%runmode))
      write(MsgOut,'(A,I10)')    &
          '  nreplica            = ', vib_info%nreplica
      write(MsgOut,'(A,A)')      &
          '  vibatm_select_index = ', trim(vib_info%vibatm_select_index)
      if (vib_info%output_minfo_atm /= '') write(MsgOut,'(A,A)')      &
          '  output_minfo_atm    = ', trim(vib_info%output_minfo_atm)
      write(MsgOut,'(A,A)')      &
          '  minfo_folder        = ', trim(vib_info%minfo_folder)

      select case (vib_info%runmode)

      case (RunModeHARM)

        ! Harmonic vibrational analysis
        write(MsgOut,'(A,E10.2)')&
          '  diff_stepsize       = ', vib_info%diff_stepsize
        ! write(MsgOut,'(A,F10.5)')       &
        !   '  cutoff              = ', vib_info%cutoff

      case (RunModeQFF)

        write(MsgOut,'(A,A)')    &
          '  gridfile            = ', trim(vib_info%gridfile)
!        write(MsgOut,'(A,A10)')  &
!          '  gridfolderID        = ', trim(vib_info%gridfolderID)

      case (RunModeGRID)

        write(MsgOut,'(A,A)')    &
          '  gridfile            = ', trim(vib_info%gridfile)
        write(MsgOut,'(A,A)')    &
          '  datafile            = ', trim(vib_info%datafile)
!        write(MsgOut,'(A,A10)')  &
!          '  gridfolderID        = ', trim(vib_info%gridfolderID)

      end select

!      if (vib_info%dryrun) then
!        write(MsgOut, '(A)')     &
!          '  dryrun              =        yes'
!      end if

      write(MsgOut,'(A)') ' '

    end if

    return

  end subroutine read_ctrl_vibration

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_vibration
  !> @brief        setup vibration information
  !> @authors      KY
  !! @param[in]    vib_info  : VIBRATION section in control parameters
  !! @param[in]    sel_info  : SELECTION section in control parameters
  !! @param[in]    molecule  : molecular information
  !! @param[in]    qmmm      : QM/MM information
  !! @param[out]   vibration : vibration information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_vibration(vib_info, sel_info, molecule, qmmm, vibration)

    ! formal arguments
    type(s_vib_info),        intent(in)    :: vib_info
    type(s_sel_info),        intent(in)    :: sel_info
    type(s_molecule),        intent(in)    :: molecule
    type(s_qmmm),            intent(inout) :: qmmm
    type(s_vibration),       intent(inout) :: vibration


    call setup_vibatoms(vib_info, sel_info, molecule, vibration)

    vibration%nreplica         = vib_info%nreplica
    vibration%diff_stepsize    = vib_info%diff_stepsize
    vibration%minfo_folder     = vib_info%minfo_folder

    select case (vib_info%runmode)
    case(RunModeHARM)
      vibration%gengrid          = .false.
      vibration%grid_ene_only    = .false.
      !vibration%cutoff           = vib_info%cutoff

    case(RunModeQFF )
      vibration%gengrid          = .true.
      vibration%grid_ene_only    = .false.

    case(RunModeGRID)
      vibration%gengrid          = .true.
      vibration%grid_ene_only    = .true.

    end select

    vibration%gridfile         = vib_info%gridfile
    vibration%datafile         = vib_info%datafile
    vibration%gridfolderID     = vib_info%gridfolderID
    vibration%dryrun           = vib_info%dryrun
    if(qmmm%do_qmmm .and. vib_info%dryrun) then
      qmmm%dryrun = .true.
    endif

  end subroutine setup_vibration

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_vibatoms
  !> @brief        define atoms subject to vibrational analysis
  !! @authors      KY
  !! @param[in]    vib_info    : vibration input information
  !! @param[in]    sel_info    : selector input information
  !! @param[in]    molecule    : molecular information
  !! @param[out]   vibration    : vibration parameters
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_vibatoms(vib_info, sel_info, molecule, vibration)

    ! formal arguments
    type(s_vib_info),        intent(in) :: vib_info
    type(s_sel_info),        intent(in) :: sel_info
    type(s_molecule),        intent(in) :: molecule
    type(s_vibration),    intent(inout) :: vibration

    ! local variables
    integer                :: igroup, ng_vib, ng_minfo, ng_minfo_bk, natom, &
                              i, j, nn, offset, temp, iatom
    integer, allocatable   :: glist_vib(:), glist_minfo(:), glist_minfo_bk(:)
    type(s_selatoms), allocatable :: selatoms(:)

    integer, parameter :: max_atm_print = 100

    ! Number of atoms for vibrational analysis
    !
    ng_vib = split_num(trim(vib_info%vibatm_select_index))

    allocate(glist_vib(ng_vib))
    call split(ng_vib, ng_vib, vib_info%vibatm_select_index, glist_vib)

    allocate(selatoms(ng_vib))

    natom = 0
    do i = 1, ng_vib
      igroup = glist_vib(i)
      call select_atom(molecule, sel_info%groups(igroup), selatoms(i))
      natom = natom + size(selatoms(i)%idx)
    end do

    vibration%vib_natoms = natom

    ! List of vibanal atoms
    !
    allocate(vibration%vibatom_id(vibration%vib_natoms))

    offset = 0
    do i = 1, ng_vib
      igroup = glist_vib(i)
      natom = size(selatoms(i)%idx)
      vibration%vibatom_id(offset+1:offset+natom) = selatoms(i)%idx(1:natom)
      offset = offset + natom
    end do

    deallocate(selatoms)

    ! sort vibanal atom indices in ascending order
    !
    do i = vibration%vib_natoms, 2, -1
      do j = 1, i - 1
        if (vibration%vibatom_id(j) > vibration%vibatom_id(j+1)) then
          temp = vibration%vibatom_id(j)
          vibration%vibatom_id(j)   = vibration%vibatom_id(j+1)
          vibration%vibatom_id(j+1) = temp
        end if
      end do
    end do

    ! Number of atoms for minfo
    !
    ng_minfo = split_num(trim(vib_info%output_minfo_atm))

    if (ng_minfo > 0) then
      allocate(glist_minfo(ng_minfo))
      call split(ng_minfo, ng_minfo, vib_info%output_minfo_atm, glist_minfo)

      allocate(selatoms(ng_minfo))

      natom = 0
      do i = 1, ng_minfo
        igroup = glist_minfo(i)
        call select_atom(molecule, sel_info%groups(igroup), selatoms(i))
        natom = natom + size(selatoms(i)%idx)
      end do

      vibration%minfo_natoms = natom

      ! List of minfo subatoms
      !
      allocate(vibration%minfoatom_id(vibration%minfo_natoms))

      offset = 0
      do i = 1, ng_minfo
        igroup = glist_minfo(i)
        natom = size(selatoms(i)%idx)
        vibration%minfoatom_id(offset+1:offset+natom) = selatoms(i)%idx(1:natom)
        offset = offset + natom
      end do

      deallocate(selatoms)

      ! sort atom indices in ascending order
      !
      do i = vibration%minfo_natoms, 2, -1
        do j = 1, i - 1
          if (vibration%minfoatom_id(j) > vibration%minfoatom_id(j+1)) then
            temp = vibration%minfoatom_id(j)
            vibration%minfoatom_id(j)   = vibration%minfoatom_id(j+1)
            vibration%minfoatom_id(j+1) = temp
          end if
        end do
      end do

    else
      vibration%minfo_natoms = 0

    end if

    ! Punch out vibatom info.
    !
    if (main_rank) then
      write(MsgOut,'(a)') &
        "Setup_Vibration_Atoms> Atoms subject to vibrational analysis"

      do i = 1, vibration%vib_natoms
        iatom   = vibration%vibatom_id(i)
        write(MsgOut,'(2i6,x,a4,x,i6,x,a4,x,a6,x,a6)') &
          i,                             &
          iatom,                         &
          molecule%segment_name(iatom),  &
          molecule%residue_no(iatom),    &
          molecule%residue_name(iatom),  &
          molecule%atom_name(iatom),     &
          molecule%atom_cls_name(iatom)
      end do

      write(MsgOut,'(a,i0)') "  number of VIB atoms = ", vibration%vib_natoms
      write(MsgOut, '(a)') ' '

      if (vibration%minfo_natoms > 0) then
        write(MsgOut,'(a)') &
          "Setup_Vibration_Atoms> Atoms punched to minfo file in addition"

        if (vibration%minfo_natoms < max_atm_print) then
          do i = 1, vibration%minfo_natoms
            iatom   = vibration%minfoatom_id(i)
            write(MsgOut,'(2i6,x,a4,x,i6,x,a4,x,a6,x,a6)') &
              i,                             &
              iatom,                         &
              molecule%segment_name(iatom),  &
              molecule%residue_no(iatom),    &
              molecule%residue_name(iatom),  &
              molecule%atom_name(iatom),     &
              molecule%atom_cls_name(iatom)
          end do
        end if

        write(MsgOut,'(a,i0)') "  number of atoms = ", vibration%minfo_natoms
        write(MsgOut, '(a)') ' '
      end if
    end if

    deallocate(glist_vib)

    return

  end subroutine setup_vibatoms

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    run_vib
  !> @brief        perform vibrational analysis
  !! @authors      YA, KY
  !! @param[inout] molecule    : molecule information
  !! @param[inout] enefunc     : potential energy functions
  !! @param[inout] dynvars     : dynamic variables
  !! @param[inout] vibration   : vibration information
  !! @param[inout] pairlist    : non-bond pair list
  !! @param[inout] boundary    : boundary condition
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine run_vib(molecule, enefunc, dynvars, vibration, &
                     output, pairlist, boundary)

    ! formal arguments
    type(s_molecule), target, intent(inout) :: molecule
    type(s_enefunc),          intent(inout) :: enefunc
    type(s_dynvars),  target, intent(inout) :: dynvars
    type(s_vibration),        intent(inout) :: vibration
    type(s_output),           intent(inout) :: output
    type(s_pairlist),         intent(inout) :: pairlist
    type(s_boundary),         intent(inout) :: boundary


    ! Open output files
    !
    call open_output(output)

    if(.not. vibration%gengrid) then
      ! Harmonic vibrational analysis
      call harmonic_analysis(molecule, enefunc, dynvars, vibration, &
                        output, pairlist, boundary)
    else
      ! Calc ene/grad at grid points
      call calc_gridpoints(molecule, enefunc, dynvars, vibration, &
                        output, pairlist, boundary)
    endif

    ! close output files
    !
    call close_output(output)

    return

  end subroutine run_vib

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    harmonic_analysis
  !> @brief        harmonic vibrational analysis
  !> @authors      KY
  !! @param[inout] molecule  : molecule information
  !! @param[inout] enefunc   : potential energy functions information
  !! @param[inout] dynvars   : dynamic variables information
  !! @param[inout] vibration : vibration information
  !! @param[inout] output    : output information
  !! @param[inout] pairlist  : non-bond pair list information
  !! @param[inout] boundary  : boundary conditions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine harmonic_analysis(molecule, enefunc, dynvars, vibration, &
                          output, pairlist, boundary)

    ! formal arguments
    type(s_molecule),  target, intent(inout) :: molecule
    type(s_enefunc),           intent(inout) :: enefunc
    type(s_dynvars),   target, intent(inout) :: dynvars
    type(s_vibration), target, intent(inout) :: vibration
    type(s_output),            intent(inout) :: output
    type(s_pairlist),          intent(inout) :: pairlist
    type(s_boundary),          intent(inout) :: boundary

    ! local variables
    integer               :: i, j, k, i1, j1, k1, k2, nat, nat3, iatom
    integer               :: count, icount, ncolomn, nc
    real(wp)              :: energy_ref, energy1, energy2, rmsg, delta, x0
    real(wp)              :: dipole_ref(3), dipole1(3), dipole2(3)
    real(wp), allocatable :: grad1(:), grad2(:), sq_mass3(:),   &
                             hessian(:,:), hess_pack(:), hess_mw_pack(:),  &
                             dipole_derivative(:,:), dd_mw(:,:), dd2(:,:), &
                             infrared(:)
    real(wp), pointer     :: coord(:,:), mass(:)
    integer , pointer     :: vibatom_id(:)
    real(wp), allocatable :: grad_ref(:), vec(:,:), omega(:)
    integer               :: replicaid, ierr

    character(256) :: folder, minfo_folder, basename, fname
    character      :: num*5
    character(1)   :: xyz(3)
    logical        :: savefile, ex

    ! use pointers
    !
    nat        =  vibration%vib_natoms
    vibatom_id => vibration%vibatom_id
    coord      => dynvars%coord
    mass       => molecule%mass

    ! allocate local variables
    nat3 = nat*3
    allocate(grad_ref(nat3), grad1(nat3), grad2(nat3), hessian(nat3,nat3))
    allocate(dipole_derivative(3,nat3))

    ! replica id
    replicaid = my_country_no + 1

    ! for print 
    xyz(1) = 'X'
    xyz(2) = 'Y'
    xyz(3) = 'Z'

    ! Print message
    if(main_rank) then
      write(MsgOut,'(''Enter vibrational analysis'')')
      write(MsgOut,*)

      write(MsgOut,'(''  Cartesian coordinates of atoms for vib. analysis'')')
      do i = 1, nat
        iatom = vibatom_id(i)
        write(MsgOut,'(2x,i4,2x,i9,x,a4,x,i6,x,a4,x,a6,4x,3f18.10)') &
          i, iatom, molecule%segment_name(iatom), molecule%residue_no(iatom), &
          molecule%residue_name(iatom), molecule%atom_name(iatom), &
          coord(1:3,iatom)
      end do
      write(MsgOut,*)
    end if

    if(enefunc%qmmm%do_qmmm) then
      folder   = 'qmmm_vib'
      icount = 0
    end if

    ! create a folder to save minfo files
    minfo_folder = vibration%minfo_folder
    call system('mkdir -p '//trim(minfo_folder)//' >& /dev/null')

    ! Generate Hessian matrix by numerical differentiations of gradients
    !
    if(main_rank) then
      write(MsgOut,'(''  Generate Hessian matrix by num. diff. of gradients  '')')
      write(MsgOut,'(''    Loop over atoms'')')
    end if

    count = 0
    dynvars%step = count

    ! ------------------------------------------------------------
    ! Energy and gradient at the current geometry using replica1

    if(replicaid == 1) then
      write(num,'(i5.5)') count
      basename = 'vib'//num
      fname = trim(minfo_folder)//'/'//trim(basename)

      inquire(file=trim(fname)//'.minfo', exist = ex)
      if (ex) then
        call read_minfo_grad(fname, nat, energy_ref, grad_ref, dipole_ref)
        dynvars%energy%total = energy_ref * CONV_UNIT_ENE
        k1 = 1
        do i1 = 1, nat
          do j1 = 1, 3
            dynvars%force(j1,vibatom_id(i1)) = -grad_ref(k1) * CONV_UNIT_FORCE
            k1 = k1 + 1
          end do
        end do
        enefunc%qmmm%qm_dipole = dipole_ref

      else
        call compute_energy(molecule, enefunc, pairlist, boundary, .true., &
                            .false.,                &
                            dynvars%coord,          &
                            dynvars%energy,         &
                            dynvars%temporary,      &
                            dynvars%force,          &
                            dynvars%virial,         &
                            dynvars%virial_extern)

        if(replica_main_rank) &
          write(MsgOut,'(6x,"Done for    0   input",7x,"replicaID = ",i5, &
                         6x,"energy = ",f20.6)') &
                replicaid, dynvars%energy%total 

        call output_vib(output, molecule, enefunc, vibration, boundary, dynvars)

        energy_ref =  dynvars%energy%total / CONV_UNIT_ENE

        k1 = 1
        do i1 = 1, nat
          do j1 = 1, 3
            grad_ref(k1) = -dynvars%force(j1,vibatom_id(i1)) / CONV_UNIT_FORCE
            k1 = k1 + 1
          end do
        end do
        dipole_ref = enefunc%qmmm%qm_dipole

        call print_minfo_grad(fname, nat, energy_ref, grad_ref, dipole_ref)

      end if

      rmsg = 0.0_wp
      do i = 1, nat
        iatom = vibatom_id(i)
        do j = 1,3
           rmsg = rmsg + dynvars%force(j,iatom)*dynvars%force(j,iatom)
        end do
      end do
      rmsg = sqrt(rmsg/real(nat*3))

    end if

    k = 1
    do i = 1, nat
      iatom = vibatom_id(i)
      delta = vibration%diff_stepsize/sqrt(mass(iatom))

      do j = 1, 3

        x0 = coord(j,iatom)

        ! ------------------------------------
        ! Compute the energy / force of +delta
        count = count + 1

        if (mod(count,vibration%nreplica) == replicaid-1) then

          dynvars%step = count
          coord(j,iatom) = x0 + delta
  
          write(num,'(i5.5)') count
          basename = 'vib'//num
          fname = trim(minfo_folder)//'/'//trim(basename)
  
          inquire(file=trim(fname)//'.minfo', exist = ex)
          if (.not. ex) then
            call compute_energy(molecule, enefunc, pairlist, boundary,  &
                              .true.,                &
                              .false.,               &
                              dynvars%coord,         &
                              dynvars%energy,        &
                              dynvars%temporary,     &
                              dynvars%force,         &
                              dynvars%virial,        &
                              dynvars%virial_extern)
  
            if(replica_main_rank) &
              write(MsgOut,'(6x,"Done for ",i4,x,a6,"+",a1,6x,"replicaID = ",i5)') &
                    i, molecule%atom_name(iatom), xyz(j), replicaid

            call output_vib(output, molecule, enefunc, vibration, boundary, dynvars)

            energy1 =  dynvars%energy%total / CONV_UNIT_ENE
            k1 = 1
            do i1 = 1, nat
              do j1 = 1, 3
                grad1(k1) = -dynvars%force(j1,vibatom_id(i1)) / CONV_UNIT_FORCE
                k1 = k1 + 1
              end do
            end do
            dipole1 = enefunc%qmmm%qm_dipole
  
            call print_minfo_grad(fname, nat, energy1, grad1, dipole1)
  
          end if
        end if

        ! ------------------------------------
        ! Compute the energy / force of -delta
        !
        count = count + 1

        if (mod(count,vibration%nreplica) == replicaid-1) then

          dynvars%step = count
          coord(j,iatom) = x0 - delta
  
          write(num,'(i5.5)') count
          basename = 'vib'//num
          fname = trim(minfo_folder)//'/'//trim(basename)
  
          inquire(file=trim(fname)//'.minfo', exist = ex)
          if (.not. ex) then
            call compute_energy(molecule, enefunc, pairlist, boundary,  &
                              .true.,                &
                              .false.,               &
                              dynvars%coord,         &
                              dynvars%energy,        &
                              dynvars%temporary,     &
                              dynvars%force,         &
                              dynvars%virial,        &
                              dynvars%virial_extern)
  
            if(replica_main_rank) &
              write(MsgOut,'(6x,"Done for ",i4,x,a6,"-",a1,6x,"replicaID = ",i5)') &
                    i, molecule%atom_name(iatom), xyz(j), replicaid

            call output_vib(output, molecule, enefunc, vibration, boundary, dynvars)
  
            energy2 =  dynvars%energy%total / CONV_UNIT_ENE
            k1 = 1
            do i1 = 1, nat
              do j1 = 1, 3
                grad2(k1) = -dynvars%force(j1,vibatom_id(i1)) / CONV_UNIT_FORCE
                k1 = k1 + 1
              end do
            end do
            dipole2 = enefunc%qmmm%qm_dipole
  
            call print_minfo_grad(fname, nat, energy2, grad2, dipole2)
  
          end if
        end if

        coord(j,iatom) = x0

        k = k + 1

      end do
    end do

#ifdef MPI
    call mpi_barrier(mpi_comm_world, ierr)
#endif

    if(main_rank) then
      write(MsgOut,*)
      write(MsgOut,'(''  RMSD of the gradient at the input geometry = '',d15.6, &
                   & '' [kcal mol-1 Angs-1]'')') rmsg

      if(rmsg > 0.35_wp) then
        write(MsgOut,'(4x,40(''=''))')
        write(MsgOut,'(4x,''!! Warning !!'')')
        write(MsgOut,'(4x,''RMSG is too large for the following vibrational '' ,/, &
                    &  4x,''analysis to be valid. Further miminization until '',/, &
                    &  4x,''RMSG < 0.35 is highly recommended.'')')
        write(MsgOut,'(4x,40(''=''))')
      end if
      write(MsgOut,*)
    end if

    if(vibration%dryrun) then
      deallocate(grad_ref, grad1, grad2, hessian, dipole_derivative)
      return
    endif

    ! Now calculate the Hessian and Dipole deriv.
    if (replicaid == 1) then
      count = 0
      k = 1
      do i = 1, nat
        iatom = vibatom_id(i)
        delta = vibration%diff_stepsize/sqrt(mass(iatom))

        do j = 1, 3

          ! ------------------------------------
          ! Compute the energy / force of +delta
          count = count + 1

          write(num,'(i5.5)') count
          basename = 'vib'//num
          fname = trim(minfo_folder)//'/'//trim(basename)
          call read_minfo_grad(fname, nat, energy1, grad1, dipole1)

          ! ------------------------------------
          ! Compute the energy / force of -delta
          !
          count = count + 1

          write(num,'(i5.5)') count
          basename = 'vib'//num
          fname = trim(minfo_folder)//'/'//trim(basename)
          call read_minfo_grad(fname, nat, energy2, grad2, dipole2)

          ! Calculate the Hessian matrix and dipole derivatives
          hessian(:,k) = (grad1 - grad2)/delta*HALF * CONV_UNIT_LEN
          dipole_derivative(:,k) = (dipole1 - dipole2)/delta*HALF * CONV_UNIT_LEN

          k = k + 1

        end do
      end do
    end if

    if(main_rank) write(MsgOut,*)

#ifdef MPI
    call mpi_bcast(hessian, nat3*nat3,        &
                     mpi_wp_real, 0, mpi_comm_world, ierror)
    call mpi_bcast(dipole_derivative, 3*nat3, &
                     mpi_wp_real, 0, mpi_comm_world, ierror)
#endif

    ! Symmetrize and pack the Hessian matrix 
    allocate(hess_pack(nat3*(nat3+1)/2))
    k = 1
    do k1 = 1, nat3
      do k2 = 1, k1
         hess_pack(k) = (hessian(k2,k1) + hessian(k1,k2))*HALF
         k = k + 1
      end do
    end do

    ! Mass-weight the Hessian matrix and dipole derivatives
    allocate(hess_mw_pack(nat3*(nat3+1)/2), dd_mw(3,nat3), sq_mass3(nat3))
    k = 1
    do i = 1, nat
      do j = 1, 3
         sq_mass3(k) = sqrt(mass(vibatom_id(i))*ELMASS)
         k = k + 1
      end do
    end do

    k = 1
    do k1 = 1, nat3
      do k2 = 1, k1
        hess_mw_pack(k) = hess_pack(k) / sq_mass3(k1) / sq_mass3(k2)
        k = k + 1
      end do
    end do

    do k = 1, nat3
      dd_mw(:,k) = dipole_derivative(:,k) / sq_mass3(k)
    end do

    allocate(vec(nat3,nat3), omega(nat3))

    call diagonalize(nat3, nat3, hess_mw_pack, vec, omega)

    do k = 1, nat3
      if(omega(k) > 0.0_wp) then
        omega(k) = sqrt( omega(k)) * HARTREE_WAVENUM
      else
        omega(k) =-sqrt(-omega(k)) * HARTREE_WAVENUM
      endif
    end do

    ! Transform the dipole_derivatives in terms of normal coordinates
    allocate(dd2(3,nat3), infrared(nat3))
    do i = 1, nat3
      dd2(:,i) = 0.0_wp
      do j = 1, nat3
        dd2(:,i) = dd2(:,i) + dd_mw(:,j)*vec(j,i)
      end do
    end do

    do i = 1, nat3
      infrared(i) = 0.0_wp
      do j = 1, 3
         infrared(i) = infrared(i) + dd2(j,i)*dd2(j,i)
      end do
      infrared(i) = infrared(i) * PI  / 3.0_wp / VLIGHT_IN_AU / VLIGHT_IN_AU &
                  * CONV_UNIT_LEN * 1.0e-13_wp * AVOGADRO
    end do

    ! print normal modes
    ncolomn = 5
    nc = (nat3-mod(nat3,ncolomn))/ncolomn

    if(main_rank) then
      write(MsgOut,'(''  Harmonic frequencies and normal displacement vectors'')')
      k=1
      do i = 1, nc
        write(MsgOut,'(6x,''    mode  '',$)')
        do j = k, k + ncolomn - 1
          write(MsgOut,'(i12,$)') j
        end do
        write(MsgOut,*)

        write(MsgOut,'(6x,''    omega    '',$)')
        do j = k, k + ncolomn - 1
          write(MsgOut,'(f12.4,$)') omega(j)
        end do
        write(MsgOut,*)

        write(MsgOut,'(6x,''   IR int.   '',$)')
        do j = k, k + ncolomn - 1
          write(MsgOut,'(f12.4,$)') infrared(j)
        end do
        write(MsgOut,*)

        k1 = 1
        do i1 = 1, nat
          iatom = vibatom_id(i1)
          do j1 = 1, 3
            write(MsgOut,'(6x,i4,x,a6,x,a1,$)') i1,molecule%atom_name(iatom),xyz(j1)
            do j = k, k + ncolomn - 1
              write(MsgOut,'(f12.4,$)') vec(k1,j)
            end do
            write(MsgOut,*)
            k1 = k1 + 1
          end do
        end do

        write(MsgOut,*)
        k = k + ncolomn
      end do

      if(mod(nat3,ncolomn) /= 0) then
        write(MsgOut,'(6x,''    mode  '',$)')
        do j = k, nat3
          write(MsgOut,'(i12,$)') j
        end do
        write(MsgOut,*)

        write(MsgOut,'(6x,''    omega    '',$)')
        do j = k, nat3
          write(MsgOut,'(f12.4,$)') omega(j)
        end do
        write(MsgOut,*)

        write(MsgOut,'(6x,''   IR int.   '',$)')
        do j = k, nat3
          write(MsgOut,'(f12.4,$)') infrared(j)
        end do
        write(MsgOut,*)

        k1 = 1
        do i1 = 1, nat
          iatom = vibatom_id(i1)
          do j1 = 1, 3
            write(MsgOut,'(6x,i4,x,a6,x,a1,$)') i1,molecule%atom_name(iatom),xyz(j1)
            do j = k, nat3
              write(MsgOut,'(f12.4,$)') vec(k1,j)
            end do
            write(MsgOut,*)
            k1 = k1 + 1
          end do
        end do

        write(MsgOut,*)
      end if
    end if

    if(main_rank) &
      call print_minfo(vibration, molecule, dynvars, nat, energy_ref, grad_ref, &
                       hess_pack, dipole_ref, dipole_derivative, omega, vec)

    deallocate(hess_pack, hess_mw_pack, dd_mw, dd2, infrared)
    deallocate(grad_ref, grad1, grad2, hessian, dipole_derivative)
    deallocate(vec, omega)

    return

  end subroutine harmonic_analysis

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    diagonalize
  !> @brief        diagonalize a symmetric matrix using Lapack (dspevx)
  !! @authors      KY
  !! @param[in]    n : dimension of the matrix
  !! @param[in]    m : number of solution
  !! @param[in]    H(n*(n+1)/2) : The matrix to be diagonalized (upper half)
  !! @param[out]   C(n,m)       : The eigenvectors
  !! @param[out]   E(n)         : The eigenvalues
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine diagonalize(n, m, H, C, E)

    integer :: n,m,m1
    real(wp), dimension(n*(n+1)/2) :: H
    real(wp), dimension(n,m) :: C
    real(wp), dimension(n) :: E

    character :: jobz,range,uplo
    real(wp)  :: vl,vu,abstol
    integer   :: il,iu,ldz,info
    integer,  dimension(:), allocatable :: ifail,iwork
    real(wp), dimension(:), allocatable :: work

      !write(6,*) n,m
      !do i = 1, n
      !   write(6,'(11f12.6)') (H(j),j=i*(i-1)/2+1,i*(i+1)/2)
      !dnd do

      allocate(work(10*n),ifail(n),iwork(10*n))

      jobz='V'
      uplo='U'
      vl=0.0_wp
      vu=0.0_wp
      il=0
      iu=0
      if(n==m) then
         range='A'
      else
         range='I'; il=1; iu=m
      endif

      abstol=0.0_wp
      ldz=n

#ifdef LAPACK
      m1=0
      ifail=0; info=0
      call dspevx(jobz,range,uplo,n,H,vl,vu,il,iu,abstol,m1,E, &
                  C,ldz,work,iwork,ifail,info)

#else
      call error_msg('diagonalize> ERROR: This subroutine needs LAPACK.')
#endif

      deallocate(work,ifail,iwork)

      if(info==0) return

      if(info<0) then
         write(MsgOut,'(''ERROR IN '',i3,''TH PARAMETER'')') info
      else
         write(MsgOut,'(3x,i3,''EIGENVECTORS FAILED TO CONVERGE'')') info
      endif

  end subroutine diagonalize

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    print_minfo
  !> @brief        print the information to minfo format
  !! @authors      KY
  !! @param[in] vibration : vibration information
  !! @param[in] molecule : molecule information
  !! @param[in] dynvars  : dynamic variables information
  !! @param[in] nat      : number of atoms for vibrational analysis
  !! @param[in] energy   : energy
  !! @param[in] grad     : gradient(3*nat)
  !! @param[in] grad     : hessian matrix(3*nat, 3*nat)
  !! @param[in] dipole   : dipole moment (3)
  !! @param[in] dipole_derivative   : dipole moment derivatives (3, nat*3)
  !! @param[in] omega    : harmonic frequency (nat*3)
  !! @param[in] vec      : normal mode displacement vector (nat*3, nat*3)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine print_minfo(vibration, molecule, dynvars, nat, energy, grad, &
                     hess, dipole, dipole_derivative, omega, vec)

    ! formal arguments
    type(s_vibration), target, intent(inout) :: vibration
    type(s_molecule),  target, intent(inout) :: molecule
    type(s_dynvars),   target, intent(inout) :: dynvars

    integer :: nat
    real(wp) :: energy
    real(wp) :: grad(3,nat), hess(nat*3*(nat*3+1)/2)
    real(wp) :: dipole(3), dipole_derivative(3,nat*3)
    real(wp) :: omega(nat*3), vec(nat*3, nat*3)

    ! local variables
    real(wp), pointer     :: coord(:,:)
    integer , pointer     :: vibatom_id(:)

    integer               :: ifile
    integer               :: i, j, k, iatom, ncolumn, isize
    integer               :: atomic_no
    real(wp)              :: mm_mass

    integer               :: nd
    integer, allocatable  :: domain_nat(:), domain_nf(:)
    integer, allocatable  :: domain_idx(:,:)

    ! use pointers
    !
    nat        =  vibration%vib_natoms
    vibatom_id => vibration%vibatom_id
    coord      => dynvars%coord

    ifile   = 10
    ncolumn = 5

    call open_file(ifile, vibration%minfofile, IOFileOutputReplace)
    write(ifile,'(''# minfo File version 2:'')')
    write(ifile,'(''#'')')
    write(ifile,'(''[ Atomic Data ]'')')
    write(ifile,'(i5)') nat
    do i = 1, nat
      iatom = vibatom_id(i)
      mm_mass = molecule%mass(iatom)
      call qm_atomic_number(mm_mass, atomic_no)
      write(ifile,'(a6,'', '',i4,'', '',f12.4,'', '',2(f17.10,'', ''),f17.10)') &
        trim(molecule%atom_name(iatom)),  &
        atomic_no,                  &
        mm_mass,                    &
        coord(:,iatom) / CONV_UNIT_LEN
    end do
    if (vibration%minfo_natoms > 0) then 
      write(ifile,'(i5)') vibration%minfo_natoms
      do i = 1, vibration%minfo_natoms
        iatom = vibration%minfoatom_id(i)
        mm_mass = molecule%mass(iatom)
        call qm_atomic_number(mm_mass, atomic_no)
        write(ifile,'(a6,'', '',i4,'', '',f12.4,'', '',2(f17.10,'', ''),f17.10)') &
          trim(molecule%atom_name(iatom)),  &
          atomic_no,                  &
          mm_mass,                    &
          coord(:,iatom) / CONV_UNIT_LEN
      end do
    end if
    write(ifile,*)

    write(ifile,'(''[ Electronic Data ]'')')
    write(ifile,'(''Energy'')')
    write(ifile,'(f25.14)') energy

    write(ifile,'(''Gradient'')')
    write(ifile,'(i5)') nat*3
    k = 1
    do i = 1, nat
      do j = 1, 3
        write(ifile,'(es15.8,$)') grad(j,i)
        if(mod(k,ncolumn) == 0 .or. (i == nat .and. j == 3)) then
          write(ifile,*)
        else
          write(ifile,'('', '',$)')
        endif
        k = k + 1
      end do
    end do

    write(ifile,'(''Hessian'')')
    isize = nat*3*(nat*3+1)/2
    write(ifile,'(i0)') isize
    do k = 1, isize
      write(ifile,'(es15.8,$)') hess(k)
      if(mod(k,ncolumn) == 0 .or. k == isize) then
        write(ifile,*)
      else
        write(ifile,'('', '',$)')
      endif
    end do

    write(ifile,'(''Dipole Moment'')')
    isize = 3
    write(ifile,'(i5)') isize
    do k = 1, isize
      write(ifile,'(es15.8,$)') dipole(k)
      if(mod(k,ncolumn) == 0 .or. k == isize) then
        write(ifile,*)
      else
        write(ifile,'('', '',$)')
      endif
    end do

    write(ifile,'(''Dipole Derivative'')')
    isize = 3*nat*3
    write(ifile,'(i5)') isize
    k = 1
    do i = 1, nat*3
      do j = 1, 3
        write(ifile,'(es15.8,$)') dipole_derivative(j, i)
        if(mod(k,ncolumn) == 0 .or. k == isize) then
          write(ifile,*)
        else
          write(ifile,'('', '',$)')
        endif
        k = k + 1
      end do
    end do
    write(ifile,*)

    nd = 1
    allocate(domain_nat(nd), domain_nf(nd), domain_idx(nat,nd))
    domain_nat(1) = nat
    domain_nf(1)  = nat*3
    do i = 1, nat
      domain_idx(i,1) = i
    end do

    write(ifile,'(''[ Vibrational Data ]'')')
    write(ifile,'(''Number of Domain'')')
    write(ifile,'(i4)') nd
    do i = 1, nd
      write(ifile,'(''Domain '',i4)') nd
      write(ifile,'(''Atom Index'')')
      isize = domain_nat(i)
      write(ifile,'(i4)') isize
      do j = 1, isize
        write(ifile,'(i15,$)') domain_idx(j,i)
        if(mod(j,ncolumn) == 0 .or. j == isize) then
          write(ifile,*)
        else
          write(ifile,'('', '',$)')
        endif
      end do

      write(ifile,'(''Local Normal modes'')')
      write(ifile,'(''Vibrational Frequency'')')
      isize = domain_nf(i)
      write(ifile,'(i4)') isize
      do j = 1, isize
        write(ifile,'(es15.8,$)') omega(j)
        if(mod(j,ncolumn) == 0 .or. j == isize) then
          write(ifile,*)
        else
          write(ifile,'('', '',$)')
        endif
      end do

      write(ifile,'(''Vibrational vector'')')
      do j = 1, domain_nf(i)
        write(ifile,'(''Mode '',i6)') j
        isize = domain_nat(i)*3
        write(ifile,'(i4)') isize
        do k = 1, isize
          write(ifile,'(es15.8,$)') vec(k,j)
          if(mod(k,ncolumn) == 0 .or. k == isize) then
            write(ifile,*)
          else
            write(ifile,'('', '',$)')
          endif
        end do
      end do

    end do

    write(ifile,'(/)')

    call close_file(ifile)

  end subroutine print_minfo

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    print_minfo_grad
  !> @brief        print the information to minfo format
  !! @authors      KY
  !! @param[in] basename : basename of the minfo file
  !! @param[in] nat      : number of atoms for vibrational analysis
  !! @param[in] energy   : energy
  !! @param[in] grad     : gradient(3*nat)
  !! @param[in] dipole   : dipole moment (3)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine print_minfo_grad(basename, nat, energy, grad, dipole)

    ! formal arguments
    character(256) :: basename
    integer        :: nat
    real(wp)       :: energy, grad(3,nat), dipole(3)

    integer        :: ifile
    integer        :: i, j, k, ncolumn, isize

    ifile   = 10
    ncolumn = 5

    call open_file(ifile, trim(basename)//'.minfo', IOFileOutputNew)
    write(ifile,'(''# minfo File version 2:'')')
    write(ifile,'(''#'')')
    write(ifile,'(''[ Electronic Data ]'')')
    write(ifile,'(''Energy'')')
    write(ifile,'(f25.14)') energy

    write(ifile,'(''Gradient'')')
    write(ifile,'(i5)') nat*3
    k = 1
    do i = 1, nat
      do j = 1, 3
        write(ifile,'(es15.8,$)') grad(j,i)
        if(mod(k,ncolumn) == 0 .or. (i == nat .and. j == 3)) then
          write(ifile,*)
        else
          write(ifile,'('', '',$)')
        endif
        k = k + 1
      end do
    end do

    write(ifile,'(''Dipole Moment'')')
    isize = 3
    write(ifile,'(i5)') isize
    do k = 1, isize
      write(ifile,'(es15.8,$)') dipole(k)
      if(mod(k,ncolumn) == 0 .or. k == isize) then
        write(ifile,*)
      else
        write(ifile,'('', '',$)')
      endif
    end do

    write(ifile,'(/)')

    call close_file(ifile)

  end subroutine print_minfo_grad


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_minfo_grad
  !> @brief        read the information to minfo format
  !! @authors      KY
  !! @param[in]    basename : basename of the minfo file
  !! @param[in]    nat      : number of atoms for vibrational analysis
  !! @param[inout] energy   : energy
  !! @param[inout] grad     : gradient(3*nat)
  !! @param[inout] dipole   : dipole moment (3)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_minfo_grad(basename, nat, energy, grad, dipole)

    ! formal arguments
    character(256) :: basename
    integer        :: nat
    real(wp)       :: energy, grad(nat*3), dipole(3)

    integer        :: ifile
    integer        :: i, j, k, ncolumn, isize

    ifile   = 10
    ncolumn = 5

    call open_file(ifile, trim(basename)//'.minfo', IOFileInput)
    read(ifile,*)
    read(ifile,*)
    read(ifile,*)
    read(ifile,*)
    read(ifile,'(f20.14)') energy

    read(ifile,*)
    read(ifile,*)
    read(ifile,*) grad(1:nat*3)

    read(ifile,*)
    read(ifile,*)
    read(ifile,*) dipole(1:3)

    call close_file(ifile)

  end subroutine read_minfo_grad

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    calc_gridpoints
  !> @brief        calculate the energy and gradient at the grid points
  !> @authors      KY
  !! @param[inout] molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !! @param[inout] dynvars  : dynamic variables information
  !! @param[inout] vibration : vibration information
  !! @param[inout] pairlist : non-bond pair list information
  !! @param[inout] boundary : boundary conditions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine calc_gridpoints(molecule, enefunc, dynvars, vibration, &
                     output, pairlist, boundary)

    ! formal arguments
    type(s_molecule),  target, intent(inout) :: molecule
    type(s_enefunc),   target, intent(inout) :: enefunc
    type(s_dynvars),   target, intent(inout) :: dynvars
    type(s_vibration), target, intent(inout) :: vibration
    type(s_output),            intent(inout) :: output
    type(s_pairlist),          intent(inout) :: pairlist
    type(s_boundary),          intent(inout) :: boundary

    ! local
    type(s_qmmm), pointer :: qmmm
    integer               :: nat, nat3
    integer               :: i, ia, j, ierr
    real(wp), pointer     :: coord(:,:)
    integer , pointer     :: vibatom_id(:)

    character(256) :: folder, minfo_folder, datafile, basename, fname, line
    character(4)   :: fid, fnum, label
    character(7)   :: num
    logical        :: savefile
    integer        :: icount, count, replicaid

    integer        :: ifile, idatafile
    logical        :: ex, thisrun_done, restart
    real(wp)       :: energy, grad(3,vibration%vib_natoms), dipole(3)


    nat        =  vibration%vib_natoms
    nat3       =  nat*3
    vibatom_id => vibration%vibatom_id
    coord      => dynvars%coord
    qmmm       => enefunc%qmmm

    ! replica id
    replicaid = my_country_no + 1
    count = 0

    if(qmmm%do_qmmm) then
      icount        = 0
      folder        = 'qmmm_grid_'//trim(vibration%gridfolderID)
      qmmm%ene_only = vibration%grid_ene_only
    endif

    if(vibration%grid_ene_only) then
      ! set datafile
      write(num,'(i0)') my_country_no
      datafile = trim(vibration%datafile)//'_'//trim(num)

      ! ignore qm error
      qmmm%ignore_qm_error = .true.

      if(main_rank) &
        write(MsgOut,*) 'Compute energy at grid points: &
         &data written to [ '//trim(vibration%datafile)//' ]'

    else
      ! create a folder to save minfo files
      minfo_folder = 'minfo.files'
      if(main_rank) then
        call system('mkdir -p '//trim(minfo_folder)//' >& /dev/null')
      end if

      if(main_rank) &
        write(MsgOut,*) 'Compute energy at grid points: &
         &minfo files created in [ '//trim(minfo_folder)//' ]'

    endif

    call open_file(ifile, trim(vibration%gridfile), IOFileInput)

    restart = .true.
    do while(.true.)
      read(ifile,*,end=20)
      read(ifile,'(a)') basename

      do i = 1, nat
        ia = vibatom_id(i)
        read(ifile,*) label,coord(:,ia)
      end do

      count = count + 1

      if (mod(count,vibration%nreplica) == replicaid-1) then

        ! check for restart
        if (restart) then
          thisrun_done = .false.
          if (vibration%grid_ene_only) then
            fname = trim(datafile)
            inquire(file=trim(datafile),exist=ex)
            if (ex) then
              call open_file(idatafile, trim(fname), IOFileInput)
              do while(.true.)
                read(idatafile,'(a)',end=10) line
                i = index(line,',')
                if(trim(basename) == line(1:i-1)) then
                   thisrun_done = .true.
                   exit
                endif 
              end do
           10 continue
              close(idatafile)
            end if
    
          else
            fname = trim(minfo_folder)//'/'//trim(basename)
            inquire(file=trim(fname)//'.minfo',exist=thisrun_done)
    
          end if
          if(thisrun_done) then
            cycle
          else
            restart = .false.
          end if

        end if

        call compute_energy(molecule, enefunc, pairlist, boundary,  &
                          .true., &
                          .false.,               &
                          dynvars%coord,         &
                          dynvars%energy,        &
                          dynvars%temporary,     &
                          dynvars%force,         &
                          dynvars%virial,        &
                          dynvars%virial_extern)


        if (.not. vibration%dryrun) then
          energy =  dynvars%energy%total / CONV_UNIT_ENE
          do i = 1, nat
            do j = 1, 3
              grad(j,i) = -dynvars%force(j,vibatom_id(i)) / CONV_UNIT_FORCE
            end do
          end do
          dipole = qmmm%qm_dipole

          if (replica_main_rank) then
            write(MsgOut,'(6x,"Done for ",a30," :",4x,"replicaID = ",i5)') &
                  trim(basename), replicaid

            if (vibration%grid_ene_only) then
              if (qmmm%qmmm_error) then
                energy = 0.0_wp
                dipole = 0.0_wp
              end if

              idatafile = get_unit_no()
              fname     = trim(datafile)
              open(idatafile, file=trim(fname), status='unknown', access='append')
              write(idatafile,'(a,'', '',f25.13,'', '',2(es16.9,'',''),es16.9)') &
                    trim(basename), energy, dipole
              close(idatafile)
              
            else
              fname = trim(minfo_folder)//'/'//trim(basename)
              call print_minfo_grad(fname, nat, energy, grad, dipole)

            endif
          endif
        endif
      endif

    end do

 20 continue
    call close_file(ifile)

#ifdef MPI
    call mpi_barrier(mpi_comm_world, ierr)
#endif

    ! combine the grid data to one file
    if (vibration%grid_ene_only .and. main_rank) then
      
      call open_file(ifile, trim(vibration%datafile), IOFileOutputReplace)

      do i = 1, vibration%nreplica

        write(num,'(i0)') i-1
        datafile = trim(vibration%datafile)//'_'//trim(num)
        call open_file(idatafile, trim(datafile), IOFileInput)
        do while(.true.)
          read(idatafile,'(a)',end=30) line
          j = index(line,',')
          read(line(j+1:),*) energy
          if(abs(energy) > 1.e-10)  write(ifile,'(a)') trim(line)
        end do
     30 continue

        close(idatafile,status='delete')

      end do

      call close_file(ifile)

    end if

  end subroutine calc_gridpoints

end module at_vibration_mod

