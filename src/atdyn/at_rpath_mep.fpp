!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_rpath_mep_mod
!> @brief   Minimum energy path search using replicas
!! @authors Yoshinobu Akinaga (YA), Kiyoshi Yagi (KY)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_rpath_mep_mod

  use at_energy_mod
  use at_dynamics_str_mod
  use at_dynvars_str_mod
  use at_dynvars_mod
  use at_rpath_str_mod
  use at_boundary_str_mod
  use at_boundary_mod
  use at_pairlist_str_mod
  use at_pairlist_mod
  use at_enefunc_str_mod
  use at_minimize_str_mod
  use at_minimize_mod
  use at_output_str_mod
  use at_output_mod
  use molecules_str_mod
  use fileio_rst_mod
  use fileio_rstmep_mod
  use select_mod
  use select_atoms_mod
  use select_atoms_str_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  use math_libs_mod
#ifdef MPI
  use mpi
#endif

  implicit none
  private

  ! subroutines
  public  :: setup_rpath_mep
  public  :: setup_mepatoms
  public  :: setup_mepatoms_qmmm
  public  :: run_rpath_mep
  public  :: run_rpath_string
  public  :: run_rpath_neb
  private :: evolve_mep
  private :: reparametrize_mep
  private :: interpolation
  private :: energy_and_force
  private :: check_convergence
  private :: check_convergence_neb
  private :: trans_mass_weight_coord
  private :: backtrans_mass_weight_coord
  private :: neb_energy_and_force
  private :: neb_force

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_rpath_mep
  !> @brief        setup RPATH
  !! @authors      YA, KY
  !! @param[in]    mepatm_select_index : index of MEP atoms
  !! @param[in]    sel_info   : selector input information
  !! @param[in]    rst        : restart information
  !! @param[in]    rstmep     : rstmep data
  !! @param[in]    molecule   : molecule information
  !! @param[in]    qmmm       : qmmm information
  !! @param[inout] minimize   : minimize information
  !! @param[inout] dynvars    : dynamic variables
  !! @param[inout] rpath      : RPATH information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_rpath_mep(mepatm_select_index, sel_info, rst, rstmep, &
                             molecule, qmmm, minimize, dynvars, rpath)

    ! formal arguments
    character(256),          intent(in)    :: mepatm_select_index
    type(s_sel_info),        intent(in)    :: sel_info
    type(s_rst),             intent(in)    :: rst
    type(s_rstmep),          intent(in)    :: rstmep
    type(s_molecule),        intent(in)    :: molecule
    type(s_qmmm),            intent(inout) :: qmmm
    type(s_minimize),        intent(inout) :: minimize
    type(s_dynvars),         intent(inout) :: dynvars
    type(s_rpath),           intent(inout) :: rpath

    ! local variables
    integer   :: i, ii
    integer   :: replicaid
    integer   :: atomid, iatom
    integer   :: errorcode, ierr
    logical   :: err


    ! define atoms
    !
    call setup_mepatoms(mepatm_select_index, sel_info, molecule, rpath)

    if (qmmm%do_qmmm) then
      call setup_mepatoms_qmmm(molecule, qmmm, rpath)
      if (minimize%macro) then
        rpath%opt_micro = .true.
      else
        rpath%opt_micro = .false.
      end if
    else
      rpath%opt_micro = .false.
    end if
    call add_fixatm(rpath%mep_natoms, rpath%mepatom_id, molecule, minimize)

    ! allocate memory
    rpath%dimension = rpath%mep_natoms * 3
    if (qmmm%do_qmmm) then
      call alloc_rpath_mep(rpath, rpath%dimension, rpath%nreplica, qmmm%qm_natoms)
    else
      call alloc_rpath_mep(rpath, rpath%dimension, rpath%nreplica)
    end if


    ! set initial images
    replicaid = my_country_no + 1

    rpath%is_qmcharge = .false.
    if (rstmep%restart) then
      if (rstmep%mep_natoms /= rpath%mep_natoms) &
        call error_msg("Setup_Rpath_MEP> The number of MEP atoms don't match &
                       &in the ctrlfile and rstmep.")
      rpath%mep_coord(:,replicaid) = rstmep%mep_coord

      if (qmmm%do_qmmm) then
        if (rstmep%qm_natoms /= qmmm%qm_natoms) &
          call error_msg("Setup_Rpath_MEP> The number of QM atoms don't match &
                         &in the ctrlfile and rstmep.")
        qmmm%qm_charge = rstmep%qm_charge
        rpath%is_qmcharge = .true.

      end if 

    else
      ii = 0
      do i = 1, rpath%mep_natoms
        iatom = rpath%mepatom_id(i)
        rpath%mep_coord(ii+1:ii+3,replicaid) = dynvars%coord(1:3,iatom)
        ii = ii + 3
      end do

    end if

    !dbg call mpi_barrier(mpi_comm_world, ierr)
    !dbg write(MsgOut,'("ID",i4,"> ",100f8.4)') my_country_no, qmmm%qm_charge
    !dbg call mpi_barrier(mpi_comm_world, ierr)
    !dbg call mpi_abort(mpi_comm_world, errorcode, ierr)

    ! write the summary of setup
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Setup_Rpath_MEP> Rpath information'
      write(MsgOut,'(A)') ''
      write(MsgOut,'(A20,I10)') '  dimension       = ', rpath%dimension
      write(MsgOut,'(A)') ''

      ! Punch out atom info.
      !
      write(MsgOut,'(a)') "  Atoms involved in MEP search"
      do i = 1, rpath%mep_natoms
        atomid   = rpath%mepatom_id(i)
        write(MsgOut,'(i6,x,a4,x,i6,x,a4,x,a6,x,a6)') &
            atomid, &
            molecule%segment_name(atomid), &
            molecule%residue_no(atomid),   &
            molecule%residue_name(atomid), &
            molecule%atom_name(atomid),    &
            molecule%atom_cls_name(atomid)
      end do
      write(MsgOut,'(a,i0)') "  number of atoms in MEP search = ", rpath%mep_natoms
      write(MsgOut, '(a)') ' '

      if (rstmep%restart) then
        write(MsgOut, '(a)') ' '
        write(MsgOut,'(a)') '  Restart retrieved from mep files.'
        write(MsgOut, '(a)') ' '
      end if

    end if

    return

  end subroutine setup_rpath_mep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_mepatoms
  !> @brief        define MEP atoms
  !! @authors      YA
  !! @param[in]    mepatm_select_index : index of MEP atoms
  !! @param[in]    sel_info   : selector input information
  !! @param[in]    molecule   : molecule information
  !! @param[out]   rpath      : rpath parameters
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_mepatoms(mepatm_select_index, sel_info, molecule, rpath)

    ! formal arguments
    character(256),      intent(in)    :: mepatm_select_index
    type(s_sel_info),    intent(in)    :: sel_info
    type(s_molecule),    intent(in)    :: molecule
    type(s_rpath),       intent(inout) :: rpath

    ! local variables
    integer                :: igroup, ngroup, natom, i, j, offset, temp
    integer, allocatable   :: group_list(:)
    type(s_selatoms), allocatable :: selatoms(:)


    ! Number of atoms in MEP
    !
    ngroup = split_num(trim(mepatm_select_index))

    allocate(group_list(ngroup))
    call split(ngroup, ngroup, mepatm_select_index, group_list)

    allocate(selatoms(ngroup))

    natom = 0
    do i = 1, ngroup
      igroup = group_list(i)
      call select_atom(molecule, sel_info%groups(igroup), selatoms(i))
      natom = natom + size(selatoms(i)%idx)
    end do

    rpath%mep_natoms = natom

    ! List of atoms
    !
    allocate(rpath%mepatom_id(rpath%mep_natoms))

    offset = 0
    do i = 1, ngroup
      igroup = group_list(i)
      natom = size(selatoms(i)%idx)
      rpath%mepatom_id(offset+1:offset+natom) = selatoms(i)%idx(1:natom)
      offset = offset + natom
    end do

    deallocate(selatoms)
    deallocate(group_list)

    ! sort atom indices in ascending order
    !
    do i = rpath%mep_natoms, 2, -1
      do j = 1, i - 1
        if (rpath%mepatom_id(j) > rpath%mepatom_id(j+1)) then
          temp = rpath%mepatom_id(j)
          rpath%mepatom_id(j)   = rpath%mepatom_id(j+1)
          rpath%mepatom_id(j+1) = temp
        end if
      end do
    end do

    return

  end subroutine setup_mepatoms

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_mepatoms_qmmm
  !> @brief        define MEP atoms for QM/MM jobs
  !! @authors      KY
  !! @param[in]    molecule   : molecule information
  !! @param[in]    qmmm       : QMMM information
  !! @param[out]   rpath      : rpath parameters
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_mepatoms_qmmm(molecule, qmmm, rpath)

    ! formal arguments
    type(s_molecule),  intent(in)    :: molecule
    type(s_qmmm),      intent(in)    :: qmmm
    type(s_rpath),     intent(inout) :: rpath

    ! local variables
    integer              :: i, j, k, id, nad, ntmp
    logical              :: error
    logical, allocatable :: found(:)
    integer, allocatable :: add(:), tmp(:)

    allocate(found(qmmm%qm_natoms))
    error = .false.
    do i = 1, qmmm%qm_natoms
      id = qmmm%qmatom_id(i)
      found(i) = .false.
      do j = 1, rpath%mep_natoms
        if(id == rpath%mepatom_id(j)) then
          found(i) = .true.
          exit
        end if
      end do
      if(.not. found(i)) error = .true.
    end do

    if(error) then
      if(main_rank) then
        write(MsgOut,'(a)') "Setup_MEP> Fatal error while setting up MEP."
        write(MsgOut,'(a)') "Setup_MEP> These QM atoms are not in MEP:"
        do i = 1, qmmm%qm_natoms
          if(.not. found(i)) then
            id = qmmm%qmatom_id(i)
            write(MsgOut,'(i6,x,a4,x,i6,x,a4,x,a6,x,a6)') &
                id, &
                molecule%segment_name(id), &
                molecule%residue_no(id),   &
                molecule%residue_name(id), &
                molecule%atom_name(id),    &
                molecule%atom_cls_name(id)
          end if
        end do
        write(MsgOut,'(a)') "Setup_MEP> check mepatm_select_index."
      end if
      call error_msg("Setup_MEP> Abort with error.")
    end if

    deallocate(found)

    if (qmmm%num_qmmmbonds > 0) then

      ! search boundary MM atoms in MEP atoms
      allocate(found(qmmm%num_qmmmbonds), add(qmmm%num_qmmmbonds))
      nad = 0
      do i = 1, qmmm%num_qmmmbonds
        id = qmmm%qmmmbond_list(2,i)
        found(i) = .false.
        do j = 1, rpath%mep_natoms
          if(id == rpath%mepatom_id(j)) then
            found(i) = .true.
            exit
          end if
        end do
        if(.not. found(i)) then
          nad = nad + 1
          add(nad) = id
        end if
      end do

      ! add boundary MM atoms to MEP atoms
      if(nad > 0) then
        ntmp = rpath%mep_natoms
        allocate(tmp(rpath%mep_natoms))
        tmp = rpath%mepatom_id

        deallocate(rpath%mepatom_id)
        rpath%mep_natoms = rpath%mep_natoms + nad
        allocate(rpath%mepatom_id(rpath%mep_natoms))

        j = 1
        k = 1
        do i = 1, rpath%mep_natoms
          if (add(j) < tmp(k)) then
            rpath%mepatom_id(i) = add(j)
            j = j + 1
            if (j > nad) then
              rpath%mepatom_id(i+1:) = tmp(k:)
              exit
            end if

          else if (tmp(k)  < add(j)) then
            rpath%mepatom_id(i) = tmp(k)
            k = k + 1
            if (k > ntmp) then
              rpath%mepatom_id(i+1:) = add(j:)
              exit
            end if 
          end if 
       
        end do

        deallocate(tmp)
      end if

      deallocate(found, add)
    end if

    return

  end subroutine setup_mepatoms_qmmm

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    run_rpath_mep
  !> @brief        control string method
  !! @authors      TM, YK, YA, KY
  !! @param[inout] output      : output information
  !! @param[inout] molecule    : molecule information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] dynvars     : dynamics variables information
  !! @param[inout] minimize    : minimize information
  !! @param[inout] dynamics    : dynamics information
  !! @param[inout] pairlist    : pairlist information
  !! @param[inout] boundary    : boundary conditions information
  !! @param[inout] rpath       : RPATH information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine run_rpath_mep(output, molecule, enefunc, dynvars, minimize, &
                           dynamics, pairlist, boundary, rpath)

    ! formal arguments
    type(s_output),           intent(inout) :: output
    type(s_molecule),         intent(inout) :: molecule
    type(s_enefunc),  target, intent(inout) :: enefunc
    type(s_dynvars),  target, intent(inout) :: dynvars
    type(s_minimize), target, intent(inout) :: minimize
    type(s_dynamics),         intent(inout) :: dynamics
    type(s_pairlist),         intent(inout) :: pairlist
    type(s_boundary),         intent(inout) :: boundary
    type(s_rpath),            intent(inout) :: rpath

    ! local variables
    logical                :: conv
    integer                :: niter
    integer                :: replicaid
    integer                :: i, ii, iatom

    real(wp)    , pointer  :: coord(:,:)


    return

  end subroutine run_rpath_mep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    evolve_mep
  !> @brief        evolve path (MEP version)
  !! @authors      YK, YM, YA, KY
  !! @param[inout] dynvars  : dynamical variables information
  !! @param[inout] rpath    : RPATH information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine evolve_mep(rpath)

    ! formal arguments
    type(s_rpath), target, intent(in) :: rpath

    ! local variables
    integer           :: dimension
    real(wp)          :: delta
    real(wp), pointer :: force(:)
    real(wp), pointer :: image(:,:)
    integer           :: i, repid

    if(.not. replica_main_rank) return

    dimension =  rpath%dimension
    delta     =  rpath%delta
    force     => rpath%force
    image     => rpath%mep_coord

    repid = my_country_no + 1

    ! Evolve images
    !
    do i = 1, dimension
      image(i, repid) = image(i, repid) + delta * force(i)
    end do

    return

  end subroutine evolve_mep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    reparametrize_mep
  !> @brief        reparametrize path
  !! @authors      YA, KY
  !! @param[inout] rpath    : RPATH information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine reparametrize_mep(rpath)

    ! formal arguments
    type(s_rpath),    target, intent(inout) :: rpath

    ! local variables
    integer                  :: i, j, k
    integer                  :: dimno, dimno_i, dimno_j
    integer                  :: repid, repid_i, repid_j

    real(wp),    allocatable :: path(:,:), path_reparm(:,:)
    real(wp),    allocatable :: path_leng(:), path_equi(:)
    real(wp),    allocatable :: nonuniform_mesh(:)
    real(wp),    allocatable :: Vs(:)

    real(dp),        pointer :: before_gather(:), after_gather(:)
    real(wp),    allocatable :: after(:)


    before_gather => rpath%before_gather
    after_gather  => rpath%after_gather

    repid = my_country_no + 1

    allocate(path(rpath%nreplica, rpath%dimension),&
             path_reparm(rpath%nreplica, rpath%dimension),&
             path_leng(rpath%nreplica), nonuniform_mesh(rpath%nreplica))

    do dimno = 1, rpath%dimension
      before_gather(dimno) = rpath%mep_coord(dimno, repid)
    end do

    allocate(after(rpath%nreplica))

#ifdef MPI
    call mpi_gather(before_gather, rpath%dimension, mpi_real8,&
                    after_gather,  rpath%dimension, mpi_real8,&
                    0, mpi_comm_airplane, ierror)
    call mpi_gather(rpath%energy, 1, mpi_wp_real,&
                    after, 1, mpi_wp_real,&
                    0, mpi_comm_airplane, ierror)

    if (main_rank) then
      do repid_i = 1, rpath%nreplica
        do dimno = 1, rpath%dimension
          path(repid_i,dimno) = after_gather((repid_i-1)*rpath%dimension+dimno)
        end do
      end do

      ! calc arc lengths
      path_leng(1) = 0.0_wp
      do repid_i = 2, rpath%nreplica
        path_leng(repid_i) = path_leng(repid_i-1) + &
          sqrt( sum( (path(repid_i,:)-path(repid_i-1,:))**2) )
      end do
      rpath%pathlength = path_leng(rpath%nreplica)

      ! non-uniform normalized mesh along path
      !
      do repid_i = 1, rpath%nreplica
        nonuniform_mesh(repid_i) = path_leng(repid_i) &
                                 / path_leng(rpath%nreplica)
      end do

      ! interpolate images and evaluate energies at uniform mesh
      allocate(Vs(rpath%nreplica))
      do dimno = 1, rpath%dimension
        call interpolation(rpath%nreplica, nonuniform_mesh, path(1,dimno), &
                           Vs, path_reparm(1,dimno))
      end do
      deallocate(Vs)

    end if

    ! broadcast mep coordinates
    call mpi_bcast(path_reparm, rpath%dimension*rpath%nreplica,&
                   mpi_wp_real, 0, mpi_comm_world, ierror)
#endif

    do dimno = 1, rpath%dimension
      rpath%mep_coord(dimno, repid) = path_reparm(repid, dimno)
    end do

    deallocate(after)
    deallocate(path, path_reparm, path_leng, nonuniform_mesh)

    return

  end subroutine reparametrize_mep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    interpolation
  !> @brief        update images
  !! @authors      YA
  !! @param[in]    nn           : # of points
  !! @param[in]    path         : images
  !! @param[in]    qi           : non-uniform mesh
  !! @param[inout] Vs           : work space
  !! @param[out]   path_reparam : updated images
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine interpolation(nimages, qi, path, Vs, path_reparm)

    ! formal arguments
    integer,         intent(in) :: nimages
    real(wp),        intent(in) :: path(nimages)
    real(wp),        intent(in) :: qi(nimages)
    real(wp),        intent(in) :: Vs(nimages)
    real(wp),     intent(inout) :: path_reparm(nimages)

    ! local variables
    integer  :: i
    real(wp) :: tmp
    real(wp) :: coeff(4,nimages-1)

    call cubic_spline_coeff(nimages, qi, path, coeff)

    do i = 1, nimages
      tmp = real(i-1,wp) / real(nimages-1,wp)
      call cubic_spline(nimages, qi, path, coeff, tmp, path_reparm(i))
    end do

    return

  end subroutine interpolation

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    energy_and_force
  !> @brief        compute energy and force for each image
  !! @authors      YA
  !! @param[in]    rpath        : rpath info
  !! @param[in]    molecule     : molecular info
  !! @param[in]    enefunc      : enefunc info
  !! @param[in]    pairlist     : pairlist info
  !! @param[in]    boundary     : boundary info
  !! @param[in]    nonb_ene     : flag for calculate nonbonded energy
  !! @param[inout] dynvars      : dynvars info
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine energy_and_force(ncount, rpath, molecule, enefunc, pairlist, &
                              boundary, nonb_ene, dynvars)

    ! formal arguments
    integer,             intent(in) :: ncount
    type(s_rpath),    intent(inout) :: rpath
    type(s_molecule),    intent(inout) :: molecule  !YA
    type(s_enefunc),  intent(inout) :: enefunc
    type(s_pairlist),    intent(in) :: pairlist
    type(s_boundary),    intent(in) :: boundary
    logical,             intent(in) :: nonb_ene
    type(s_dynvars),  intent(inout) :: dynvars

    ! local variables
    character(256)                  :: folder, basename
    character(5)                    :: num
    logical                         :: savefile
    integer                         :: repid
    integer                         :: dimno, iatom, i
    real(wp)                        :: x, y, x2, y2, r

    ! qmmm settings
    !
    if(enefunc%qmmm%do_qmmm) then
    !  folder   = 'qmmm_mep'
    !  write(num,'(i5.5)') ncount
    !  basename = 'mep'//num
    !  savefile = (mod(ncount,rpath%qmsave_period) == 0)
    !  call set_runtime_qmmm(enefunc%qmmm, folder, basename, savefile)
    !
      enefunc%qmmm%qm_classical = .FALSE.
      if(rpath%opt_micro) enefunc%qmmm%qm_get_esp = .TRUE.
    end if

    ! compute energy and force
    !
    call compute_energy(molecule, enefunc, pairlist, boundary, .true., &
                        .false.,           &
                        dynvars%coord,     &
                        dynvars%energy,    &
                        dynvars%temporary, &
                        dynvars%force,     &
                        dynvars%virial,    &
                        dynvars%virial_extern)

    if (.not. replica_main_rank) return

    ! energy and force
    !
    rpath%energy = dynvars%energy%total
    do i = 1, rpath%mep_natoms
      iatom = rpath%mepatom_id(i)
      rpath%force(i*3-2:i*3) = dynvars%force(1:3, iatom)
    end do

    ! clear force for terminals
    !
    if(rpath%fix_terminal) then
      repid = my_country_no + 1
      if((repid == 1) .or. (repid == rpath%nreplica)) then
        rpath%force = 0.0_wp
      end if
    end if

    return

  end subroutine energy_and_force

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    check_convergence
  !> @brief        
  !! @authors      YA
  !! @param[in]    rpath        : rpath info
  !! @param[in]    molecule     : molecular info
  !! @param[in]    enefunc      : enefunc info
  !! @param[in]    pairlist     : pairlist info
  !! @param[in]    boundary     : boundary info
  !! @param[in]    nonb_ene     : flag for calculate nonbonded energy
  !! @param[inout] dynvars      : dynvars info
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine check_convergence(rpath, niter, conv)

    ! formal arguments
    type(s_rpath),    intent(inout) :: rpath
    integer,             intent(in) :: niter
    logical,            intent(out) :: conv

    ! local variables
    logical                         :: my_conv, conv_e, conv_path
    integer                         :: repid, i
    real(wp)                        :: disp, delta, dmax, sum
    real(wp),           allocatable :: work(:), work1(:)


    repid = my_country_no + 1

    ! Energy convergence
    !
    allocate(work(rpath%nreplica))
    allocate(work1(rpath%nreplica))
    delta = rpath%energy - rpath%energy_prev
    rpath%energy_prev = rpath%energy
    call mpi_allgather(delta, 1, mpi_wp_real, work, 1, mpi_wp_real, &
      mpi_comm_airplane, ierror)
    call mpi_allgather(rpath%energy, 1, mpi_wp_real, work1, 1, mpi_wp_real, &
      mpi_comm_airplane, ierror)
    if (main_rank) then

      dmax = work(1)
      sum = abs(work(1))
      do i = 2, rpath%nreplica
        sum = sum + abs(work(i))
        if (abs(dmax) < abs(work(i))) dmax = work(i)
      end do
      conv_e = abs(dmax) < rpath%tol_energy
      
      write(MsgOut, '(/,"Image ", &
                        "        Energy (kcal/mol)", &
                        "              Relative E.", &
                        "             Energy Conv.")')
      write(MsgOut, '("------", &
                      "-------------------------", &
                      "-------------------------", &
                      "-------------------------")')
      do i = 1, rpath%nreplica
        write(MsgOut, '(i5,X,3F25.10)') i, work1(i), work1(i) - work1(1), work(i)
      end do
      write(MsgOut, '("------", &
        "-------------------------", &
        "-------------------------", &
        "-------------------------")')

      
      write(MsgOut, '(3X, "Energy Conv. (Max) = ", F25.15)') dmax
      write(MsgOut, '(3X, "Path length: current value / variation = ", F15.10, "/", F15.10)') &
        rpath%pathlength, rpath%pathlength - rpath%pathlength_prev

      ! Path-length convergence
      !
      conv_path = abs(rpath%pathlength - rpath%pathlength_prev) < rpath%tol_path
      rpath%pathlength_prev = rpath%pathlength
      conv = conv_e .and. conv_path
      
    end if
    call mpi_bcast(conv, 1, mpi_logical, 0, mpi_comm_world, ierror)

    deallocate(work)
    deallocate(work1)
    
    return

  end subroutine check_convergence

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    check_convergence_neb
  !> @brief        
  !! @authors      YA, KY
  !! @param[in]    rpath        : rpath info
  !! @param[in]    molecule     : molecular info
  !! @param[in]    enefunc      : enefunc info
  !! @param[in]    pairlist     : pairlist info
  !! @param[in]    boundary     : boundary info
  !! @param[in]    nonb_ene     : flag for calculate nonbonded energy
  !! @param[inout] dynvars      : dynvars info
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine check_convergence_neb(rpath, niter, conv)

    ! formal arguments
    type(s_rpath),    intent(inout) :: rpath
    integer,             intent(in) :: niter
    logical,            intent(out) :: conv

    ! local variables
    character(3)                    :: stat
    logical                         :: conv_image(rpath%nreplica)
    integer                         :: repid, i, image
    real(wp)                        :: disp, delta, absg, dmax, rmsg_max, maxg_max
    real(wp),           allocatable :: work(:), rmsg(:), maxg(:)


    repid = my_country_no + 1

    allocate(rmsg(rpath%nreplica))
    allocate(maxg(rpath%nreplica))
    allocate(work(rpath%nreplica))

    delta = rpath%energy - rpath%energy_prev
    rpath%energy_prev = rpath%energy
    call mpi_allgather(delta, 1, mpi_wp_real, work, 1, mpi_wp_real, &
      mpi_comm_airplane, ierror)

    ! RMSG and MaxG
    !
    rmsg(:) = 0.0_wp
    maxg(:) = 0.0_wp
    do image = 1, rpath%nreplica
      do i = 1, rpath%mep_natoms * 3
        rmsg(image) = rmsg(image) + rpath%mep_force(i,image)**2
        absg = abs(rpath%mep_force(i,image))
        if (absg > maxg(image)) maxg(image) = absg
      end do
      rmsg(image) = sqrt(rmsg(image) / real(rpath%mep_natoms*3,wp))
      conv_image(image) = (rmsg(image) < rpath%tol_rmsg) .and. (maxg(image) < rpath%tol_maxg)
    end do

    rmsg_max = maxval(rmsg(:))
    maxg_max = maxval(maxg(:))
    conv = (rmsg_max < rpath%tol_rmsg) .and. (maxg_max < rpath%tol_maxg)

    dmax = 0.0_wp
    do i = 1, rpath%nreplica
      if (abs(work(i)) > dmax) dmax = abs(work(i))
    end do

    ! Output
    !
    if (main_rank) then
      
      write(MsgOut, '(/,"Image ", &
                        "     Energy (kcal/mol)", &
                        "           Relative E.", &
                        "                DeltaE", &
                        "                  RMSG", &
                        "                  MaxG", &
                        "           Convergence")')
      write(MsgOut, '("------", &
                      "----------------------", &
                      "----------------------", &
                      "----------------------", &
                      "----------------------", &
                      "----------------------", &
                      "----------------------")')
      do i = 1, rpath%nreplica
        if (conv_image(i)) then
          stat = "yes"
        else
          stat = " no"
        end if
        write(MsgOut, '(i5,X,5F22.10,19x,a3)') i, &
                             rpath%mep_energy(i), &
       rpath%mep_energy(i) - rpath%mep_energy(1), &
                                         work(i), &
                                         rmsg(i), &
                                         maxg(i), &
                                         stat
      end do
      write(MsgOut, '("------", &
                      "----------------------", &
                      "----------------------", &
                      "----------------------", &
                      "----------------------", &
                      "----------------------", &
                      "----------------------")')

      write(MsgOut, '(3X, "Energy Conv. (Max) = ", F25.15)') dmax
      write(MsgOut, '(3X, "RMSG   Conv. (Max) = ", F25.15)') rmsg_max
      write(MsgOut, '(3X, "MAXG   Conv. (Max) = ", F25.15)') maxg_max

    end if

    ! Switch-on CINEB
    if (rpath%climbing_image .and. rmsg_max < rpath%tol_rmsg_cineb ) then
      rpath%do_cineb = .true.
    end if
    
    deallocate(work)
    deallocate(rmsg)
    deallocate(maxg)
    
    return

  end subroutine check_convergence_neb

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    trans_mass_weight_coord
  !> @brief        
  !! @authors      YA, KY
  !! @param[in]    molecule     : molecular info
  !! @param[inout] rpath        : rpath info
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine trans_mass_weight_coord(molecule, rpath)

    ! formal arguments
    type(s_molecule),    intent(in) :: molecule
    type(s_rpath),    intent(inout) :: rpath

    ! local variables
    integer                         :: repid
    integer                         :: i, ii, iatom, offset
    real(wp)                        :: massfac

    repid = my_country_no + 1

    do i = 1, rpath%mep_natoms
      iatom = rpath%mepatom_id(i)
      massfac = sqrt(molecule%mass(iatom))

      offset = (i - 1) * 3
      rpath%mep_coord(offset+1:offset+3,repid) = &
        rpath%mep_coord(offset+1:offset+3,repid) * massfac
      rpath%force(offset+1:offset+3) =  &
        rpath%force(offset+1:offset+3) / massfac
    end do

    return

  end subroutine trans_mass_weight_coord

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    backtrans_mass_weight_coord
  !> @brief        
  !! @authors      YA, KY
  !! @param[in]    molecule     : molecular info
  !! @param[inout] rpath        : rpath info
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine backtrans_mass_weight_coord(molecule, rpath)

    ! formal arguments
    type(s_molecule),    intent(in) :: molecule
    type(s_rpath),    intent(inout) :: rpath

    ! local variables
    integer                         :: repid
    integer                         :: i, ii, iatom, offset
    real(wp)                        :: massfac

    repid = my_country_no + 1

    do i = 1, rpath%mep_natoms
      iatom = rpath%mepatom_id(i)
      massfac = sqrt(molecule%mass(iatom))

      offset = (i - 1) * 3
      rpath%mep_coord(offset+1:offset+3,repid) = &
        rpath%mep_coord(offset+1:offset+3,repid) / massfac
      rpath%force(offset+1:offset+3) = &
        rpath%force(offset+1:offset+3) * massfac
    end do

    return

  end subroutine backtrans_mass_weight_coord

!YA_Bgn
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    run_rpath_string
  !> @brief        string method
  !! @authors      TM, YK, YA, KY
  !! @param[inout] output      : output information
  !! @param[inout] molecule    : molecule information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] dynvars     : dynamics variables information
  !! @param[inout] minimize    : minimize information
  !! @param[inout] dynamics    : dynamics information
  !! @param[inout] pairlist    : pairlist information
  !! @param[inout] boundary    : boundary conditions information
  !! @param[inout] rpath       : RPATH information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine run_rpath_string(output, molecule, enefunc, dynvars, minimize, &
                           dynamics, pairlist, boundary, rpath, conv, niter)

    ! formal arguments
    type(s_output),           intent(inout) :: output
    type(s_molecule),         intent(inout) :: molecule
    type(s_enefunc),  target, intent(inout) :: enefunc
    type(s_dynvars),  target, intent(inout) :: dynvars
    type(s_minimize), target, intent(inout) :: minimize
    type(s_dynamics),         intent(inout) :: dynamics
    type(s_pairlist),         intent(inout) :: pairlist
    type(s_boundary),         intent(inout) :: boundary
    type(s_rpath),            intent(inout) :: rpath
    logical,                    intent(out) :: conv
    integer,                    intent(out) :: niter

    ! local variables
    integer                :: i, j, k, n
    integer                :: iloop_start, iloop_end
    integer                :: replicaid
    integer                :: ii, jj, iatom, id

    integer                :: natom_micro
    integer, pointer       :: optatom_micro_id(:)
    real(wp)               :: energy0corr, e0
    real(wp), allocatable  :: coord0(:,:), force0corr(:,:)

    real(wp)    , pointer  :: coord(:,:), force(:,:)
    type(s_qmmm), pointer  :: qmmm

    character(256)         :: folder, basename
    character(5)           :: num
    logical                :: savefile

    integer                :: repid

    conv  = .true.
    niter = 0

    return

  end subroutine run_rpath_string

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    run_rpath_neb
  !> @brief        NEB method
  !! @authors      YA
  !! @param[inout] output      : output information
  !! @param[inout] molecule    : molecule information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] dynvars     : dynamics variables information
  !! @param[inout] minimize    : minimize information
  !! @param[inout] dynamics    : dynamics information
  !! @param[inout] pairlist    : pairlist information
  !! @param[inout] boundary    : boundary conditions information
  !! @param[inout] rpath       : RPATH information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine run_rpath_neb(output, molecule, enefunc, dynvars, minimize, &
                           dynamics, pairlist, boundary, rpath, conv, niter)

    ! formal arguments
    type(s_output),           intent(inout) :: output
    type(s_molecule),         intent(inout) :: molecule
    type(s_enefunc),          intent(inout) :: enefunc
    type(s_dynvars),          intent(inout) :: dynvars
    type(s_minimize),         intent(inout) :: minimize
    type(s_dynamics),         intent(inout) :: dynamics
    type(s_pairlist),         intent(inout) :: pairlist
    type(s_boundary),         intent(inout) :: boundary
    type(s_rpath),            intent(inout) :: rpath
    logical,                    intent(out) :: conv
    integer,                    intent(out) :: niter

    ! local variables
    conv  = .true.
    niter = 0

    return

  end subroutine run_rpath_neb

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    neb_energy_and_force
  !> @brief        compute energy and NEB force for each image
  !! @authors      YA
  !! @param[in]    rpath        : rpath info
  !! @param[in]    molecule     : molecular info
  !! @param[in]    enefunc      : enefunc info
  !! @param[in]    pairlist     : pairlist info
  !! @param[in]    boundary     : boundary info
  !! @param[in]    nonb_ene     : flag for calculate nonbonded energy
  !! @param[inout] dynvars      : dynvars info
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine neb_energy_and_force(ncount, rpath, molecule, enefunc, pairlist, boundary, &
                                  nonb_ene, dynvars, update)

    ! formal arguments
    integer,             intent(in) :: ncount
    type(s_rpath),    intent(inout) :: rpath
    type(s_molecule),    intent(inout) :: molecule
    type(s_enefunc),  intent(inout) :: enefunc
    type(s_pairlist),    intent(in) :: pairlist
    type(s_boundary),    intent(in) :: boundary
    logical,             intent(in) :: nonb_ene
    type(s_dynvars),  intent(inout) :: dynvars
    logical,             intent(in) :: update

    ! local variables
    character(256)                  :: folder, basename
    character(5)                    :: num
    logical                         :: savefile
    integer                         :: repid
    integer                         :: dimno, iatom, i, ii
    real(wp)                        :: x, y, x2, y2, r, t

    ! qmmm settings
    !
    if(enefunc%qmmm%do_qmmm) then
      !folder   = 'qmmm_mep'
      !write(num,'(i5.5)') ncount
      !basename = 'mep'//num
      !savefile = (mod(ncount,rpath%qmsave_period) == 0)
      !call set_runtime_qmmm(enefunc%qmmm, folder, basename, savefile)

      enefunc%qmmm%qm_classical = .FALSE.
      if(rpath%opt_micro) enefunc%qmmm%qm_get_esp = .TRUE.
    end if

    ! compute energy and force
    !
    call compute_energy(molecule, enefunc, pairlist, boundary, .true., &
                        .false.,           &
                        dynvars%coord,     &
                        dynvars%energy,    &
                        dynvars%temporary, &
                        dynvars%force,     &
                        dynvars%virial,    &
                        dynvars%virial_extern)
!TMP_YA_Bgn
!TMP    x = dynvars%coord(1,1)
!TMP    y = dynvars%coord(2,1)
!TMP    t = x*x + y*y
!TMP    dynvars%energy%total = (1.0_wp - t)*(1.0_wp - t) + y*y / t
!TMP    dynvars%force(:,:) = 0.0_wp
!TMP    dynvars%force(1,1) = -(-4.0_wp*x*(1.0_wp - t) - y*y * 2.0_wp * x/(t*t))
!TMP    dynvars%force(2,1) = -(-4.0_wp*y*(1.0_wp - t) - y*y * 2.0_wp * y/(t*t) + 2*y/t)
!TMP_YA_End

    if (.not. replica_main_rank) return

    rpath%energy = dynvars%energy%total
    ii = 0
    do i = 1, rpath%mep_natoms
      iatom = rpath%mepatom_id(i)
      rpath%force(ii+1:ii+3) = dynvars%force(1:3,iatom)
      ii = ii + 3
    end do

    ! NEB force
    !
    call neb_force(rpath, update)

    return

  end subroutine neb_energy_and_force

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    neb_force
  !> @brief        compute NEB force
  !! @authors      YA
  !! @param[in]    rpath        : rpath info
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine neb_force(rpath, update)

    ! formal arguments
    type(s_rpath), target,    intent(inout) :: rpath
    logical,                     intent(in) :: update

    ! local variables
    logical           :: isCI(rpath%nreplica)
    integer           :: natom, image, repid
    real(wp), parameter :: zero = 0.0_wp
    real(wp), pointer :: mep_energy(:), mep_coord(:,:), mep_force(:,:)
    real(wp), allocatable :: tmp1(:), tmp2(:), tau(:,:), work(:)
    real(wp)          :: norm, norm1, norm2, deltaE_f, deltaE_b, fac1, fac2, fac
    real(wp)          :: ddot


    natom = rpath%mep_natoms
    repid = my_country_no + 1

    ! Pointers
    !
    mep_energy => rpath%mep_energy
    mep_coord  => rpath%mep_coord
    mep_force  => rpath%mep_force

    ! collect image information
    if (update) then
      allocate(work(3*natom))
      work(:) = mep_coord(:,repid)
      call mpi_allgather(rpath%energy, 1, mpi_wp_real, mep_energy, 1, mpi_wp_real, &
        mpi_comm_airplane, ierror)
      call mpi_allgather(work, 3*natom, mpi_wp_real, mep_coord, 3*natom, &
        mpi_wp_real, mpi_comm_airplane, ierror)
      call mpi_allgather(rpath%force, 3*natom, mpi_wp_real, mep_force, 3*natom, &
        mpi_wp_real, mpi_comm_airplane, ierror)
      deallocate(work)
    end if
    
    ! identify climbing image
    isCI(:) = .FALSE.
    !if (rpath%climbing_image) then
    if (rpath%do_cineb) then
      do image = 2, rpath%nreplica - 1
        if (mep_energy(image) > mep_energy(image-1) .and. &
          mep_energy(image) > mep_energy(image+1)) then
          isCI(image) = .TRUE.
          if (main_rank) then
            write(*, '(" Climbing image: ",i0)') image
          end if
        end if
      end do
    end if

    ! tangential vector
    allocate(tau(3*natom,rpath%nreplica))
    allocate(tmp1(3*natom))
    allocate(tmp2(3*natom))

    tau(:,:) = zero
    do image = 2, rpath%nreplica - 1
      deltaE_f = mep_energy(image+1) - mep_energy(image)
      deltaE_b = mep_energy(image) - mep_energy(image-1)
      tmp1(:) = mep_coord(:,image+1) - mep_coord(:,image)
      tmp2(:) = mep_coord(:,image) - mep_coord(:,image-1)
      if (deltaE_b > zero .and. deltaE_f > zero) then
        tau(:,image) = tmp1(:)
      else if (deltaE_b < zero .and. deltaE_f < zero) then
        tau(:,image) = tmp2(:)
      else
        fac1 = max(abs(deltaE_b),abs(deltaE_f))
        fac2 = min(abs(deltaE_b),abs(deltaE_f))
        if (mep_energy(image+1) > mep_energy(image-1)) then
          tau(:,image) = fac1 * tmp1(:) + fac2 * tmp2(:)
        else
          tau(:,image) = fac2 * tmp1(:) + fac1 * tmp2(:)
        end if
      end if
      ! normalization
      norm = ddot(natom*3, tau(1,image), 1, tau(1,image), 1)
      tau(:,image) = tau(:,image) / sqrt(norm)
    end do

    ! NEB force: perpeudicular component
    do image = 2, rpath%nreplica - 1
      norm = ddot(natom*3, tau(1,image), 1, rpath%mep_force(1,image),1)
      if (isCI(image)) then
        mep_force(:,image) = mep_force(:,image) - 2.0_wp * norm * tau(:,image)
      else
        mep_force(:,image) = mep_force(:,image) - norm * tau(:,image)
      end if
    end do

    ! NEB force: spring force
    do image = 2, rpath%nreplica - 1
      if (isCI(image)) cycle
      tmp1(:) = mep_coord(:,image+1) - mep_coord(:,image)
      tmp2(:) = mep_coord(:,image) - mep_coord(:,image-1)
      norm1 = ddot(natom*3, tmp1, 1, tmp1, 1)
      norm2 = ddot(natom*3, tmp2, 1, tmp2, 1)
      fac = rpath%k_spring * (sqrt(norm1) - sqrt(norm2))
      mep_force(:,image) = mep_force(:,image) + fac * tau(:,image)      
    end do

    deallocate(tmp1)
    deallocate(tmp2)
    deallocate(tau)

    ! clear force for terminals
    !
    if (rpath%fix_terminal) then
      mep_force(:,1) = zero
      mep_force(:,rpath%nreplica) = zero
    end if

    return

  end subroutine neb_force

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    gradient_correction
  !> @brief        energy and gradient correction for micro iteration
  !! @authors      YA
  !! @param[in]    rpath        : rpath info
  !! @param[in]    molecule     : molecular info
  !! @param[in]    enefunc      : enefunc info
  !! @param[in]    pairlist     : pairlist info
  !! @param[in]    boundary     : boundary info
  !! @param[in]    nonb_ene     : flag for calculate nonbonded energy
  !! @param[inout] dynvars      : dynvars info
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine gradient_correction(minimize, molecule, enefunc, pairlist, boundary, dynvars, &
                                 rpath, energy0corr, coord0, force0corr)

    ! formal arguments
    type(s_minimize), target,   intent(in) :: minimize
    type(s_molecule),        intent(inout) :: molecule
    type(s_enefunc), target, intent(inout) :: enefunc
    type(s_pairlist),           intent(in) :: pairlist
    type(s_boundary),           intent(in) :: boundary
    type(s_dynvars), target, intent(inout) :: dynvars
    type(s_rpath),              intent(inout) :: rpath
    real(wp),                   intent(out) :: energy0corr, coord0(3,*), force0corr(3,*)

    ! local variables
    integer               :: i, repid, natom_micro
    integer, pointer      :: optatom_micro_id(:)
    real(wp), pointer     :: coord(:,:), force(:,:)
    type(s_qmmm), pointer :: qmmm


    natom_micro = minimize%num_optatoms_micro

    ! Pointers
    !
    coord => dynvars%coord
    force => dynvars%force
    qmmm => enefunc%qmmm
    optatom_micro_id => minimize%optatom_micro_id

    ! energy and gradient correction terms for micro_iteration
    !
    energy0corr = dynvars%energy%total
    do i = 1, natom_micro
      coord0(1:3,i)     = coord(1:3,optatom_micro_id(i))
      force0corr(1:3,i) = force(1:3,optatom_micro_id(i))
    end do

    qmmm%qm_classical = .TRUE.
    call compute_energy(molecule, enefunc, pairlist, boundary, .true., &
                        enefunc%nonb_limiter,                     &
                        coord, dynvars%energy, dynvars%temporary, &
                        force, dynvars%virial, dynvars%virial_extern)

    energy0corr = energy0corr - dynvars%energy%total
    do i = 1, natom_micro
      force0corr(1:3,i) = force0corr(1:3,i) - force(1:3,optatom_micro_id(i))
    end do

    ! QM internal energy
    !
    repid = my_country_no + 1
    rpath%qm_energy(repid) = energy0corr

    return

  end subroutine gradient_correction
!YA_End

end module at_rpath_mep_mod

