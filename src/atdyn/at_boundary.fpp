!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_boundary_mod
!> @brief   utilities for boundary conditions
!! @authors Takaharu Mori (TM), Takashi Imai (TI), Jaewoon Jung (JJ), 
!!          Norio Takase (NT), Motoshi Kamiya (MK), Kiyoshi Yagi (KY),
!!          Yuji Sugita (YS)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_boundary_mod

  use at_boundary_str_mod
  use molecules_str_mod
  use molecules_mod
  use fileio_rst_mod
  use fileio_spot_mod
  use fileio_control_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  use select_mod
  use select_atoms_mod
  use select_atoms_str_mod
#ifdef MPI
  use mpi
#endif

  implicit none
  private

  real(wp), parameter   :: TableWaterMargin = 2.0_wp
  real(wp), parameter   :: FakeDefault      = 999999_wp

  ! structures
  type, public :: s_pbc_info
    real(wp)            :: box_size_x       = FakeDefault
    real(wp)            :: box_size_y       = FakeDefault
    real(wp)            :: box_size_z       = FakeDefault
    real(wp)            :: pairlist_grid    = 3.0_wp
    logical             :: wrap_all         = .false.
  end type s_pbc_info

  type, public :: s_spot_info
    ! spherical potential
    logical                         :: spherical_pot  = .false.
    real(wp)                        :: const          = 10.0_wp
    integer                         :: exponent       = 2
    logical                         :: read_spot      = .false.
    real(wp)                        :: mod_ratio      = 1.0_wp
    integer                         :: nindex         = 0
    integer                         :: nfunctions     = 0
    character(MaxLine), allocatable :: center(:) 
    real(wp), allocatable           :: radius(:)
    logical                         :: fixatom        = .true.
    real(wp)                        :: fix_layer      = 1.0_wp
    character(MaxLine)              :: nospot_select_index = ''
    logical                         :: restart        = .true.
  end type s_spot_info

  type, public :: s_boundary_info
    integer             :: type             = BoundaryTypePBC
    type(s_pbc_info)    :: pbc_info
    type(s_spot_info)   :: spot_info
    real(wp)            :: origin_x         = 0.0_wp
    real(wp)            :: origin_y         = 0.0_wp
    real(wp)            :: origin_z         = 0.0_wp
  end type s_boundary_info

  ! subroutines
  public  :: show_ctrl_boundary
  public  :: read_ctrl_boundary
  public  :: setup_boundary
  public  :: update_boundary
  private :: update_boundary_pbc
  public  :: update_spot_atomlist
  public  :: wrap_molecules
  private :: wrap_all
  public  :: clear_fixatm_component

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_boundary
  !> @brief        show BOUNDARY section usage
  !! @authors      NT, TM, KY
  !! @param[in]    show_all : show all usage or not
  !! @param[in]    run_mode : run mode string : "md", "min", "remd", "rpath"
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine show_ctrl_boundary(show_all, run_mode)

    ! formal arguments
    logical,                 intent(in)    :: show_all
    character(*),            intent(in)    :: run_mode


    if (show_all) then

      select case (run_mode)

      case ('md', 'min', 'vib', 'remd', 'rpath')

        write(MsgOut,'(A)') '[BOUNDARY]'
        write(MsgOut,'(A)') 'type          = PBC       # [NOBC,PBC]'
        write(MsgOut,'(A)') 'box_size_x    = 999999    # box size (x)     in [PBC]'
        write(MsgOut,'(A)') 'box_size_y    = 999999    # box size (y)     in [PBC]'
        write(MsgOut,'(A)') 'box_size_z    = 999999    # box size (z)     in [PBC]'
        write(MsgOut,'(A)') ' '


      end select

    else

      select case (run_mode)

      case ('md', 'min', 'vib', 'remd', 'rpath')

        write(MsgOut,'(A)') '[BOUNDARY]'
        write(MsgOut,'(A)') 'type          = PBC       # [NOBC,PBC]'
        write(MsgOut,'(A)') 'box_size_x    = 999999    # box size (x) in [PBC]'
        write(MsgOut,'(A)') 'box_size_y    = 999999    # box size (y) in [PBC]'
        write(MsgOut,'(A)') 'box_size_z    = 999999    # box size (z) in [PBC]'
        write(MsgOut,'(A)') ' '


      end select

    end if


    return

  end subroutine show_ctrl_boundary
  
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_boundary
  !> @brief        read BOUNDARY section in the control file
  !! @authors      YS, TI, JJ, TM, KY
  !! @param[in]    handle     : unit number
  !! @param[out]   bound_info : BOUNDARY section control parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ctrl_boundary(handle, bound_info)

    ! parameters
    character(*),            parameter     :: Section = 'Boundary'

    ! formal arguments
    integer,                        intent(in)    :: handle
    type(s_boundary_info), target,  intent(inout) :: bound_info

    ! local variables
    integer                    :: i, nfunc, ierr
    character(20)              :: key1, key2
    type(s_spot_info), pointer :: spot_info

    ! use pointers
    !
    spot_info => bound_info%spot_info

    ! read parameters from control file
    ! 
    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_type(handle, Section, 'type',          &
                            bound_info%type, BoundaryTypeTypes)

    ! read parameters for PBC from control file
    !
    call read_ctrlfile_real(handle, Section, 'box_size_x',  &
                            bound_info%pbc_info%box_size_x)
    call read_ctrlfile_real(handle, Section, 'box_size_y',  &
                            bound_info%pbc_info%box_size_y)
    call read_ctrlfile_real(handle, Section, 'box_size_z',  &
                            bound_info%pbc_info%box_size_z)
    call read_ctrlfile_logical(handle, Section, 'wrapall',  &
                            bound_info%pbc_info%wrap_all)
    call read_ctrlfile_real(handle, Section, 'pairlist_grid',  &
                              bound_info%pbc_info%pairlist_grid)

    call end_ctrlfile_section(handle)


    ! write parameters to MsgOut
    !
    if (main_rank) then

      write(MsgOut,'(A)') 'Read_Ctrl_Boundary> Parameters of Boundary Condition'
      write(MsgOut,'(A20,A10)') &
            '  type            = ', BoundaryTypeTypes(bound_info%type)

      if (abs(bound_info%pbc_info%box_size_x - FakeDefault) > EPS .and. &
          abs(bound_info%pbc_info%box_size_y - FakeDefault) > EPS .and. &
          abs(bound_info%pbc_info%box_size_z - FakeDefault) > EPS) &
      write(MsgOut,'(A20,3F10.3)')                                  &
            '  box_size(x,y,z) = ', bound_info%pbc_info%box_size_x, &
                                    bound_info%pbc_info%box_size_y, &
                                    bound_info%pbc_info%box_size_z
      write(MsgOut,'(A20,F10.3)')                                   &
            '  pairlist_grid   = ', bound_info%pbc_info%pairlist_grid

      write(MsgOut,'(A)') ' '

    end if

    return

  end subroutine read_ctrl_boundary

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_boundary
  !> @brief        set essential variables for boundary condition
  !! @authors      YS, TM, JJ, TI, KY
  !! @param[in]    bound_info  : BOUNDARY section control parameters information
  !! @param[in]    use_table   : flag for use table or not
  !! @param[in]    pairlistdist: pair-list distance
  !! @param[in]    sel_info    : SELECTOR section in control parameters
  !! @param[in]    molecule    : molecule
  !! @param[in]    rst         : restart file information
  !! @param[in]    spot        : spherical potential information
  !! @param[out]   boundary    : boundary conditions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_boundary(bound_info, use_table, pairlistdist, sel_info, &
                            molecule, rst, spot, boundary)

    ! formal arguments
    type(s_boundary_info),   intent(in)    :: bound_info
    logical,                 intent(in)    :: use_table
    real(wp),                intent(in)    :: pairlistdist
    type(s_sel_info),        intent(in)    :: sel_info
    type(s_molecule),        intent(inout) :: molecule
    type(s_rst),             intent(in)    :: rst
    type(s_spot),            intent(in)    :: spot
    type(s_boundary),        intent(inout) :: boundary

    ! local variables
    integer    :: nf, i, iatom

    ! initialize 
    !
    call init_boundary(boundary)

    ! allocate fixatm
    allocate(boundary%fixatm(molecule%num_atoms))
    boundary%fixatm = .false.

    ! setup variables
    !
    select case (bound_info%type)
    

    case (BoundaryTypePBC)

      boundary%type     = bound_info%type
      boundary%origin_x = bound_info%origin_x
      boundary%origin_y = bound_info%origin_y
      boundary%origin_z = bound_info%origin_z
      boundary%wrap_all = bound_info%pbc_info%wrap_all

      if (rst%rstfile_type == RstfileTypeUndef) then
        if (abs(bound_info%pbc_info%box_size_x - FakeDefault) < EPS) &
          call error_msg('Setup_Boundary> box_size_x is not specified in ctrl')
        if (abs(bound_info%pbc_info%box_size_y - FakeDefault) < EPS) &
          call error_msg('Setup_Boundary> box_size_y is not specified in ctrl')
        if (abs(bound_info%pbc_info%box_size_z - FakeDefault) < EPS) &
          call error_msg('Setup_Boundary> box_size_z is not specified in ctrl')

        boundary%box_size_x = bound_info%pbc_info%box_size_x
        boundary%box_size_y = bound_info%pbc_info%box_size_y
        boundary%box_size_z = bound_info%pbc_info%box_size_z
      else
        boundary%box_size_x = rst%box_size_x
        boundary%box_size_y = rst%box_size_y
        boundary%box_size_z = rst%box_size_z
      end if
      boundary%box_size_x_ref = boundary%box_size_x
      boundary%box_size_y_ref = boundary%box_size_y
      boundary%box_size_z_ref = boundary%box_size_z

      boundary%pairlist_grid = bound_info%pbc_info%pairlist_grid

    case (BoundaryTypeNOBC)

      boundary%type     = bound_info%type
      boundary%origin_x = bound_info%origin_x
      boundary%origin_y = bound_info%origin_y
      boundary%origin_z = bound_info%origin_z

    end select

    ! write setup info
    !
    if (main_rank) then
      write(MsgOut,'(A)') &
            'Setup_Boundary> Setup Variables for Boundary Condition'

      if (boundary%type == BoundaryTypePBC) then
        write(MsgOut,'(A20,3F10.3)')                       &
              '  box_size(x,y,z) = ', boundary%box_size_x, &
                                      boundary%box_size_y, &
                                      boundary%box_size_z
      end if

      write(MsgOut,'(A20,3F10.3)')                         &
              '  origin(x,y,z)   = ', boundary%origin_x,   &
                                      boundary%origin_y,   &
                                      boundary%origin_z
      write(MsgOut,'(A)') ' '
      write(MsgOut,'(A)') ' '

    end if

    ! update DOF
    if(boundary%num_fixatm > 0) then
      call update_num_deg_freedom('After setup of fixed atom', &
                                  -3*boundary%num_fixatm,      &
                                  molecule%num_deg_freedom)
    end if

    ! update boundary
    !
    call update_boundary(use_table, pairlistdist, boundary)

    return

  end subroutine setup_boundary

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_boundary
  !> @brief        update variables for boundary condition
  !! @authors      TM
  !! @param[in]    use_table    : flag for use table or not
  !! @param[in]    pairlistdist : pair list distance
  !! @param[inout] boundary     : information of boundary condition
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_boundary(use_table, pairlistdist, boundary)

    ! formal arguments
    logical,                 intent(in)    :: use_table
    real(wp),                intent(in)    :: pairlistdist
    type(s_boundary),        intent(inout) :: boundary


    select case (boundary%type)

    case (BoundaryTypeNOBC)

      ! do nothing

    case (BoundaryTypePBC)

      call update_boundary_pbc(use_table, pairlistdist, boundary)

    end select

    return

  end subroutine update_boundary

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_boundary_pbc
  !> @brief        update variables for boundary condition
  !! @authors      TI, JJ, TM, MK
  !! @param[in]    use_table    : flag for use table or not
  !! @param[in]    pairlistdist : pair list distance
  !! @param[inout] boundary     : information of boundary condition
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_boundary_pbc(use_table, pairlistdist, boundary)

    ! formal arguments
    logical,                 intent(in)    :: use_table
    real(wp),                intent(in)    :: pairlistdist
    type(s_boundary),        intent(inout) :: boundary

    ! local variables
    real(wp)                 :: bsize_x, bsize_y, bsize_z
    real(wp)                 :: csize_x, csize_y, csize_z
    integer                  :: i, j, k
    integer                  :: inb, jnb, knb, inbs, jnbs, knbs, lc, lcnb, ln
    integer                  :: ncell_x, ncell_y, ncell_z, ncell
    integer                  :: ncell_x_check, ncell_y_check, ncell_z_check

    ! for small cell division
    real(wp)                 :: pairdist_w_margin, pairdist_w_margin2
    real(wp)                 :: pairdist_check
    real(wp)                 :: ci, cj, ck
    integer                  :: cellpair_max_x, cellpair_max_y, cellpair_max_z
    integer                  :: ik, jk, kk


    boundary%use_cell_linked_list = .true.
    bsize_x = boundary%box_size_x
    bsize_y = boundary%box_size_y
    bsize_z = boundary%box_size_z

    if (use_table) then

      ncell_x = int(bsize_x/boundary%pairlist_grid)
      ncell_y = int(bsize_y/boundary%pairlist_grid)
      ncell_z = int(bsize_z/boundary%pairlist_grid)

      pairdist_check = max(boundary%pairlist_grid,pairlistdist)

      ncell_x_check = int(bsize_x/pairdist_check)
      ncell_y_check = int(bsize_y/pairdist_check)
      ncell_z_check = int(bsize_z/pairdist_check)

      if (ncell_x_check < 3 .or. ncell_y_check < 3 .or. ncell_z_check < 3) &
        call error_msg('Update_Boundary_Pbc> too small boxsize/pairdist.'//&
                       ' larger boxsize or shorter pairdist should be used (see "Chapter: Trouble shooting" in the user manual).')

    else

      ncell_x = int(bsize_x/pairlistdist)
      ncell_y = int(bsize_y/pairlistdist)
      ncell_z = int(bsize_z/pairlistdist)

      if (ncell_x < 3 .or. ncell_y < 3 .or. ncell_z < 3) then
        if (pairlistdist < 0.5_wp*bsize_x .and. &
            pairlistdist < 0.5_wp*bsize_y .and. &
            pairlistdist < 0.5_wp*bsize_z) then
          boundary%use_cell_linked_list = .false.
          return
        else
          call error_msg('Update_Boundary_Pbc> too small boxsize/pairdist.'//&
                         ' larger boxsize or shorter pairdist should be used (see "Chapter: Trouble shooting" in the user manual).')
        end if
      end if

    end if

    csize_x = bsize_x/real(ncell_x,wp)
    csize_y = bsize_y/real(ncell_y,wp)
    csize_z = bsize_z/real(ncell_z,wp)
    ncell   = ncell_x*ncell_y*ncell_z

    ! MK (cellpairs; TableWaterMargin for water pairs)
    pairdist_w_margin = pairlistdist + TableWaterMargin
    pairdist_w_margin2 = pairdist_w_margin * pairdist_w_margin
    cellpair_max_x = pairdist_w_margin / csize_x + 1
    cellpair_max_y = pairdist_w_margin / csize_y + 1
    cellpair_max_z = pairdist_w_margin / csize_z + 1
    boundary%num_neighbor_cells = 0

    do k = -cellpair_max_z, cellpair_max_z

      ck = csize_z * real( max( 0, abs(k) - 1 ), wp )
      ck = ck * ck

      do j = -cellpair_max_y, cellpair_max_y

        cj = csize_y * real( max( 0, abs(j) - 1 ), wp )
        cj = cj * cj

        do i = -cellpair_max_x, cellpair_max_x

          ci = csize_x * real( max( 0, abs(i) - 1 ), wp )
          ci = ci * ci

          if ( (ci+cj+ck) < pairdist_w_margin2 ) then
            boundary%num_neighbor_cells = boundary%num_neighbor_cells + 1
          endif

        end do
      end do
    end do

    boundary%num_cells_x = ncell_x
    boundary%num_cells_y = ncell_y
    boundary%num_cells_z = ncell_z
    boundary%cell_size_x = csize_x
    boundary%cell_size_y = csize_y
    boundary%cell_size_z = csize_z
    boundary%num_cells   = ncell

    ! prepare cell neighbor list
    !
    call alloc_boundary(boundary, BoundaryCells, ncell)

    ln = 0
    do k = -cellpair_max_z, cellpair_max_z

      ck = csize_z * real( max( 0, abs(k) - 1 ), wp )
      ck = ck * ck

      do j = -cellpair_max_y, cellpair_max_y

        cj = csize_y * real( max( 0, abs(j) - 1 ), wp )
        cj = cj * cj

        do i = -cellpair_max_x, cellpair_max_x

          ci = csize_x * real( max( 0, abs(i) - 1 ), wp )
          ci = ci * ci

          if ( (ci+cj+ck) < pairdist_w_margin2 ) then
            ln = ln + 1
            boundary%neighbor_cell_common_x(ln) = i
            boundary%neighbor_cell_common_y(ln) = j
            boundary%neighbor_cell_common_z(ln) = k
          endif

        end do
      end do
    end do

    do k = 0, ncell_z-1
      do j = 0, ncell_y-1
        do i = 0, ncell_x-1
              
          lc = 1 + i + j*ncell_x + k*ncell_x*ncell_y

          do ln = 1, boundary%num_neighbor_cells

            inb = i + boundary%neighbor_cell_common_x(ln)
            if (inb < 0) then
              inbs = ncell_x + inb
            else if (inb >= ncell_x) then
              inbs = inb - ncell_x
            else
              inbs = inb
            end if

            jnb = j + boundary%neighbor_cell_common_y(ln)
            if (jnb < 0) then
              jnbs = ncell_y + jnb
            else if (jnb >= ncell_y) then
              jnbs = jnb - ncell_y
            else
              jnbs = jnb
            end if

            knb = k + boundary%neighbor_cell_common_z(ln)
            if (knb < 0) then
              knbs = ncell_z + knb
            else if (knb >= ncell_z) then
              knbs = knb - ncell_z
            else
              knbs = knb
            end if

            lcnb = 1 + inbs + jnbs*ncell_x + knbs*ncell_x*ncell_y
            boundary%neighbor_cells(ln,lc) = lcnb

          end do

        end do
      end do
    end do

    return

  end subroutine update_boundary_pbc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_spot_atomlist
  !> @brief        update list of atoms for spherical potential
  !!               This routine is currently not used.
  !! @authors      KY
  !! @param[in]    natom    : number of atoms
  !! @param[in]    coord    : coordinates
  !! @param[inout] boundary : boundary conditions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_spot_atomlist(natom, coord, boundary)

    ! formal arguments
    integer,                 intent(in)    :: natom
    real(wp),                intent(in)    :: coord(:,:)
    type(s_boundary),        intent(inout) :: boundary

    ! local variables
    integer  :: nfunc, nc
    integer  :: i, j, n
    integer  :: id, my_id
    real(wp) :: rn, rmin, din(3), rin
#ifdef OMP
    integer      :: omp_get_thread_num, omp_get_max_threads
#endif

#ifdef OMP
    nthread = omp_get_max_threads()
#else
    nthread = 1
#endif

    nfunc = boundary%nfunctions
    nc    = 0

    !$omp parallel                     &
    !$omp private(i, j, n,             &
    !$omp         rn, rmin, din, rin,  &
    !$omp         id, my_id)           &
    !$omp reduction(+:nc) 

#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    my_id = my_city_rank * nthread + id

    do i = 1, natom
      if (mod(i-1,nproc_city*nthread) /= my_id) cycle

      ! skip fix atoms
      if (boundary%fixatm(i)) then
        boundary%atomlist(i) = -1
        cycle
      end if

      ! find the nearest center
      rmin = -1.0_wp
      do n = 1, nfunc
        rn = boundary%radius(n)

        din = coord(:,i) - boundary%center(:,n)
        rin = sqrt(din(1)*din(1) + din(2)*din(2) + din(3)*din(3))

        if (rin < rn) then
          boundary%atomlist(i) = -1
          exit

        else
          if (rmin < 0.0_wp .or. (rin - rn) < rmin) then
            rmin  = rin - rn
            boundary%atomlist(i) = n
          end if

        end if

      end do

    end do
    !$omp end parallel

  end subroutine update_spot_atomlist

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    wrap molecules
  !> @brief        wrap molecules
  !! @authors      TM
  !! @param[in]    molecule : molecule information
  !! @param[in]    boundary : boundary conditions information
  !! @param[inout] coord    : coordinates
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine wrap_molecules(molecule, boundary, coord)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_boundary),        intent(in)    :: boundary
    real(wp),                intent(inout) :: coord(:,:)


    if (boundary%wrap_all) then
      call wrap_all(molecule, boundary, coord)
    else
      ! do nothing
    end if

    return

  end subroutine wrap_molecules

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    wrap_all
  !> @brief        wrap molecules
  !! @authors      TM
  !! @param[in]    molecule : molecule information
  !! @param[in]    boundary : boundary conditions information
  !! @param[inout] coord    : coordinates
  !! @note         both coord_ref and coord are wrapped for restart
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine wrap_all(molecule, boundary, coord)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_boundary),        intent(in)    :: boundary
    real(wp),                intent(inout) :: coord(:,:)

    ! local variables
    real(wp)                 :: com(3), box(3), ori(3)
    integer                  :: i, j, initial, final


    box(1) = boundary%box_size_x
    box(2) = boundary%box_size_y
    box(3) = boundary%box_size_z
    ori(1) = boundary%origin_x
    ori(2) = boundary%origin_y
    ori(3) = boundary%origin_z

    do i = 1, molecule%num_molecules

      initial  = molecule%molecule_atom_no(i)
      if (i /= molecule%num_molecules) then
        final = molecule%molecule_atom_no(i+1) - 1
      else
        final = molecule%num_atoms
      end if

      ! compute center of mass of a molecule
      !
      com(1:3) = 0.0_wp
      do j = initial, final
        coord(1:3,j) = coord(1:3,j) - ori(1:3)
        com(1:3) = com(1:3) + coord(1:3,j)*molecule%mass(j)
      end do
      com(1:3) = com(1:3)/molecule%molecule_mass(i)

      ! move molecule into the unit cell
      !
      do j = initial, final
        coord(1:3,j) = coord(1:3,j) - box(1:3)*nint(com(1:3)/box(1:3))
        coord(1:3,j) = coord(1:3,j) + ori(1:3)
      end do

    end do

    return

  end subroutine wrap_all

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    clear_component
  !> @brief        clear to zero the specified component
  !! @authors      KY
  !! @param[in]    boundary : boundary conditions information
  !! @param[in]    natom    : number of atoms
  !! @param[inout] comp     : velocity or force
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine clear_fixatm_component(boundary, natom, comp)

    ! formal arguments
    type(s_boundary),    intent(in)    :: boundary
    integer,             intent(in)    :: natom
    real(wp),            intent(inout) :: comp(3,natom)

    ! local variable
    integer  :: i

    if (boundary%num_fixatm == 0) return

    !$omp parallel 
    !$omp do 
    do i = 1, natom
      if (boundary%fixatm(i)) comp(:,i) = 0.0_wp
    end do
    !$omp end do
    !$omp end parallel

  end subroutine clear_fixatm_component

end module at_boundary_mod
