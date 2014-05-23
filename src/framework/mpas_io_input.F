! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module mpas_io_input

   use mpas_grid_types
   use mpas_dmpar
   use mpas_block_decomp
   use mpas_block_creator
   use mpas_sort
   use mpas_configure
   use mpas_timekeeping
   use mpas_io_streams
   use mpas_io_units

   integer, parameter :: STREAM_INPUT=1, STREAM_SFC=2, STREAM_RESTART=3

   type io_input_object
      character (len=StrKIND) :: filename
      integer :: rd_ncid
      integer :: stream

      integer :: time

      type (MPAS_Stream_type) :: io_stream

   end type io_input_object

   integer :: readCellStart, readCellEnd, nReadCells
   integer :: readEdgeStart, readEdgeEnd, nReadEdges
   integer :: readVertexStart, readVertexEnd, nReadVertices

   contains

!   subroutine mpas_input_state_for_domain(domain)!{{{
!
!       STEVE_edit: this routine runs mpas_input_state_for_domain_iter
!       n_mgslevel times, this and mpas_input_state_for_domain_iter replaces
!       mpas_input_state_for_domain 
!
   subroutine mpas_input_state_for_domain(domain)!{{{

      implicit none

      type (domain_type), pointer :: domain

      integer :: mgslevel
      integer :: n_mgslevel

      n_mgslevel = config_num_mgslevels

      ! Loop over the the number of mgslevels and run input state

      do mgslevel=1,n_mgslevel
        call mpas_input_state_for_domain_iter(domain,mgslevel)
      end do
      !STEVE_edit routine added to test the multigrid
      !call test_mgs(domain)

!	write(0,*)'weights =', domain % blocklist % mgs_vars % restrictStencil % array(1,1:5)
!	write(0,*)'weights =', domain % blocklist % mgs_vars % restrictWeights % array(1,1:5)
!	if (n_mgslevel == 2) then
!	  write(0,*)'mgslevel = 2 and sten =', domain % blocklist % coarser % mgs_vars % restrictStencil % array(1,1:5)
!	  write(0,*)'weights =', domain % blocklist % coarser % mgs_vars % restrictWeights % array(1,1:5)
!	end if
! check resict is being read in
      
   end subroutine mpas_input_state_for_domain!}}}

! comented out as it doent work for init (as iv not declared mgs_vars in Reg file)
!!Steve testing the multigrid hierarchy
!   subroutine test_mgs(domain)!{{{
!       implicit none
!
!       type (domain_type), pointer :: domain
!
!       type (block_type), pointer :: block_ptr, block_ptr_coarser
!
!	integer :: n1, n2, nCellsSolve, nCells
!
!	write(0,*)'STEVE: testing the multigrid hierarchy'
!
!       	block_ptr => domain % blocklist
!       	! pass arrays from pressure field to the mgs_vars pressure field
!       	write(0,*)'size(helm_inc)', size(block_ptr % mgs_vars % helm_inc % array,1), size(block_ptr % mgs_vars % helm_inc % array,2)
!       	write(0,*)'size(pressure_base)', size(block_ptr % diag % pressure_base % array,1), size(block_ptr % diag % pressure_base % array,2)
!
!	write(0,*) ' 1:4 vals = ' ,block_ptr % diag % pressure_base % array(1,1:4),' - - ',block_ptr % mgs_vars % helm_inc % array(1,1:4)
!	block_ptr % mgs_vars % helm_inc % array = block_ptr % diag % pressure_base % array
!        write(0,*) ' 1:4 vals = ' ,block_ptr % diag % pressure_base % array(1,1:4),' - - ',block_ptr % mgs_vars % helm_inc % array(1,1:4)
!
!
!	! pass pressure down to lower levels
!	block_ptr_coarser => block_ptr % coarser
!	n1 = size(block_ptr_coarser % mgs_vars % helm_inc % array,1)
!	n2 = size(block_ptr_coarser % mgs_vars % helm_inc % array,2)
!
!	block_ptr_coarser % mgs_vars % helm_inc % array(1:n1,1:n2) = block_ptr % mgs_vars % helm_inc % array(1:n1,1:n2)
!
!	write(0,*) ' 1:4 vals coarse1 = ' ,block_ptr_coarser % mgs_vars % helm_inc % array(1,1:4)
!	write(0,*) ' 1:4 vals coarse2 = ' ,block_ptr % coarser % mgs_vars % helm_inc % array(1,1:4)
!
!	! check the domains are communicateing
!
!	! Somthing like this:
!	!mpas_dmpar_exch_halo_field(field , haloLayersIn)
!
!	nCells = block_ptr_coarser % mesh % nCells
!	nCellsSolve = block_ptr_coarser % mesh % nCellsSolve
!
!        write(0,*) 'Tetsing comms is working'
!
!	write(0,*) '1. owned/hallo',block_ptr_coarser % mesh % indexToCellID % array(nCellsSolve-5:nCellsSolve+5),sum(block_ptr_coarser % mesh % indexToCellID % array(nCellsSolve+1:nCells))
!	block_ptr_coarser % mesh % indexToCellID % array(nCellsSolve+1:nCells) = 0
!
!	write(0,*) '2. owned/hallo',block_ptr_coarser % mesh % indexToCellID % array(nCellsSolve-5:nCellsSolve+5),sum(block_ptr_coarser % mesh % indexToCellID % array(nCellsSolve+1:nCells))
!	call mpas_dmpar_exch_halo_field(block_ptr_coarser % mesh % indexToCellID)
!
!	write(0,*) '3. owned/hallo',block_ptr_coarser % mesh % indexToCellID % array(nCellsSolve-5:nCellsSolve+5),sum(block_ptr_coarser % mesh % indexToCellID % array(nCellsSolve+1:nCells))
!
!   end subroutine test_mgs!{{{

!  subroutine mpas_input_state_for_domain_iter(domain,mgslevel)!{{{
!
!       STEVE_edit: this routine is called from mpas_input_state_for_domain
!       n_mgslevel times. This routine is based on the original to handle mulple grids
!
   subroutine mpas_input_state_for_domain_iter(domain,mgslevel)!{{{

      implicit none
   
      type (domain_type), pointer :: domain
!STEVE_edit: add varible to define which level we are working at
      integer, intent(in) :: mgslevel

      type (block_type), pointer :: block_ptr
      type (block_type), pointer :: readingBlock
   
      integer :: i, j, k
      type (io_input_object) :: input_obj
#include "dim_decls.inc"

      character (len=StrKIND) :: c_on_a_sphere
      real (kind=RKIND) :: r_sphere_radius

      integer :: ierr
      integer, dimension(:), pointer :: readIndices
      type (MPAS_IO_Handle_type) :: inputHandle
   
      type (field1dInteger), pointer :: indexToCellIDField
      type (field1dInteger), pointer :: indexToEdgeIDField
      type (field1dInteger), pointer :: indexToVertexIDField
      type (field1dInteger), pointer :: nEdgesOnCellField
      type (field2dInteger), pointer :: cellsOnCellField
      type (field2dInteger), pointer :: edgesOnCellField
      type (field2dInteger), pointer :: verticesOnCellField
      type (field2dInteger), pointer :: cellsOnEdgeField
      type (field2dInteger), pointer :: cellsOnVertexField

      type (field1dReal), pointer :: xCellField,   yCellField,   zCellField
      type (field1dReal), pointer :: xEdgeField,   yEdgeField,   zEdgeField
      type (field1dReal), pointer :: xVertexField, yVertexField, zVertexField

      type (field1DChar) :: xtime

      type (field1dInteger), pointer :: nCellsSolveField
      type (field1dInteger), pointer :: nVerticesSolveField
      type (field1dInteger), pointer :: nEdgesSolveField

      type (field1DInteger), pointer :: indexToCellID_Block
      type (field1DInteger), pointer :: nEdgesOnCell_Block
      type (field2DInteger), pointer :: cellsOnCell_Block
      type (field2DInteger), pointer :: verticesOnCell_Block
      type (field2DInteger), pointer :: edgesOnCell_Block

      type (field1DInteger), pointer :: indexToVertexID_Block
      type (field2DInteger), pointer :: cellsOnVertex_Block

      type (field1DInteger), pointer :: indexToEdgeID_Block
      type (field2DInteger), pointer :: cellsOnEdge_Block

      type (field1DReal), pointer :: xCell, yCell, zCell
      type (field1DReal), pointer :: xEdge, yEdge, zEdge
      type (field1DReal), pointer :: xVertex, yVertex, zVertex
   
      integer, dimension(:), pointer :: local_cell_list
      integer, dimension(:), pointer :: block_id, block_start, block_count
      type (graph) :: partial_global_graph_info

      type (MPAS_Time_type) :: startTime
      character(len=StrKIND) :: timeStamp, restartTimeStamp
      character(len=StrKIND) :: filename

      integer :: nHalos
!STEVE_edit: dumChar to spec a temp char to pass into mgslevelChar
!STEVE_edit: mgslevelChar to define a char = _mgslevel
      integer :: nGrids
      character(len=StrKIND) :: dumChar
      character(len=3) :: mgslevelChar

      nHalos = config_num_halos

!STEVE_edit: make a char to use when reading the .nc file
      if (mgslevel.eq.1) then
        !for the first time there is no level defined in the nc fle so make an
        !emtpy string
        mgslevelChar=''
      elseif(mgslevel.lt.10) then
        write(mgslevelChar,'(I1)') mgslevel
        mgslevelChar='_'//trim(mgslevelChar)
      elseif(mgslevel.lt.100) then
        write(mgslevelChar,'(I2)') mgslevel
        mgslevelChar='_'//trim(mgslevelChar)
      else
        write(*,*) 'I doubt that there should be more than 100 layers'
        write(*,*) '... in mpas_input_state_for_domain_iter'
      endif


      if (config_do_restart) then
        ! this get followed by set is to ensure that the time is in standard format
        if(trim(config_start_time) == 'file') then
          open(22,file=trim(config_restart_timestamp_name),form='formatted',status='old')
          read(22,*) restartTimeStamp
          close(22)

        else
          restartTimeStamp = config_start_time
        end if

        write(stderrUnit,*) 'RestartTimeStamp ', trim(restartTimeStamp)
        call mpas_set_time(curr_time=startTime, dateTimeString=restartTimeStamp)
        call mpas_get_time(curr_time=startTime, dateTimeString=timeStamp)
        call mpas_insert_string_suffix(trim(config_restart_name), timeStamp, filename)

        input_obj % filename = trim(filename)
        input_obj % stream = STREAM_RESTART
      else
!STEVE_edit: 1st level read input file: either init or gridi otherwise read the
!grid file
        if (mgslevel.eq.1) then
          input_obj % filename = trim(config_input_name)
        else
          input_obj % filename = trim(config_mg_input_name)
        end if
        input_obj % stream = STREAM_INPUT
      end if
      inputHandle = MPAS_io_open(trim(input_obj % filename), MPAS_IO_READ, MPAS_IO_PNETCDF, ierr)
      if (ierr /= MPAS_IO_NOERR) then
        write(stderrUnit,*) ' '
        if (input_obj % stream == STREAM_RESTART) then
          write(stderrUnit,*) 'Error opening restart file ''', trim(input_obj % filename), ''''
        else if (input_obj % stream == STREAM_INPUT) then
          write(stderrUnit,*) 'Error opening input file ''', trim(input_obj % filename), ''''
        else if (input_obj % stream == STREAM_SFC) then
          write(stderrUnit,*) 'Error opening sfc file ''', trim(input_obj % filename), ''''
        end if
        write(stderrUnit,*) ' '
        call mpas_dmpar_abort(domain % dminfo)
      end if

!STEVE_edit: include some error catching to check the number of available grids
      if (mgslevel.eq.2) then! only catches the error if more than 1-layer - info storred in the multigrid *.nc file
        call MPAS_io_inq_dim(inputHandle, 'nGrids', nGrids, ierr)
        if (nGrids < config_num_mgslevels) then
          write(stderrUnit,*) 'config_num_mgslevels = ',config_num_mgslevels
          write(stderrUnit,*) 'should not exceed the total number of grids available = ',nGrids
          write(stderrUnit,*) ' '
          call mpas_dmpar_abort(domain % dminfo)
        endif
	print*,'nGrids', nGrids, mgslevel
      end if


      !
      ! Read global number of cells/edges/vertices
      !
! Steve_edit: I have edited the parser to include the additional mgs meshes in the following include
! using *_mgsLevel defined above as mgslevelChar
#include "read_dims.inc"
   
      !
      ! Determine the range of cells/edges/vertices that a processor will initially read
      !   from the input file
      !
      call mpas_dmpar_get_index_range(domain % dminfo, 1, nCells, readCellStart, readCellEnd)   
      nReadCells = readCellEnd - readCellStart + 1
   
      call mpas_dmpar_get_index_range(domain % dminfo, 1, nEdges, readEdgeStart, readEdgeEnd)   
      nReadEdges = readEdgeEnd - readEdgeStart + 1
   
      call mpas_dmpar_get_index_range(domain % dminfo, 1, nVertices, readVertexStart, readVertexEnd)   
      nReadVertices = readVertexEnd - readVertexStart + 1

      allocate(readingBlock)
      readingBlock % domain => domain
      readingBlock % blockID = domain % dminfo % my_proc_id
      readingBlock % localBlockID = 0
!STEVE_edit: Added the mgslevel and mgslevelChar to readingBlock
      readingBlock % mgsLevel = mgslevel
      readingBlock % mgslevelChar = mgslevelChar
      !
      ! Allocate and read fields that we will need in order to ultimately work out
      !   which cells/edges/vertices are owned by each block, and which are ghost
      !

      call mpas_io_setup_cell_block_fields(inputHandle, nreadCells, readCellStart, readingBlock, maxEdges, indexTocellIDField, xCellField, &
                                           yCellField, zCellField, nEdgesOnCellField, cellsOnCellField, edgesOnCellField, verticesOnCellField)

      call mpas_io_setup_edge_block_fields(inputHandle, nReadEdges, readEdgeStart, readingBlock, indexToEdgeIDField, xEdgeField, yEdgeField, zEdgeField, cellsOnEdgeField)

      call mpas_io_setup_vertex_block_fields(inputHandle, nReadVertices, readVertexStart, readingBlock, vertexDegree, indexToVertexIDField, &
                                             xVertexField, yVertexField, zVertexField, cellsOnVertexField)
      !
      ! Set up a graph derived data type describing the connectivity for the cells 
      !   that were read by this process
      ! A partial description is passed to the block decomp module by each process,
      !   and the block decomp module returns with a list of global cell indices
      !   that belong to the block on this process
      !
      partial_global_graph_info % nVertices = nReadCells
      partial_global_graph_info % nVerticesTotal = nCells
      partial_global_graph_info % maxDegree = maxEdges
      partial_global_graph_info % ghostStart = nVertices+1
      allocate(partial_global_graph_info % vertexID(nReadCells))
      allocate(partial_global_graph_info % nAdjacent(nReadCells))
      allocate(partial_global_graph_info % adjacencyList(maxEdges, nReadCells))
   
      partial_global_graph_info % vertexID(:) = indexToCellIDField % array(:)
      partial_global_graph_info % nAdjacent(:) = nEdgesOnCellField % array(:)
      partial_global_graph_info % adjacencyList(:,:) = cellsOnCellField % array(:,:)
   
      ! TODO: Ensure (by renaming or exchanging) that initial cell range on each proc is contiguous
      !       This situation may occur when reading a restart file with cells/edges/vertices written
      !       in a scrambled order
   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! NEED TO INCLUDE !!!!!!!!!!!!!!!!!!!!!  
!STEVE_edit: I need to include a routine here that generates the graph info
!somehow
!nCells is available
!Need coarseToFineCell for the current mesh


      ! Determine which cells are owned by this process
      call mpas_block_decomp_cells_for_proc(domain % dminfo, partial_global_graph_info, local_cell_list, block_id, block_start, block_count)

      deallocate(partial_global_graph_info % vertexID)
      deallocate(partial_global_graph_info % nAdjacent)
      deallocate(partial_global_graph_info % adjacencyList)

      call mpas_block_creator_setup_blocks_and_0halo_cells(domain, indexToCellID_Block, local_cell_list, block_id, block_start, block_count)
      call mpas_block_creator_build_0halo_cell_fields(indexToCellIDField, nEdgesOnCellField, cellsOnCellField, verticesOnCellField, edgesOnCellField, indexToCellID_Block, nEdgesOnCell_Block, cellsOnCell_Block, verticesOnCell_Block, edgesOnCell_Block)

      call mpas_block_creator_build_0_and_1halo_edge_fields(indexToEdgeIDField, cellsOnEdgeField, indexToCellID_Block, nEdgesOnCell_Block, edgesOnCell_Block, indexToEdgeID_Block, cellsOnEdge_Block, nEdgesSolveField)
      call mpas_block_creator_build_0_and_1halo_edge_fields(indexToVertexIDField, cellsOnVertexField, indexToCellID_Block, nEdgesOnCell_Block, verticesOnCell_Block, indexToVertexID_Block, cellsOnVertex_Block, nVerticesSolveField)

      call mpas_block_creator_build_cell_halos(indexToCellID_Block, nEdgesOnCell_Block, cellsOnCell_Block, verticesOnCell_Block, edgesOnCell_Block, nCellsSolveField)

      call mpas_block_creator_build_edge_halos(indexToCellID_Block, nEdgesOnCell_Block, nCellsSolveField, edgesOnCell_Block, indexToEdgeID_Block, cellsOnEdge_Block, nEdgesSolveField)
      call mpas_block_creator_build_edge_halos(indexToCellID_Block, nEdgesOnCell_Block, nCellsSolveField, verticesOnCell_Block, indexToVertexID_Block, cellsOnVertex_Block, nVerticesSolveField)

     ! Allocate blocks, and copy indexTo arrays into blocks

     call mpas_block_creator_finalize_block_init(domain % blocklist, &
#include "dim_dummy_args.inc"
                             , nCellsSolveField, nEdgesSolveField, nVerticesSolveField, indexToCellID_Block, indexToEdgeID_Block, indexToVertexID_Block, mgslevel)

     call mpas_io_input_init(input_obj, domain % blocklist, domain % dminfo)

     if (mgslevel==1) then
!STEVE_edit: just working on mgslevel = 1
!a local blocklist so that the input pointer is left unchanged

        call MPAS_readStreamAtt(input_obj % io_stream, 'sphere_radius', r_sphere_radius, ierr)
        if (ierr /= MPAS_STREAM_NOERR) then

          write(stderrUnit,*) 'Warning: Attribute sphere_radius not found in '//trim(input_obj % filename)
          write(stderrUnit,*) '   Setting sphere_radius to 1.0'
          domain % blocklist % mesh % sphere_radius = 1.0
        else
          domain % blocklist % mesh % sphere_radius = r_sphere_radius
        end if

        call MPAS_readStreamAtt(input_obj % io_stream, 'on_a_sphere', c_on_a_sphere, ierr)
        if (ierr /= MPAS_STREAM_NOERR) then
          write(stderrUnit,*) 'Warning: Attribute on_a_sphere not found in '//trim(input_obj % filename)
          write(stderrUnit,*) '   Setting on_a_sphere to ''YES'''
          domain % blocklist % mesh % on_a_sphere = .true.
        else
          if (index(c_on_a_sphere, 'YES') /= 0) then
            domain % blocklist % mesh % on_a_sphere = .true.
          else
            domain % blocklist % mesh % on_a_sphere = .false.
          end if
        end if
!STEVE_edit: Note: these directives have been added to ver 2.0
#ifndef MPAS_CESM
        call MPAS_readStreamAtt(input_obj % io_stream, 'history', domain % history, ierr)
        if (ierr /= MPAS_STREAM_NOERR) then
          write(stderrUnit,*) 'Warning: Attribute History not found in '//trim(input_obj % filename)
          write(stderrUnit,*) '   Setting History to '''''
          domain % history = ""
        else
          ! Remove C String NULL characters, replace C String newlines with semicolons
          do i = 1, len(domain % history)
            if(iachar(domain % history(i:i)) == 0) then
              domain % history(i:i) = " "
            else if(iachar(domain % history(i:i)) == 10) then
              domain % history(i:i) = ";"
            end if
          end do
        end if
#else
        domain % history = "cesm_run"
#endif
      end if


!STEVE_edit: No longer pointing at next because I want to populate all blocks so I start at the top level
      !block_ptr => domain % blocklist % next
      block_ptr => domain % blocklist
!STEVE_edit: Loop down to the coarse level
      do while (associated(block_ptr % coarser))
        block_ptr => block_ptr % coarser
      end do

      do while (associated(block_ptr))
        block_ptr % mesh % sphere_radius = domain % blocklist % mesh % sphere_radius
        block_ptr % mesh % on_a_sphere = domain % blocklist % mesh % on_a_sphere

        ! Link the sendList and recvList pointers in each field type to the appropriate lists 
        !   in parinfo, e.g., cellsToSend and cellsToRecv; in future, it can also be extended to 
        !   link blocks of fields to eachother
!STEVE_edit: because I start at top instead of next, I nedd a condition to stop
!mpas_create_field_links being run
        if (associated(block_ptr % prev)) then
          call mpas_create_field_links(block_ptr)
        end if

        block_ptr => block_ptr % next
      end do

      if (.not. config_do_restart) then
        input_obj % time = 1
      else
        !
        ! If doing a restart, we need to decide which time slice to read from the 
        !   restart file
        !
        input_obj % time = MPAS_seekStream(input_obj % io_stream, restartTimeStamp, MPAS_STREAM_EXACT_TIME, timeStamp, ierr)
        if (ierr == MPAS_IO_ERR) then
          write(stderrUnit,*) 'Error: restart file '//trim(filename)//' did not contain time '//trim(restartTimeStamp)
          call mpas_dmpar_abort(domain % dminfo)
        end if

!       input_obj % time = MPAS_seekStream(input_obj % io_stream, config_start_time, MPAS_STREAM_EXACT_TIME, timeStamp, ierr)
!       if (ierr == MPAS_IO_ERR) then
!         write(stderrUnit,*) 'Error: restart file '//trim(filename)//' did not contain time '//trim(config_start_time)
!         call mpas_dmpar_abort(domain % dminfo)
!       end if
!write(stderrUnit,*) 'MGD DEBUGGING time = ', input_obj % time
        write(stderrUnit,*) 'Restarting model from time ', trim(timeStamp)
      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      ! Do the actual work of reading all fields in from the input or restart file
      ! For each field:
      !   1) Each process reads a contiguous range of cell/edge/vertex indices, which
      !      may not correspond with the cells/edges/vertices that are owned by the
      !      process
      !   2) All processes then send the global indices that were read to the 
      !      processes that own those indices based on 
      !      {send,recv}{Cell,Edge,Vertex}List
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      call mpas_read_and_distribute_fields(input_obj)

      call mpas_io_input_finalize(input_obj, domain % dminfo)

      call MPAS_io_close(inputHandle, ierr)

!STEVE_edit: adding code to poupulate dimension fieilds that are not available in the mg_*nc file or are in there but incorect. Ultimately, all of the information could be storred in init file - I need to do better than this!
!STEVE_edit: Loop down to the coarse level starting from the top
      block_ptr => domain % blocklist
      do while (associated(block_ptr % coarser))
         block_ptr => block_ptr % coarser
      end do

      !
      ! Exchange halos for all of the fields that were read from the input file
      !
      call mpas_exch_input_field_halos(domain, input_obj)

      call mpas_block_creator_reindex_block_fields(domain % blocklist)

      call mpas_dmpar_destroy_mulithalo_exchange_list(indexToCellIDField % sendList)
      call mpas_dmpar_destroy_mulithalo_exchange_list(indexToCellIDField % recvList)
      call mpas_dmpar_destroy_mulithalo_exchange_list(indexToCellIDField % copyList)

      call mpas_dmpar_destroy_mulithalo_exchange_list(indexToEdgeIDField % sendList)
      call mpas_dmpar_destroy_mulithalo_exchange_list(indexToEdgeIDField % recvList)
      call mpas_dmpar_destroy_mulithalo_exchange_list(indexToEdgeIDField % copyList)

      call mpas_dmpar_destroy_mulithalo_exchange_list(indexToVertexIDField % sendList)
      call mpas_dmpar_destroy_mulithalo_exchange_list(indexToVertexIDField % recvList)
      call mpas_dmpar_destroy_mulithalo_exchange_list(indexToVertexIDField % copyList)

      call mpas_deallocate_field(indexToCellIDField)
      call mpas_deallocate_field(indexToEdgeIDField)
      call mpas_deallocate_field(indexToVertexIDField)
      call mpas_deallocate_field(cellsOnCellField)

      call mpas_deallocate_field(edgesOnCellField)
      call mpas_deallocate_field(verticesOnCellField)
      call mpas_deallocate_field(cellsOnEdgeField)
      call mpas_deallocate_field(cellsOnVertexField)

      call mpas_deallocate_field(indexToCellID_Block)
      call mpas_deallocate_field(nEdgesOnCell_Block)
      call mpas_deallocate_field(cellsOnCell_Block)
      call mpas_deallocate_field(verticesOnCell_Block)
      call mpas_deallocate_field(edgesOnCell_Block)
      call mpas_deallocate_field(indexToVertexID_Block)
      call mpas_deallocate_field(cellsOnVertex_Block)
      call mpas_deallocate_field(indexToEdgeID_Block)
      call mpas_deallocate_field(cellsOnEdge_Block)

      call mpas_deallocate_field(nCellsSolveField)
      call mpas_deallocate_field(nVerticesSolveField)
      call mpas_deallocate_field(nEdgesSolveField)

      deallocate(local_cell_list)
      deallocate(block_id)
      deallocate(block_start)
      deallocate(block_count)
      deallocate(readingBlock)

   end subroutine mpas_input_state_for_domain_iter!}}}

   !CR:TODO: an identical subroutine is found in module_io_output - merge
   subroutine mpas_insert_string_suffix(stream, suffix, filename)!{{{

      implicit none

      character (len=*), intent(in) :: stream
      character (len=*), intent(in) :: suffix
      character (len=*), intent(out) :: filename
      integer :: length, i

      filename = trim(stream) // '.' // trim(suffix)

      length = len_trim(stream)
      do i=length-1,1,-1
         if(stream(i:i) == '.') then
            filename = trim(stream(:i)) // trim(suffix) // trim(stream(i:))
            exit
         end if
      end do

      do i=1,len_trim(filename)
         if (filename(i:i) == ':') filename(i:i) = '.'
      end do

   end subroutine mpas_insert_string_suffix!}}}

   subroutine mpas_read_and_distribute_fields(input_obj)!{{{
      
      implicit none

      type (io_input_object), intent(inout) :: input_obj

      integer :: ierr

      call MPAS_readStream(input_obj % io_stream, input_obj % time, ierr)


   end subroutine mpas_read_and_distribute_fields!}}}

!   subroutine mpas_io_input_init()
!
!   Edited to deal with multiple grids
!
   subroutine mpas_io_input_init(input_obj, blocklistIN, dminfo)!{{{
! change the code so that it does not run this subR for levels below finE 
      implicit none

      type (io_input_object), intent(inout) :: input_obj
!      type (block_type), intent(in) :: blocklist
!STEVE_edit: Ive changed the input block type useing (I could set intent inout?)
      type (block_type), pointer :: blocklistIN
      type (dm_info), intent(in) :: dminfo

!STEVE_edit: new local blocklist - to handle the multi-grid
      type (block_type),pointer  :: blocklist
 
      integer :: nferr

!STEVE_edit: I dont understand why this doesnt include looping of % next as there
!could be multiple meshes on each proc
!STEVE_edit: Im useing a blocklistIN so that the position in the blocklist is not
!moved from the head. This may not matter but its a simple case to swich back:
     blocklist => blocklistIN

!STEVE_edit: looping to the % coarser level
      do while (associated(blocklist % coarser))
        blocklist => blocklist % coarser
      end do
 
      call MPAS_createStream(input_obj % io_stream, trim(input_obj % filename), MPAS_IO_PNETCDF, MPAS_IO_READ, 1, nferr)
!STEVE_edit: this could be included in the above call?
!STEVE_edit: this passes the mgsLevel into the stream reading
      
      input_obj % io_stream % mgslevel = blocklist % mgsLevel

      if (nferr /= MPAS_STREAM_NOERR) then
         write(stderrUnit,*) ' '
         if (input_obj % stream == STREAM_RESTART) then
            write(stderrUnit,*) 'Error opening restart file ''', trim(input_obj % filename), ''''
         else if (input_obj % stream == STREAM_INPUT) then
            write(stderrUnit,*) 'Error opening input file ''', trim(input_obj % filename), ''''
         else if (input_obj % stream == STREAM_SFC) then
            write(stderrUnit,*) 'Error opening sfc file ''', trim(input_obj % filename), ''''
         end if
         write(stderrUnit,*) ' '
         call mpas_dmpar_abort(dminfo)
      end if
!STEVE_edit: use an updated include script
!this is reading all the varibles but I dont need to do this for the
!subgrids so I sould look to makeing an additional include script.
#include "add_input_fields.inc"

   end subroutine mpas_io_input_init!}}}

!  subroutine mpas_exch_input_field_halos()
!
!  Edited to deal with multiple grids
!
   subroutine mpas_exch_input_field_halos(domain, input_obj)!{{{

      implicit none

      type (domain_type), intent(inout) :: domain
      type (io_input_object), intent(inout) :: input_obj
!STEVE_edit: add a local blocklist
      type(block_type), pointer ::blocklist

!STEVE_edit: original code was using domian % blocklist but now I am using local
!blocklist and this has been cahnged in the *inc files below. I works on the %
!coarser block
      blocklist => domain % blocklist
      do while (associated(blocklist % coarser))
        blocklist => blocklist % coarser
      end do

#include "exchange_input_field_halos.inc"

#include "non_decomp_copy_input_fields.inc"

   end subroutine mpas_exch_input_field_halos!}}}

   subroutine mpas_io_input_finalize(input_obj, dminfo)!{{{
 
      implicit none
 
      type (io_input_object), intent(inout) :: input_obj
      type (dm_info), intent(in) :: dminfo

      integer :: nferr
 
      call MPAS_closeStream(input_obj % io_stream, nferr)
 
   end subroutine mpas_io_input_finalize!}}}

!  subroutine mpas_io_setup_cell_block_fields
!
!       STEVE_edit: edited thir routine to handle multi-grid
!
   subroutine mpas_io_setup_cell_block_fields(inputHandle, nReadCells, readCellStart, readingBlock, maxEdges, indexToCellID, xCell, yCell, zCell, nEdgesOnCell, cellsOnCell, edgesOnCell, verticesOnCell)!{{{
     type (MPAS_IO_Handle_type) :: inputHandle
     integer, intent(in) :: nReadCells
     integer, intent(in) :: readCellStart
     integer, intent(in) :: maxEdges
     type (block_type), pointer :: readingBlock
     type (field1dInteger), pointer :: indexToCellID
     type (field1dReal), pointer :: xCell
     type (field1dReal), pointer :: yCell
     type (field1dReal), pointer :: zCell
     type (field1dInteger), pointer :: nEdgesOnCell
     type (field2dInteger), pointer :: cellsOnCell
     type (field2dInteger), pointer :: edgesOnCell
     type (field2dInteger), pointer :: verticesOnCell

     integer :: i, nHalos
     integer, dimension(:), pointer :: readIndices
!STEVE_edit: dumChar is used to make a string to read the correct varible from grid 
     character(len=StrKIND) :: dumChar

     nHalos = config_num_halos
  
     !
     ! Allocate and read fields that we will need in order to ultimately work out
     !   which cells/edges/vertices are owned by each block, and which are ghost
     !

     ! Global cell indices
!STEVE_edit: I am using a varible dumChar to pass the varible to be read as there
!are multiple grids. For this grid I have kept the original code (comented out)
!for refferance
     dumChar='indexToCellID'//trim(readingBlock % mgslevelChar)
     allocate(indexToCellID)
     allocate(indexToCellID % ioinfo)
!     indexToCellID % ioinfo % fieldName = 'indexToCellID'
     indexToCellID % ioinfo % fieldName = trim(dumChar)
     indexToCellID % ioinfo % start(1) = readCellStart
     indexToCellID % ioinfo % count(1) = nReadCells
     allocate(indexToCellID % array(nReadCells))
     allocate(readIndices(nReadCells))
     do i=1,nReadCells
        readIndices(i) = i + readCellStart - 1
     end do
!     call MPAS_io_inq_var(inputHandle, 'indexToCellID', ierr=ierr)
!     call MPAS_io_set_var_indices(inputHandle, 'indexToCellID', readIndices, ierr=ierr)
!     call mpas_io_get_var(inputHandle, 'indexToCellID', indexToCellID % array, ierr)
     call MPAS_io_inq_var(inputHandle, trim(dumChar), ierr=ierr)
     call MPAS_io_set_var_indices(inputHandle, trim(dumChar), readIndices,ierr=ierr)
     call mpas_io_get_var(inputHandle, trim(dumChar), indexToCellID % array,ierr)

     indexToCellID % dimSizes(1) = nReadCells
     indexToCellID % block => readingBlock
     call mpas_dmpar_init_multihalo_exchange_list(indexToCellID % sendList, nHalos)
     call mpas_dmpar_init_multihalo_exchange_list(indexToCellID % recvList, nHalos)
     call mpas_dmpar_init_multihalo_exchange_list(indexToCellID % copyList, nHalos)
     nullify(indexToCellID % next)
   

     ! Number of cell/edges/vertices adjacent to each cell
     dumChar='nEdgesOnCell'//trim(readingBlock % mgslevelChar)
     allocate(nEdgesOnCell)
     allocate(nEdgesOnCell % ioinfo)
!     nEdgesOnCell % ioinfo % fieldName = 'nEdgesOnCell'
     nEdgesOnCell % ioinfo % fieldName = trim(dumChar)
     nEdgesOnCell % ioinfo % start(1) = readCellStart
     nEdgesOnCell % ioinfo % count(1) = nReadCells
     allocate(nEdgesOnCell % array(nReadCells))
!     call MPAS_io_inq_var(inputHandle, 'nEdgesOnCell', ierr=ierr)
!     call MPAS_io_set_var_indices(inputHandle, 'nEdgesOnCell', readIndices, ierr=ierr)
!     call mpas_io_get_var(inputHandle, 'nEdgesOnCell', nEdgesOnCell % array, ierr)
     call MPAS_io_inq_var(inputHandle, trim(dumChar), ierr=ierr)
     call MPAS_io_set_var_indices(inputHandle, trim(dumChar), readIndices,ierr=ierr)
     call mpas_io_get_var(inputHandle, trim(dumChar), nEdgesOnCell % array,ierr)
     nEdgesOnCell % dimSizes(1) = nReadCells
     nEdgesOnCell % block => readingBlock
     nEdgesOnCell % sendList => indexToCellID % sendList
     nEdgesOnCell % recvList => indexToCellID % recvList
     nEdgesOnCell % copyList => indexToCellID % copyList
     nullify(nEdgesOnCell % next)
   
     ! Global indices of cells adjacent to each cell
     dumChar='cellsOnCell'//trim(readingBlock % mgslevelChar)
     allocate(cellsOnCell)
     allocate(cellsOnCell % ioinfo)
!     cellsOnCell % ioinfo % fieldName = 'cellsOnCell'
     cellsOnCell % ioinfo % fieldName = trim(dumChar)
     cellsOnCell % ioinfo % start(1) = 1
     cellsOnCell % ioinfo % start(2) = readCellStart
     cellsOnCell % ioinfo % count(1) = maxEdges
     cellsOnCell % ioinfo % count(2) = nReadCells
     allocate(cellsOnCell % array(maxEdges,nReadCells))
!     call MPAS_io_inq_var(inputHandle, 'cellsOnCell', ierr=ierr)
!     call MPAS_io_set_var_indices(inputHandle, 'cellsOnCell', readIndices, ierr=ierr)
!     call mpas_io_get_var(inputHandle, 'cellsOnCell', cellsOnCell % array, ierr)
     call MPAS_io_inq_var(inputHandle, trim(dumChar), ierr=ierr)
     call MPAS_io_set_var_indices(inputHandle, trim(dumChar), readIndices,ierr=ierr)
     call mpas_io_get_var(inputHandle, trim(dumChar), cellsOnCell % array,ierr)
     cellsOnCell % dimSizes(1) = maxEdges
     cellsOnCell % dimSizes(2) = nReadCells
     cellsOnCell % block => readingBlock
     cellsOnCell % sendList => indexToCellID % sendList
     cellsOnCell % recvList => indexToCellID % recvList
     cellsOnCell % copyList => indexToCellID % copyList
     nullify(cellsOnCell % next)
   
     ! Global indices of edges adjacent to each cell
     dumChar='edgesOnCell'//trim(readingBlock % mgslevelChar)
     allocate(edgesOnCell)
     allocate(edgesOnCell % ioinfo)
!     edgesOnCell % ioinfo % fieldName = 'edgesOnCell'
     edgesOnCell % ioinfo % fieldName = trim(dumChar)
     edgesOnCell % ioinfo % start(1) = 1
     edgesOnCell % ioinfo % start(2) = readCellStart
     edgesOnCell % ioinfo % count(1) = maxEdges
     edgesOnCell % ioinfo % count(2) = nReadCells
     allocate(edgesOnCell % array(maxEdges,nReadCells))
!     call MPAS_io_inq_var(inputHandle, 'edgesOnCell', ierr=ierr)
!     call MPAS_io_set_var_indices(inputHandle, 'edgesOnCell', readIndices, ierr=ierr)
!     call mpas_io_get_var(inputHandle, 'edgesOnCell', edgesOnCell % array, ierr)
     call MPAS_io_inq_var(inputHandle, trim(dumChar), ierr=ierr)
     call MPAS_io_set_var_indices(inputHandle, trim(dumChar), readIndices,ierr=ierr)
     call mpas_io_get_var(inputHandle, trim(dumChar), edgesOnCell % array,ierr)
     edgesOnCell % dimSizes(1) = maxEdges
     edgesOnCell % dimSizes(2) = nReadCells
     edgesOnCell % block => readingBlock
     edgesOnCell % sendList => indexToCellID % sendList
     edgesOnCell % recvList => indexToCellID % recvList
     edgesOnCell % copyList => indexToCellID % copyList
     nullify(edgesOnCell % next)
   
     ! Global indices of vertices adjacent to each cell
     dumChar='verticesOnCell'//trim(readingBlock % mgslevelChar)
     allocate(verticesOnCell)
     allocate(verticesOnCell % ioinfo)
!     verticesOnCell % ioinfo % fieldName = 'verticesOnCell'
     verticesOnCell % ioinfo % fieldName = trim(dumChar)
     verticesOnCell % ioinfo % start(1) = 1
     verticesOnCell % ioinfo % start(2) = readCellStart
     verticesOnCell % ioinfo % count(1) = maxEdges
     verticesOnCell % ioinfo % count(2) = nReadCells
     allocate(verticesOnCell % array(maxEdges,nReadCells))
!     call MPAS_io_inq_var(inputHandle, 'verticesOnCell', ierr=ierr)
!     call MPAS_io_set_var_indices(inputHandle, 'verticesOnCell', readIndices, ierr=ierr)
!     call mpas_io_get_var(inputHandle, 'verticesOnCell', verticesOnCell % array, ierr)
     call MPAS_io_inq_var(inputHandle, trim(dumChar), ierr=ierr)
     call MPAS_io_set_var_indices(inputHandle, trim(dumChar), readIndices,ierr=ierr)
     call mpas_io_get_var(inputHandle, trim(dumChar), verticesOnCell % array,ierr)
     verticesOnCell % dimSizes(1) = maxEdges
     verticesOnCell % dimSizes(2) = nReadCells
     verticesOnCell % block => readingBlock
     verticesOnCell % sendList => indexToCellID % sendList
     verticesOnCell % recvList => indexToCellID % recvList
     verticesOnCell % copyList => indexToCellID % copyList
     nullify(verticesOnCell % next)

     deallocate(readIndices)
   
   end subroutine mpas_io_setup_cell_block_fields!}}}

!  subroutine mpas_io_setup_edge_block_fields()
!
!       STEVE_edit: edited thir routine to handle multi-grid
!
   subroutine mpas_io_setup_edge_block_fields(inputHandle, nReadEdges, readEdgeStart, readingBlock, indexToEdgeID, xEdge, yEdge, zEdge, cellsOnEdge)!{{{
     type (MPAS_IO_Handle_type) :: inputHandle
     integer, intent(in) :: nReadEdges
     integer, intent(in) :: readEdgeStart
     type (block_type), pointer :: readingBlock
     type (field1dInteger), pointer :: indexToEdgeID
     type (field1dReal), pointer :: xEdge
     type (field1dReal), pointer :: yEdge
     type (field1dReal), pointer :: zEdge
     type (field2dInteger), pointer :: cellsOnEdge

     integer :: i, nHalos
     integer, dimension(:), pointer :: readIndices
!STEVE_edit: dumChar is used to make a string to read the correct varible from
!grid 
     character(len=StrKIND) :: dumChar

     nHalos = config_num_halos
  
     !
     ! Allocate and read fields that we will need in order to ultimately work out
     !   which cells/edges/vertices are owned by each block, and which are ghost
     !

     allocate(readIndices(nReadEdges))

     ! Global edge indices
     dumChar='indexToEdgeID'//trim(readingBlock % mgslevelChar)
     allocate(indexToEdgeID)
     allocate(indexToEdgeID % ioinfo)
!     indexToEdgeID % ioinfo % fieldName = 'indexToEdgeID'
     indexToEdgeID % ioinfo % fieldName = trim(dumChar)
     indexToEdgeID % ioinfo % start(1) = readEdgeStart
     indexToEdgeID % ioinfo % count(1) = nReadEdges
     allocate(indexToEdgeID % array(nReadEdges))
     allocate(indexToEdgeID % array(nReadEdges))
     do i=1,nReadEdges
        readIndices(i) = i + readEdgeStart - 1
     end do
!     call MPAS_io_inq_var(inputHandle, 'indexToEdgeID', ierr=ierr)
!     call MPAS_io_set_var_indices(inputHandle, 'indexToEdgeID', readIndices, ierr=ierr)
!     call mpas_io_get_var(inputHandle, 'indexToEdgeID', indexToEdgeID % array, ierr)
     call MPAS_io_inq_var(inputHandle, trim(dumChar), ierr=ierr)
     call MPAS_io_set_var_indices(inputHandle, trim(dumChar), readIndices,ierr=ierr)
     call mpas_io_get_var(inputHandle, trim(dumChar), indexToEdgeID % array,ierr)
     indexToEdgeID % dimSizes(1) = nREadEdges
     indexToEdgeID % block => readingBlock
     call mpas_dmpar_init_multihalo_exchange_list(indexToEdgeID % sendList, nHalos+1)
     call mpas_dmpar_init_multihalo_exchange_list(indexToEdgeID % recvList, nHalos+1)
     call mpas_dmpar_init_multihalo_exchange_list(indexToEdgeID % copyList, nHalos+1)
     nullify(indexToEdgeID % next)
   
   
     ! Global indices of cells adjacent to each edge
     !    used for determining which edges are owned by a block, where 
     !    iEdge is owned iff cellsOnEdge(1,iEdge) is an owned cell
     dumChar='cellsOnEdge'//trim(readingBlock % mgslevelChar)
     allocate(cellsOnEdge)
     allocate(cellsOnEdge % ioinfo)
!     cellsOnEdge % ioinfo % fieldName = 'cellsOnEdge'
     cellsOnEdge % ioinfo % fieldName = trim(dumChar)
     cellsOnEdge % ioinfo % start(1) = 1
     cellsOnEdge % ioinfo % start(2) = readEdgeStart
     cellsOnEdge % ioinfo % count(1) = 2
     cellsOnEdge % ioinfo % count(2) = nReadEdges
     allocate(cellsOnEdge % array(2,nReadEdges))
!     call MPAS_io_inq_var(inputHandle, 'cellsOnEdge', ierr=ierr)
!     call MPAS_io_set_var_indices(inputHandle, 'cellsOnEdge', readIndices, ierr=ierr)
!     call mpas_io_get_var(inputHandle, 'cellsOnEdge', cellsOnEdge % array, ierr)
     call MPAS_io_inq_var(inputHandle, trim(dumChar), ierr=ierr)
     call MPAS_io_set_var_indices(inputHandle, trim(dumChar), readIndices,ierr=ierr)
     call mpas_io_get_var(inputHandle, trim(dumChar), cellsOnEdge % array,ierr)
     cellsOnEdge % dimSizes(1) = 2
     cellsOnEdge % dimSizes(2) = nReadEdges
     cellsOnEdge % block => readingBlock
     cellsOnEdge % sendList => indexToEdgeID % sendList
     cellsOnEdge % recvList => indexToEdgeID % recvList
     cellsOnEdge % copyList => indexToEdgeID % copyList
     nullify(cellsOnEdge % next)

     deallocate(readIndices)
   
   end subroutine mpas_io_setup_edge_block_fields!}}}

!  subroutine mpas_io_setup_vertex_block_fields()
!
!       STEVE_edit: edited thir routine to handle multi-grid
!
   subroutine mpas_io_setup_vertex_block_fields(inputHandle, nReadVertices, readVertexStart, readingBlock, vertexDegree, indexToVertexID, xVertex, yVertex, zVertex, cellsOnVertex)!{{{
     type (MPAS_IO_Handle_type) :: inputHandle
     integer, intent(in) :: nReadVertices
     integer, intent(in) :: readVertexStart
     integer, intent(in) :: vertexDegree
     type (block_type), pointer :: readingBlock
     type (field1dInteger), pointer :: indexToVertexID
     type (field1dReal), pointer :: xVertex
     type (field1dReal), pointer :: yVertex
     type (field1dReal), pointer :: zVertex
     type (field2dInteger), pointer :: cellsOnVertex

     integer :: i, nHalos
     integer, dimension(:), pointer :: readIndices
!STEVE_edit: dumChar is used to make a string to read the correct varible from
!grid 
     character(len=StrKIND) :: dumChar

     nHalos = config_num_halos
  
     ! Global vertex indices
     dumChar='indexToVertexID'//trim(readingBlock % mgslevelChar)
     allocate(indexToVertexID)
     allocate(indexToVertexID % ioinfo)
!     indexToVertexID % ioinfo % fieldName = 'indexToVertexID'
     indexToVertexID % ioinfo % fieldName = trim(dumChar)
     indexToVertexID % ioinfo % start(1) = readVertexStart
     indexToVertexID % ioinfo % count(1) = nReadVertices
     allocate(indexToVertexID % array(nReadVertices))
     allocate(readIndices(nReadVertices))
     do i=1,nReadVertices
        readIndices(i) = i + readVertexStart - 1
     end do
!     call MPAS_io_inq_var(inputHandle, 'indexToVertexID', ierr=ierr)
!     call MPAS_io_set_var_indices(inputHandle, 'indexToVertexID', readIndices, ierr=ierr)
!     call mpas_io_get_var(inputHandle, 'indexToVertexID', indexToVertexID % array, ierr)
     call MPAS_io_inq_var(inputHandle, trim(dumChar), ierr=ierr)
     call MPAS_io_set_var_indices(inputHandle, trim(dumChar), readIndices,ierr=ierr)
     call mpas_io_get_var(inputHandle, trim(dumChar), indexToVertexID % array,ierr)
     indexToVertexID % dimSizes(1) = nReadVertices
     indexToVertexID % block => readingBlock
     call mpas_dmpar_init_multihalo_exchange_list(indexToVertexID % sendList, nHalos+1)
     call mpas_dmpar_init_multihalo_exchange_list(indexToVertexID % recvList, nHalos+1)
     call mpas_dmpar_init_multihalo_exchange_list(indexToVertexID % copyList, nHalos+1)
     nullify(indexToVertexID % next)   

   
     ! Global indices of cells adjacent to each vertex
     !    used for determining which vertices are owned by a block, where 
     !    iVtx is owned iff cellsOnVertex(1,iVtx) is an owned cell
     dumChar='cellsOnVertex'//trim(readingBlock % mgslevelChar)
     allocate(cellsOnVertex)
     allocate(cellsOnVertex % ioinfo)
!     cellsOnVertex % ioinfo % fieldName = 'cellsOnVertex'
     cellsOnVertex % ioinfo % fieldName = trim(dumChar)
     cellsOnVertex % ioinfo % start(1) = 1
     cellsOnVertex % ioinfo % start(2) = readVertexStart
     cellsOnVertex % ioinfo % count(1) = vertexDegree
     cellsOnVertex % ioinfo % count(2) = nReadVertices
     allocate(cellsOnVertex % array(vertexDegree,nReadVertices))
!     call MPAS_io_inq_var(inputHandle, 'cellsOnVertex', ierr=ierr)
!     call MPAS_io_set_var_indices(inputHandle, 'cellsOnVertex', readIndices, ierr=ierr)
!     call mpas_io_get_var(inputHandle, 'cellsOnVertex', cellsOnVertex % array, ierr)
     call MPAS_io_inq_var(inputHandle, trim(dumChar), ierr=ierr)
     call MPAS_io_set_var_indices(inputHandle, trim(dumChar), readIndices,ierr=ierr)
     call mpas_io_get_var(inputHandle, trim(dumChar), cellsOnVertex % array,ierr)
     cellsOnVertex % dimSizes(1) = vertexDegree
     cellsOnVertex % dimSizes(2) = nReadVertices
     cellsOnVertex % block => readingBlock
     cellsOnVertex % sendList => indexToVertexID % sendList
     cellsOnVertex % recvList => indexToVertexID % recvList
     cellsOnVertex % copyList => indexToVertexID % copyList
     nullify(cellsOnVertex % next)

     deallocate(readIndices)

   end subroutine mpas_io_setup_vertex_block_fields!}}}
 
end module mpas_io_input
