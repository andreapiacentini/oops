module qg_geom_mod
   use aq_constants_mod
   use fckit_module
   use atlas_module

   implicit none

   private
   public :: qg_geom
   
   type qg_geom
      type(atlas_StructuredGrid)                  :: grid
      type(fckit_mpi_comm)                        :: fmpi
      character(len=:), allocatable               :: domain
      character(len=:), allocatable               :: model
      character(len=:), allocatable               :: orientation
      type(atlas_functionspace_StructuredColumns) :: fs_2d
      type(atlas_functionspace_StructuredColumns) :: fs_3d
      integer(atlas_kind_idx)                     :: levels
      integer(atlas_kind_idx)                     :: mod_levels
      integer(atlas_kind_idx)                     :: halo
      integer(atlas_kind_idx)                     :: bbox_imin, bbox_imax, bbox_isz
      integer(atlas_kind_idx)                     :: bbox_jmin, bbox_jmax, bbox_jsz
      real(aq_real), allocatable                  :: halo_mask(:)
      real(aq_real)                               :: deltax, deltay
      real(aq_real)                               :: x0, y0
      integer(atlas_kind_idx)                     :: nx, ny, nz
   contains
      procedure, public :: create => aq_geom_create
      procedure, public :: clone  => aq_geom_clone
      procedure, public :: info   => aq_geom_info
      procedure, public :: delete => aq_geom_delete
      procedure, public :: final  => aq_geom_delete
   end type qg_geom

contains
   
   subroutine aq_geom_create(self, config, fckit_mpi)
      class(qg_geom), intent(inout)    :: self
      type(fckit_Configuration)        :: config
      type(fckit_mpi_comm), intent(in) :: fckit_mpi
      !
      character(len=:), allocatable :: str
      !
      type(atlas_Partitioner) :: partitioner
      type(atlas_GridDistribution) :: griddistribution
      real(oops_real) :: dom_sw_corner(2)
      integer, allocatable :: ila_work(:)
      integer(atlas_kind_idx) :: ib_i, ib_j
      !
      self%fmpi = fckit_mpi
      call config%get_or_die("nx",self%nx)
      call config%get_or_die("ny",self%ny)
      call config%get_or_die("dx",self%deltax)
      call config%get_or_die("dy",self%deltay)
      call config%get_or_die("x0",self%x0)
      call config%get_or_die("y0",self%y0)
      call config%get_or_die("domain",self%domain)
      dom_sw_corner = [self%x0, self%y0]
      self%grid = atlas_RegionalGrid(self%nx,self%ny,self%deltax,self%deltay,dom_sw_corner)
      partitioner = atlas_Partitioner(type="checkerboard")
      griddistribution = partitioner%partition(self%grid)
      !
      self%halo = 0
      if (config%has("halo")) call config%get_or_die("halo",self%halo)
      !
      self%model = "MOCAGE"
      if (config%has("model")) call config%get_or_die("model",self%model)
      !
      self%orientation = "up"
      if (config%has("orientation")) then
         call config%get_or_die("orientation",str)
         select case(str)
            case("up")
               self%orientation = str
            case("UP")
               self%orientation = "up"
            case("down")
               self%orientation = str
            case("DOWN")
               self%orientation = "down"
            case default
               call abor1_ftn("Geom orientation can only be 'up' or 'down'")
            end select
      end if
      !
      call config%get_or_die("levels",self%levels)
      self%mod_levels = self%levels
      !
      self%nz = self%levels
      !
      self%fs_2d = atlas_functionspace_StructuredColumns(self%grid, &
         &                                               griddistribution, &
         &                                               halo=self%halo)
      !
      self%fs_3d = atlas_functionspace_StructuredColumns(self%grid, &
         &                                               griddistribution, &
         &                                               halo=self%halo, &
         &                                               levels=self%levels)
      !
      call griddistribution%final()
      call partitioner%final()
      !
      ! Bounding box global indexes. Notice that C++ atlas has already a similar function
      self%bbox_jmin = self%fs_3d%j_begin()
      self%bbox_jmax = self%fs_3d%j_end()
      allocate(ila_work(self%bbox_jmin:self%bbox_jmax))
      do ib_j = self%bbox_jmin, self%bbox_jmax
         ila_work(ib_j) = self%fs_3d%i_begin(ib_j)
      end do
      self%bbox_imin = minval(ila_work)
      do ib_j = self%bbox_jmin, self%bbox_jmax
         ila_work(ib_j) = self%fs_3d%i_end(ib_j)
      end do
      self%bbox_imax = maxval(ila_work)
      deallocate(ila_work)
      !
      if ( self%halo>0 ) then
         allocate( self%halo_mask(self%fs_3d%size()) )
         self%halo_mask(:) = 0.0
         do ib_j = self%fs_3d%j_begin(), self%fs_3d%j_end()
            do ib_i = self%fs_3d%i_begin(ib_j), self%fs_3d%i_end(ib_j)
               self%halo_mask(self%fs_3d%index(ib_i,ib_j)) = 1.0
            end do
         end do
      end if
      !
   end subroutine aq_geom_create

   subroutine aq_geom_clone(self, other)
      class(qg_geom), intent(inout) :: self
      class(qg_geom), intent(in)    :: other
      !
      self%grid = other%grid
      self%nx = other%nx
      self%ny = other%ny
      self%nz = other%nz
      self%deltax = other%deltax
      self%deltay = other%deltay
      self%x0 = other%x0
      self%y0 = other%y0
      self%fmpi = other%fmpi
      self%domain = other%domain
      self%model = other%model
      self%orientation = other%orientation
      self%levels = other%levels
      self%mod_levels = other%mod_levels
      self%halo = other%halo
      self%bbox_imin = other%bbox_imin 
      self%bbox_imax = other%bbox_imax 
      self%bbox_isz  = other%bbox_isz  
      self%bbox_jmin = other%bbox_jmin 
      self%bbox_jmax = other%bbox_jmax 
      self%bbox_jsz  = other%bbox_jsz
      self%fs_2d = other%fs_2d
      self%fs_3d = other%fs_3d
      if ( self%halo>0 ) then
         allocate( self%halo_mask(size(other%halo_mask)) )
         self%halo_mask(:) = other%halo_mask(:)
      end if
      !
   end subroutine aq_geom_clone
   
   subroutine aq_geom_delete(self)
      class(qg_geom), intent(inout)  :: self
      call self%fs_2d%final()
      call self%fs_3d%final()
      call self%grid%final()
      if (allocated(self%halo_mask)) deallocate(self%halo_mask)
   end subroutine aq_geom_delete

   subroutine aq_geom_info(self, nx, ny, nz, deltax, deltay, x0, y0, halo, domain, orientation, model)
      class(qg_geom), intent(in)                  :: self
      integer(aq_int), intent(out)                :: nx, ny, nz 
      real(aq_real), intent(out)                  :: deltax, deltay
      real(aq_real), intent(out)                  :: x0, y0
      integer(aq_int), intent(out)                :: halo
      character(len=:), allocatable , intent(out) :: domain
      character(len=:), allocatable , intent(out) :: orientation
      character(len=:), allocatable , intent(out) :: model
      !
      nx = self%nx
      ny = self%ny
      nz = self%nz
      x0 = self%x0
      y0 = self%y0
      deltax = self%deltax
      deltay = self%deltay
      halo = self%halo
      domain = self%domain
      orientation = self%orientation
      model = self%model
      !
   end subroutine aq_geom_info

end module qg_geom_mod
