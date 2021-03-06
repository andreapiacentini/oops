! (C) Copyright 2009-2016 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module qg_geom_interface

use atlas_module, only: atlas_fieldset, atlas_functionspace_pointcloud
use fckit_configuration_module, only: fckit_configuration
use fckit_log_module,only: fckit_log
use fckit_mpi_module,only: fckit_mpi_comm
use kinds
use iso_c_binding
!AP use qg_projection_mod
use qg_geom_mod
use aq_constants_mod, only : aq_strlen

implicit none

private
public :: qg_geom_registry

#define LISTED_TYPE qg_geom

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: qg_geom_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

!> Linked list implementation
#include "oops/util/linkedList_c.f"

!> Setup geometry
subroutine qg_geom_setup_c(c_key_self,c_conf,c_comm) bind(c,name='qg_geom_setup_f90')

! Passed variables
integer(c_int),intent(inout) :: c_key_self !< Geometry
type(c_ptr),value,intent(in) :: c_conf     !< Configuration
type(c_ptr), value, intent(in)    :: c_comm

! Local variables
type(fckit_configuration) :: f_conf
type(fckit_mpi_comm)      :: f_comm
type(qg_geom),pointer :: self

! Interface
f_conf = fckit_configuration(c_conf)
f_comm = fckit_mpi_comm(c_comm)
call qg_geom_registry%init()
call qg_geom_registry%add(c_key_self)
call qg_geom_registry%get(c_key_self,self)

! Call Fortran
call self%create(f_conf,f_comm)

end subroutine qg_geom_setup_c
! ------------------------------------------------------------------------------
!> Set ATLAS lon/lat field
! subroutine qg_geom_set_atlas_lonlat_c(c_key_self,c_afieldset) bind(c,name='qg_geom_set_atlas_lonlat_f90')

! ! Passed variables
! integer(c_int),intent(in) :: c_key_self     !< Geometry
! type(c_ptr),intent(in),value :: c_afieldset !< ATLAS fieldset pointer

! ! Local variables
! type(qg_geom),pointer :: self
! type(atlas_fieldset) :: afieldset

! ! Interface
! call qg_geom_registry%get(c_key_self,self)
! afieldset = atlas_fieldset(c_afieldset)

! ! Call Fortran
! call qg_geom_set_atlas_lonlat(self,afieldset)

! end subroutine qg_geom_set_atlas_lonlat_c
! ------------------------------------------------------------------------------
!> Set ATLAS function space pointer
! subroutine qg_geom_set_atlas_functionspace_pointer_c(c_key_self,c_afunctionspace) &
!  & bind(c,name='qg_geom_set_atlas_functionspace_pointer_f90')

! ! Passed variables
! integer(c_int),intent(in) :: c_key_self          !< Geometry
! type(c_ptr),intent(in),value :: c_afunctionspace !< ATLAS function space pointer

! ! Local variables
! type(qg_geom),pointer :: self

! ! Interface
! call qg_geom_registry%get(c_key_self,self)
! self%afunctionspace = atlas_functionspace_pointcloud(c_afunctionspace)

! end subroutine qg_geom_set_atlas_functionspace_pointer_c
! ------------------------------------------------------------------------------
!> Fill ATLAS fieldset
! subroutine qg_geom_fill_atlas_fieldset_c(c_key_self,c_afieldset) &
!  & bind(c,name='qg_geom_fill_atlas_fieldset_f90')

! ! Passed variables
! integer(c_int),intent(in) :: c_key_self     !< Geometry
! type(c_ptr),intent(in),value :: c_afieldset !< ATLAS fieldset pointer

! ! Local variables
! type(qg_geom),pointer :: self
! type(atlas_fieldset) :: afieldset

! ! Interface
! call qg_geom_registry%get(c_key_self,self)
! afieldset = atlas_fieldset(c_afieldset)

! ! Call Fortran
! call qg_geom_fill_atlas_fieldset(self,afieldset)

! end subroutine qg_geom_fill_atlas_fieldset_c
! ------------------------------------------------------------------------------
!> Clone geometry
subroutine qg_geom_clone_c(c_key_self,c_key_other) bind(c,name='qg_geom_clone_f90')

! Passed variables
integer(c_int),intent(inout) :: c_key_self !< Geometry
integer(c_int),intent(in) :: c_key_other   !< Other geometry

! Local variables
type(qg_geom),pointer :: self,other

! Interface
call qg_geom_registry%get(c_key_other,other)
call qg_geom_registry%init()
call qg_geom_registry%add(c_key_self)
call qg_geom_registry%get(c_key_self ,self )

! Call Fortran
call self%clone(other)

end subroutine qg_geom_clone_c
! ------------------------------------------------------------------------------
!> Delete geometry
subroutine qg_geom_delete_c(c_key_self) bind(c,name='qg_geom_delete_f90')

! Passed variables
integer(c_int),intent(inout) :: c_key_self !< Geometry

! Local variables
type(qg_geom),pointer :: self

! Interface
call qg_geom_registry%get(c_key_self,self)

! Call Fortran
call self%delete()

! Clear interface
call qg_geom_registry%remove(c_key_self)

end subroutine qg_geom_delete_c
! ------------------------------------------------------------------------------
!> Get geometry info
subroutine qg_geom_info_c(c_key_self,c_nx,c_ny,c_nz,c_deltax,c_deltay,&
   & c_orientation_ptr, c_domain_ptr, c_model_ptr) bind(c,name='qg_geom_info_f90')

! Passed variables
integer(c_int),intent(in) :: c_key_self  !< Geometry
integer(c_int),intent(inout) :: c_nx     !< Number of points in the zonal direction
integer(c_int),intent(inout) :: c_ny     !< Number of points in the meridional direction
integer(c_int),intent(inout) :: c_nz     !< Number of vertical levels
real(c_double),intent(inout) :: c_deltax !< Zonal cell size
real(c_double),intent(inout) :: c_deltay !< Meridional cell size
character(kind=c_char), target, intent(inout) :: c_orientation_ptr(*)
character(kind=c_char), target, intent(inout) :: c_domain_ptr(*)
character(kind=c_char), target, intent(inout) :: c_model_ptr(*)
! Local variables
integer :: i, length
character(kind=c_char), pointer :: c_orientation(:)
character(kind=c_char), pointer :: c_domain(:)
character(kind=c_char), pointer :: c_model(:)
type(qg_geom),pointer :: self
real(c_double) :: x0, y0
integer(c_int) :: halo
character(len=:), allocatable :: domain, orientation, model

! Interface
call qg_geom_registry%get(c_key_self,self)

! Call Fortran
call self%info(c_nx,c_ny,c_nz,c_deltax,c_deltay,x0, y0, halo, domain, orientation, model)
length = len_trim(orientation)
call c_f_pointer(c_loc(c_orientation_ptr), c_orientation, [aq_strlen])
do i = 1, length
   c_orientation(i)=orientation(i:i)
end do
c_orientation(length+1) = C_NULL_CHAR
length = len_trim(domain)
call c_f_pointer(c_loc(c_domain_ptr), c_domain, [aq_strlen])
do i = 1, length
   c_domain(i)=domain(i:i)
end do
c_domain(length+1) = C_NULL_CHAR
length = len_trim(model)
call c_f_pointer(c_loc(c_model_ptr), c_model, [aq_strlen])
do i = 1, length
   c_model(i)=model(i:i)
end do
c_model(length+1) = C_NULL_CHAR

end subroutine qg_geom_info_c

! ------------------------------------------------------------------------------
end module qg_geom_interface
