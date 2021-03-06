! (C) Copyright 2009-2016 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module qg_fields_interface

use atlas_module, only: atlas_fieldset
use datetime_mod
use fckit_configuration_module, only: fckit_configuration
use iso_c_binding
use kinds
use oops_variables_mod
use qg_fields_mod
use qg_geom_mod
use qg_geom_interface, only: qg_geom_registry
use qg_geom_iter_mod
!AP use qg_locs_mod

implicit none

private
public :: qg_fields_registry

#define LISTED_TYPE qg_fields

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: qg_fields_registry
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

!> Linked list implementation
#include "oops/util/linkedList_c.f"

!> Create fields from geometry and variables
subroutine qg_fields_create_c(c_key_self,c_key_geom,c_vars) bind(c,name='qg_fields_create_f90')

implicit none

! Passed variables
integer(c_int),intent(inout) :: c_key_self       !< Fields
integer(c_int),intent(in) :: c_key_geom          !< Geometry
type(c_ptr),value,intent(in) :: c_vars           !< List of variables

! Local variables
type(qg_fields),pointer :: self
type(qg_geom),pointer :: geom
type(oops_variables) :: vars

! Interface
call qg_fields_registry%init()
call qg_fields_registry%add(c_key_self)
call qg_fields_registry%get(c_key_self,self)
call qg_geom_registry%get(c_key_geom,geom)
vars = oops_variables(c_vars)

! Call Fortran
call self%create(geom,vars)

end subroutine qg_fields_create_c
! ------------------------------------------------------------------------------
!> Create fields from another one
subroutine qg_fields_create_from_other_c(c_key_self,c_key_other) bind(c,name='qg_fields_create_from_other_f90')

implicit none

! Passed variables
integer(c_int),intent(inout) :: c_key_self  !< Fields
integer(c_int),intent(in)    :: c_key_other !< Other fields

! Local variables
type(qg_fields),pointer :: self
type(qg_fields),pointer :: other

! Interface
call qg_fields_registry%get(c_key_other,other)
call qg_fields_registry%init()
call qg_fields_registry%add(c_key_self)
call qg_fields_registry%get(c_key_self,self)

! Call Fortran
call self%create_from(other)

end subroutine qg_fields_create_from_other_c
! ------------------------------------------------------------------------------
!> Delete fields
subroutine qg_fields_delete_c(c_key_self) bind(c,name='qg_fields_delete_f90')

implicit none

! Passed variables
integer(c_int),intent(inout) :: c_key_self !< Fields

! Local variables
type(qg_fields),pointer :: self

! Interface
call qg_fields_registry%get(c_key_self,self)

! Call Fortran
call self%delete()

! Clear interface
call qg_fields_registry%remove(c_key_self)

end subroutine qg_fields_delete_c
! ------------------------------------------------------------------------------
!> Set fields to zero
subroutine qg_fields_zero_c(c_key_self) bind(c,name='qg_fields_zero_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< Fields

! Local variables
type(qg_fields),pointer :: self

! Interface
call qg_fields_registry%get(c_key_self,self)

! Call Fortran
call self%zero()

end subroutine qg_fields_zero_c
! ------------------------------------------------------------------------------
!> Set fields to ones
subroutine qg_fields_ones_c(c_key_self) bind(c,name='qg_fields_ones_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< Fields

! Local variables
type(qg_fields),pointer :: self

! Interface
call qg_fields_registry%get(c_key_self,self)

! Call Fortran
call self%ones()

end subroutine qg_fields_ones_c
! ------------------------------------------------------------------------------
!> Set fields to Diracs
subroutine qg_fields_dirac_c(c_key_self,c_conf) bind(c,name='qg_fields_dirac_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< Fields
type(c_ptr),value,intent(in) :: c_conf  !< Configuration

! Local variables
type(fckit_configuration) :: f_conf
type(qg_fields),pointer :: self

! Interface
f_conf = fckit_configuration(c_conf)
call qg_fields_registry%get(c_key_self,self)

! Call Fortran
call self%dirac(f_conf)

end subroutine qg_fields_dirac_c
! ------------------------------------------------------------------------------
!> Generate random fields
subroutine qg_fields_random_c(c_key_self) bind(c,name='qg_fields_random_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< Fields

! Local variables
type(qg_fields),pointer :: self

! Interface
call qg_fields_registry%get(c_key_self,self)

! Call Fortran
call self%random()

end subroutine qg_fields_random_c
! ------------------------------------------------------------------------------
!> Copy fields
subroutine qg_fields_copy_c(c_key_self,c_key_other) bind(c,name='qg_fields_copy_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self  !< Fields
integer(c_int),intent(in) :: c_key_other !< Other fields

! Local variables
type(qg_fields),pointer :: self
type(qg_fields),pointer :: other

! Interface
call qg_fields_registry%get(c_key_self,self)
call qg_fields_registry%get(c_key_other,other)

! Call Fortran
call self%copy(other)

end subroutine qg_fields_copy_c
! ------------------------------------------------------------------------------
!> Add fields
subroutine qg_fields_self_add_c(c_key_self,c_key_rhs) bind(c,name='qg_fields_self_add_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< Fields
integer(c_int),intent(in) :: c_key_rhs  !< Right-hand side

! Local variables
type(qg_fields),pointer :: self
type(qg_fields),pointer :: rhs

! Interface
call qg_fields_registry%get(c_key_self,self)
call qg_fields_registry%get(c_key_rhs,rhs)

! Call Fortran
call self%self_add(rhs)

end subroutine qg_fields_self_add_c
! ------------------------------------------------------------------------------
!> Subtract fields
subroutine qg_fields_self_sub_c(c_key_self,c_key_rhs) bind(c,name='qg_fields_self_sub_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< Fields
integer(c_int),intent(in) :: c_key_rhs  !< Right-hand side

! Local variables
type(qg_fields),pointer :: self
type(qg_fields),pointer :: rhs

! Interface
call qg_fields_registry%get(c_key_self,self)
call qg_fields_registry%get(c_key_rhs,rhs)

! Call Fortran
call self%self_sub(rhs)

end subroutine qg_fields_self_sub_c
! ------------------------------------------------------------------------------
!> Multiply fields by a scalar
subroutine qg_fields_self_mul_c(c_key_self,c_zz) bind(c,name='qg_fields_self_mul_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< Fields
real(c_double),intent(in) :: c_zz       !< Multiplier

! Local variables
type(qg_fields),pointer :: self

! Interface
call qg_fields_registry%get(c_key_self,self)

! Call Fortran
call self%self_mul(c_zz)

end subroutine qg_fields_self_mul_c
! ------------------------------------------------------------------------------
!> Apply axpy operator to fields
subroutine qg_fields_axpy_c(c_key_self,c_zz,c_key_rhs) bind(c,name='qg_fields_axpy_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< Fields
real(c_double),intent(in) :: c_zz       !< Multiplier
integer(c_int),intent(in) :: c_key_rhs  !< Right-hand side

! Local variables
type(qg_fields),pointer :: self
type(qg_fields),pointer :: rhs

! Interface
call qg_fields_registry%get(c_key_self,self)
call qg_fields_registry%get(c_key_rhs,rhs)

! Call Fortran
call self%axpy(c_zz,rhs)

end subroutine qg_fields_axpy_c
! ------------------------------------------------------------------------------
!> Schur product of fields
subroutine qg_fields_self_schur_c(c_key_self,c_key_rhs) bind(c,name='qg_fields_self_schur_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< Fields
integer(c_int),intent(in) :: c_key_rhs  !< Right-hand side

! Local variables
type(qg_fields),pointer :: self
type(qg_fields),pointer :: rhs

! Interface
call qg_fields_registry%get(c_key_self,self)
call qg_fields_registry%get(c_key_rhs,rhs)

! Call Fortran
call self%schur(rhs)

end subroutine qg_fields_self_schur_c
! ------------------------------------------------------------------------------
!> Compute dot product for fields
subroutine qg_fields_dot_prod_c(c_key_fld1,c_key_fld2,c_prod) bind(c,name='qg_fields_dot_prod_f90')

implicit none

! Passed variables
integer(c_int),intent(in)    :: c_key_fld1 !< First fields
integer(c_int),intent(in)    :: c_key_fld2 !< Second fields
real(c_double),intent(inout) :: c_prod     !< Dot product

! Local variables
type(qg_fields),pointer :: fld1,fld2

! Interface
call qg_fields_registry%get(c_key_fld1,fld1)
call qg_fields_registry%get(c_key_fld2,fld2)

! Call Fortran
call fld1%dot_prod_with(fld2,c_prod)

end subroutine qg_fields_dot_prod_c
! ------------------------------------------------------------------------------
!> Add increment to fields
subroutine qg_fields_add_incr_c(c_key_self,c_key_rhs) bind(c,name='qg_fields_add_incr_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< Fields
integer(c_int),intent(in) :: c_key_rhs  !< Right-hand side

! Local variables
type(qg_fields),pointer :: self
type(qg_fields),pointer :: rhs

! Interface
call qg_fields_registry%get(c_key_self,self)
call qg_fields_registry%get(c_key_rhs,rhs)

! Call Fortran
call self%add_incr(rhs)

end subroutine qg_fields_add_incr_c
! ------------------------------------------------------------------------------
!> Compute increment from the difference of two fields
subroutine qg_fields_diff_incr_c(c_key_lhs,c_key_fld1,c_key_fld2) bind(c,name='qg_fields_diff_incr_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_lhs  !< Left-hand side
integer(c_int),intent(in) :: c_key_fld1 !< First fields
integer(c_int),intent(in) :: c_key_fld2 !< Second fields

! Local variables
type(qg_fields),pointer :: lhs
type(qg_fields),pointer :: fld1
type(qg_fields),pointer :: fld2

! Interface
call qg_fields_registry%get(c_key_lhs,lhs)
call qg_fields_registry%get(c_key_fld1,fld1)
call qg_fields_registry%get(c_key_fld2,fld2)

! Call Fortran
call lhs%diff_incr(fld1,fld2)

end subroutine qg_fields_diff_incr_c
! ------------------------------------------------------------------------------
! !> Change fields resolution
! subroutine qg_fields_change_resol_c(c_key_fld,c_key_rhs) bind(c,name='qg_fields_change_resol_f90')

! implicit none

! ! Passed variables
! integer(c_int),intent(in) :: c_key_fld !< Fields
! integer(c_int),intent(in) :: c_key_rhs !< Right-hand side

! ! Local variables
! type(qg_fields),pointer :: fld,rhs

! ! Interface
! call qg_fields_registry%get(c_key_fld,fld)
! call qg_fields_registry%get(c_key_rhs,rhs)

! ! Call Fortran
! call qg_fields_change_resol(fld,rhs)

! end subroutine qg_fields_change_resol_c
! ------------------------------------------------------------------------------
!> Get fields info
subroutine qg_fields_info_c(c_key_fld,c_conf) bind(c,name='qg_fields_info_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_fld  !< Fields
type(c_ptr),value,intent(in) :: c_conf  !< Configuration

! Local variables
type(fckit_configuration) :: f_conf
type(qg_fields),pointer :: fld

! Interface
f_conf = fckit_configuration(c_conf)
call qg_fields_registry%get(c_key_fld,fld)

! Call Fortran
call fld%info(f_conf)

end subroutine qg_fields_info_c
! ------------------------------------------------------------------------------
!> Read fields from file
subroutine qg_fields_read_file_c(c_key_fld,c_conf,c_dt) bind(c,name='qg_fields_read_file_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_fld  !< Fields
type(c_ptr),value,intent(in) :: c_conf  !< Configuration
type(c_ptr),value,intent(in) :: c_dt    !< Date and time

! Local variables
type(fckit_configuration) :: f_conf
type(qg_fields),pointer :: fld
type(datetime) :: fdate

! Interface
f_conf = fckit_configuration(c_conf)
call qg_fields_registry%get(c_key_fld,fld)
call c_f_datetime(c_dt,fdate)

! Call Fortran
call fld%read(f_conf,fdate)

end subroutine qg_fields_read_file_c
! ------------------------------------------------------------------------------
!> Write fields to file
subroutine qg_fields_write_file_c(c_key_fld,c_conf,c_dt) bind(c,name='qg_fields_write_file_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_fld !< Fields
type(c_ptr),value,intent(in) :: c_conf !< Configuration
type(c_ptr),value,intent(in) :: c_dt   !< Date and time

! Local variables
type(fckit_configuration) :: f_conf
type(qg_fields),pointer :: fld
type(datetime) :: fdate

! Interface
f_conf = fckit_configuration(c_conf)
call qg_fields_registry%get(c_key_fld,fld)
call c_f_datetime(c_dt,fdate)

! Call Fortran
call fld%write(f_conf,fdate)

end subroutine qg_fields_write_file_c
! ------------------------------------------------------------------------------
!> Analytic initialization of fields
subroutine qg_fields_analytic_init_c(c_key_fld,c_conf) bind(c,name='qg_fields_analytic_init_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_fld  !< Fields
type(c_ptr),value,intent(in) :: c_conf  !< Configuration

! Local variables
type(fckit_configuration) :: f_conf
type(qg_fields),pointer :: fld

! Interface
f_conf = fckit_configuration(c_conf)
call qg_fields_registry%get(c_key_fld,fld)

! Call Fortran
 call fld%analytic_IC()

end subroutine qg_fields_analytic_init_c
! ------------------------------------------------------------------------------
!> Fields statistics
subroutine qg_fields_gpnorm_c(c_key_fld,vpresent,vmin,vmax,vrms) bind(c,name='qg_fields_gpnorm_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_fld      !< Fields
integer(c_int),intent(inout) :: vpresent(6) !< Variables presence flag
real(c_double),intent(inout) :: vmin(6)     !< Variables minimum
real(c_double),intent(inout) :: vmax(6)     !< Variables maximum
real(c_double),intent(inout) :: vrms(6)     !< Variables RMS

! Local variables
type(qg_fields),pointer :: fld

! Interface
call qg_fields_registry%get(c_key_fld,fld)

! Call Fortran
! call fld%stats(vpresent,vmin,vmax,vrms)

end subroutine qg_fields_gpnorm_c
! ------------------------------------------------------------------------------
!> Fields RMS
subroutine qg_fields_rms_c(c_key_fld,prms) bind(c,name='qg_fields_rms_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_fld !< Fields
real(c_double),intent(inout) :: prms

! Local variables
type(qg_fields),pointer :: fld

! Interface
call qg_fields_registry%get(c_key_fld,fld)

! Call Fortran
call fld%norm(prms)

end subroutine qg_fields_rms_c
! ------------------------------------------------------------------------------
!> Create ATLAS fields
subroutine qg_fields_set_atlas_c(c_key_fld,c_vars,c_afieldset) bind (c,name='qg_fields_set_atlas_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_fld           !< Fields
type(c_ptr),value,intent(in) :: c_vars           !< List of variables
type(c_ptr),intent(in),value :: c_afieldset      !< ATLAS fieldset pointer

! Local variables
type(qg_fields),pointer :: fld
type(oops_variables) :: vars
type(atlas_fieldset) :: afieldset

! Interface
call qg_fields_registry%get(c_key_fld,fld)
vars = oops_variables(c_vars)
afieldset = atlas_fieldset(c_afieldset)

! Call Fortran
! call qg_fields_set_atlas(fld,vars,afieldset)

end subroutine qg_fields_set_atlas_c
! ------------------------------------------------------------------------------
!> Convert fields to ATLAS
subroutine qg_fields_to_atlas_c(c_key_fld,c_vars,c_afieldset) bind (c,name='qg_fields_to_atlas_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_fld           !< Fields
type(c_ptr),value,intent(in) :: c_vars           !< List of variables
type(c_ptr),intent(in),value :: c_afieldset      !< ATLAS fieldset pointer

! Local variables
type(qg_fields),pointer :: fld
type(oops_variables) :: vars
type(atlas_fieldset) :: afieldset

! Interface
call qg_fields_registry%get(c_key_fld,fld)
vars = oops_variables(c_vars)
afieldset = atlas_fieldset(c_afieldset)

! Call Fortran
! call qg_fields_to_atlas(fld,vars,afieldset)

end subroutine qg_fields_to_atlas_c
! ------------------------------------------------------------------------------
!> Get fields from ATLAS
subroutine qg_fields_from_atlas_c(c_key_fld,c_vars,c_afieldset) bind (c,name='qg_fields_from_atlas_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_fld           !< Fields
type(c_ptr),value,intent(in) :: c_vars           !< List of variables
type(c_ptr),intent(in),value :: c_afieldset      !< ATLAS fieldset pointer

! Local variables
type(qg_fields),pointer :: fld
type(oops_variables) :: vars
type(atlas_fieldset) :: afieldset

! Interface
call qg_fields_registry%get(c_key_fld,fld)
vars = oops_variables(c_vars)
afieldset = atlas_fieldset(c_afieldset)

! Call Fortran
! call qg_fields_from_atlas(fld,vars,afieldset)

end subroutine qg_fields_from_atlas_c
! ------------------------------------------------------------------------------
!> Get points from fields
subroutine qg_fields_getpoint_c(c_key_fld,c_key_iter,c_nval,c_vals) bind(c,name='qg_fields_getpoint_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_fld         !< Fields
integer(c_int),intent(in) :: c_key_iter        !< Geometry iterator
integer(c_int),intent(in) :: c_nval            !< Number of values
real(c_double),intent(inout) :: c_vals(c_nval) !< Values

! Local variables
type(qg_fields),pointer :: fld
type(qg_geom_iter),pointer :: iter

! Interface
call qg_fields_registry%get(c_key_fld,fld)
call qg_geom_iter_registry%get(c_key_iter,iter)

! Call Fortran
! call qg_fields_getpoint(fld,iter,c_nval,c_vals)

end subroutine qg_fields_getpoint_c
! ------------------------------------------------------------------------------
!> Set points for the fields
subroutine qg_fields_setpoint_c(c_key_fld,c_key_iter,c_nval,c_vals) bind(c,name='qg_fields_setpoint_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_fld         !< Fields
integer(c_int),intent(in) :: c_key_iter        !< Geometry iterator
integer(c_int),intent(in) :: c_nval            !< Number of values
real(c_double),intent(in) :: c_vals(c_nval)    !< Values

! Local variables
type(qg_fields),pointer :: fld
type(qg_geom_iter),pointer :: iter

! Interface
call qg_fields_registry%get(c_key_fld,fld)
call qg_geom_iter_registry%get(c_key_iter,iter)

! Call Fortran
! call qg_fields_setpoint(fld,iter,c_nval,c_vals)

end subroutine qg_fields_setpoint_c
! ------------------------------------------------------------------------------
!> Serial Size
subroutine qg_fields_serialsize_c(c_key_fld,c_vsize) bind(c,name='qg_fields_serialsize_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_fld            !< Fields
integer(c_int),intent(out) :: c_vsize              !< Size

! Local variables
type(qg_fields),pointer :: fld

! Interface
call qg_fields_registry%get(c_key_fld,fld)

! Call Fortran
c_vsize = fld%serialsize()

end subroutine qg_fields_serialsize_c
! ------------------------------------------------------------------------------
!> Serialize fields
subroutine qg_fields_serialize_c(c_key_fld,c_vsize,c_vect_fld) bind(c,name='qg_fields_serialize_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_fld            !< Fields
integer(c_int),intent(in) :: c_vsize              !< Size
real(c_double),intent(out) :: c_vect_fld(c_vsize) !< Vector

! Local variables
type(qg_fields),pointer :: fld

! Interface
call qg_fields_registry%get(c_key_fld,fld)

! Call Fortran
call fld%serialize(c_vect_fld)

end subroutine qg_fields_serialize_c
! ------------------------------------------------------------------------------
!> Deserialize fields
subroutine qg_fields_deserialize_c(c_key_self,c_vsize,c_vect_fld,c_index) bind(c,name='qg_fields_deserialize_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self          !< Fields
integer(c_int),intent(in) :: c_vsize             !< Size
real(c_double),intent(in) :: c_vect_fld(c_vsize) !< Vector
integer(c_int), intent(inout):: c_index          !< Index

! Local variables
type(qg_fields),pointer :: self

! Interface
call qg_fields_registry%get(c_key_self,self)

! Call Fortran
call self%deserialize(c_vect_fld,c_index)

end subroutine qg_fields_deserialize_c
! ------------------------------------------------------------------------------
end module qg_fields_interface
