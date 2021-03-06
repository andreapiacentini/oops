! (C) Copyright 2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module qg_getvalues_mod

use atlas_module, only: atlas_field
use datetime_mod
use fckit_log_module,only: fckit_log
use iso_c_binding
use kinds
use oops_variables_mod
use qg_fields_mod
use qg_geom_mod
use qg_gom_mod
use qg_locs_mod
use aq_constants_mod

!AQ interpolator
use interp_matrix_structure_mod, only : csr_format
use space_time_operator_mod, only : observ_operator
use matrix_manipulations, only : multiply_matrix_csr_vector, addmult_matrixt_csr_vector


implicit none

private
public :: qg_hmat_registry
public :: qg_getvalues_interp, qg_getvalues_build, qg_getvalues_interp_tl, qg_getvalues_interp_ad

#define LISTED_TYPE csr_format

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: qg_hmat_registry
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
! Public
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "oops/util/linkedList_c.f"
! ------------------------------------------------------------------------------
!> Build an interpolation matrix from locations and fields
subroutine qg_getvalues_build(locs,fld,t1,t2,hmat)

implicit none

! Passed variables
type(qg_locs), intent(in) :: locs      !< Locations
type(qg_fields),intent(in) :: fld      !< Fields
type(datetime),intent(in) :: t1, t2    !< times
type(csr_format),intent(inout) :: hmat !< Interpolation matrix

! Local variables
integer :: jloc
real(kind_real), pointer :: lonlat(:,:)
type(atlas_field) :: lonlat_field

!AQ interpolator
integer :: nlev = 1
character(len=aq_strlen) :: msg
integer :: ib
real(kind_real), allocatable :: lonmod(:), latmod(:)
real(kind_real), dimension(1,1) :: dummylev
real(kind_real), dimension(1) :: dummycoord
real(kind_real), dimension(:), allocatable :: dummytime

! Get locations
lonlat_field = locs%lonlat()
call lonlat_field%data(lonlat)

if (trim(fld%geom%orientation) == 'down') nlev = fld%geom%levels

if (fld%geom%fmpi%rank() == 0) then
  !AQ could be stored once for all in the geometry so to avoid the extraction at every build.
  allocate(lonmod(fld%geom%grid%nx(1)))
  allocate(latmod(fld%geom%grid%ny()))
  do ib = 1, fld%geom%grid%nx(1)
     lonmod(ib) = fld%geom%grid%x(ib,1)
  end do
  do ib = 1, fld%geom%grid%ny()
     latmod(ib) = fld%geom%grid%y(ib)
  end do
  allocate(dummytime(locs%nlocs()))
  dummytime(:) = 0_kind_real

  do jloc = 1, locs%nlocs()
     if ( locs%times(jloc) < t1 .or. locs%times(jloc) > t2 ) then
        write(msg,'(A,I0,A)') &
           & 'WARNING from qg_getvalues_build: location nr. ',jloc,' does not fit in the current time interval'
        call fckit_log%warning(msg)
     end if
  end do
  
  call observ_operator ( &
     &   fld%geom%grid%ny(), &
     &   fld%geom%grid%nx(1), &
     &   1, & ! Only on input level in the gathered surface field
     &   1, & ! Only one exact time (no time interpolation)
     &   latmod, &
     &   lonmod, &
     &   dummycoord, & ! Vert coord not relevat
     &   dummycoord, & ! Obs time not relevant
     &   locs%nlocs(), &
     &   1, & ! Levels in and out are 1
     &   1, &
     &   lonlat(2,:), &
     &   lonlat(1,:), &
     &   .false., & ! aq grid not considered as lon periodic
     &   dummytime, &
     &   dummylev, &
     &   2, & ! It is the ground interpolator option
     &   .false., & ! no input scaling
     &   .false., & ! no averaging kernel
     &   .false., & ! no time interpolation
     &   Hmat)

  write(msg,'(A)') 'Built interpolator'
  call fckit_log%debug(msg)
  
  deallocate(lonmod)
  deallocate(latmod)
  deallocate(dummytime)
endif
!Release memory
call lonlat_field%final()

end subroutine qg_getvalues_build
! ------------------------------------------------------------------------------
!> Interpolation from fields
subroutine qg_getvalues_interp(locs,fld,t1,t2,gom)

implicit none

! Passed variables
type(qg_locs), intent(in) :: locs      !< Locations
type(qg_fields),intent(in) :: fld      !< Fields
type(datetime),intent(in) :: t1, t2    !< times
type(qg_gom),intent(inout) :: gom      !< Interpolated values

! Local variables
real(kind_real), pointer :: lonlat(:,:)
type(atlas_field) :: lonlat_field

!AQ interpolator
integer :: nlev = 1
real(aq_real), allocatable, dimension(:,:) :: surf_fld
integer       :: n_vars
character(len=:), allocatable :: var_name(:)
character(len=aq_strlen) :: msg
type(csr_format) :: Hmat

! Get locations
lonlat_field = locs%lonlat()
call lonlat_field%data(lonlat)

if (trim(fld%geom%orientation) == 'down') nlev = fld%geom%levels

n_vars = gom%vars%nvars()
if (n_vars.ne.1) call abor1_ftn('Getvalues interpolates only one field (variable)')

var_name = gom%vars%varlist()
allocate(surf_fld(fld%geom%grid%nx(1),fld%geom%grid%ny()))
call fld%gather_var_at_lev(trim(var_name(1)), nlev, surf_fld, 0)

if (fld%geom%fmpi%rank() == 0) then
  write(msg,'(3A,I2,2(A,G16.8))') 'Interpolate  ',trim(var_name(1)),' at lev ',nlev,' min fld', minval(surf_fld),' max fld', maxval(surf_fld)
  call fckit_log%debug(msg)

  call qg_getvalues_build(locs,fld,t1,t2,hmat)
   
  call multiply_matrix_csr_vector( &
     &   Hmat, &
     &   pack(surf_fld,.true.), &
     &   1, &
     &   locs%nlocs(), &
     &   gom%x)

  write(msg,'(3A,I2,2(A,G16.8))') 'Interpolated ',trim(var_name(1)),' at lev ',nlev,' min Hx', minval(gom%x),' max Hx', maxval(gom%x)
  call fckit_log%debug(msg)
  
endif

call fld%geom%fmpi%broadcast(gom%x,root=0)

!Release memory
deallocate(surf_fld)
call lonlat_field%final()

end subroutine qg_getvalues_interp
! ------------------------------------------------------------------------------
!> Interpolation from fields - tangent linear
subroutine qg_getvalues_interp_tl(locs,fld,t1,t2,hmat,gom)

implicit none

! Passed variables
type(qg_locs), intent(in) :: locs      !< Locations
type(qg_fields),intent(in) :: fld      !< Fields
type(datetime),intent(in) :: t1, t2    !< times
type(csr_format),intent(in) :: hmat    !< Interpolation matrix
type(qg_gom),intent(inout) :: gom      !< Interpolated values

!AQ interpolator
integer :: nlev = 1
real(aq_real), allocatable, dimension(:,:) :: surf_fld
integer       :: n_vars
character(len=:), allocatable :: var_name(:)
character(len=aq_strlen) :: msg

if (trim(fld%geom%orientation) == 'down') nlev = fld%geom%levels

n_vars = gom%vars%nvars()
if (n_vars.ne.1) call abor1_ftn('Getvalues interpolates only one field (variable)')

var_name = gom%vars%varlist()
allocate(surf_fld(fld%geom%grid%nx(1),fld%geom%grid%ny()))
call fld%gather_var_at_lev(trim(var_name(1)), nlev, surf_fld, 0)

if (fld%geom%fmpi%rank() == 0) then
  
  call multiply_matrix_csr_vector( &
     &   Hmat, &
     &   pack(surf_fld,.true.), &
     &   1, &
     &   locs%nlocs(), &
     &   gom%x)

  write(msg,'(3A,I2,2(A,G16.8))') 'Linear interp of ',trim(var_name(1)),' at lev ',nlev,' min Hx', minval(gom%x),' max Hx', maxval(gom%x)
  call fckit_log%debug(msg)
endif

call fld%geom%fmpi%broadcast(gom%x,root=0)

! Release memory
deallocate(surf_fld)

end subroutine qg_getvalues_interp_tl
! ------------------------------------------------------------------------------
!> Interpolation from fields - adjoint
subroutine qg_getvalues_interp_ad(locs,fld,t1,t2,hmat,gom)

implicit none

! Passed variables
type(qg_locs), intent(in) :: locs      !< Locations
type(qg_fields),intent(inout) :: fld   !< Fields
type(datetime),intent(in) :: t1, t2    !< times
type(csr_format),intent(in) :: hmat    !< Interpolation matrix
type(qg_gom),intent(in) :: gom         !< Interpolated values

!AQ interpolator
integer :: nlev = 1
real(aq_real), allocatable, dimension(:) :: surf_1d
real(aq_real), allocatable, dimension(:,:) :: surf_fld
integer       :: n_vars
character(len=:), allocatable :: var_name(:)
character(len=aq_strlen) :: msg

if (trim(fld%geom%orientation) == 'down') nlev = fld%geom%levels

n_vars = gom%vars%nvars()
if (n_vars.ne.1) call abor1_ftn('Getvalues interpolates only one field (variable)')

var_name = gom%vars%varlist()
allocate(surf_fld(fld%geom%grid%nx(1),fld%geom%grid%ny()))
surf_fld(:,:) = 0_kind_real

if (fld%geom%fmpi%rank() == 0) then

   allocate(surf_1d(fld%geom%grid%nx(1)*fld%geom%grid%ny()))
   surf_1d = 0_kind_real
   call addmult_matrixt_csr_vector( &
      &   Hmat, &
      &   gom%x, &
      &   1, &
      &   locs%nlocs(), &
      &   surf_1d)
   surf_fld = unpack(surf_1d,surf_fld==0_kind_real,surf_fld)
   deallocate(surf_1d)

   write(msg,'(3A,I2,2(A,G16.8))') 'Adjoint interp of ',trim(var_name(1)),' at lev ',nlev,' min dx', minval(surf_fld),' max dx', maxval(surf_fld)
  call fckit_log%debug(msg)
endif

call fld%scatteradd_var_at_lev(trim(var_name(1)), nlev, surf_fld, 0)

! Release memory
deallocate(surf_fld)

end subroutine qg_getvalues_interp_ad
! ------------------------------------------------------------------------------
end module qg_getvalues_mod
