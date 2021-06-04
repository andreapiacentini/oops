module aq_constants_mod
   use kinds
   use mpi
   
   implicit none

   private
   public :: aq_strlen, aq_varlen
   public :: oops_int, oops_long, oops_single, oops_real
   public :: aq_int, aq_long, aq_single, aq_real, aq_double

   integer, parameter :: aq_strlen = 256
   integer, parameter :: aq_varlen = 100 ! as MAXVARLEN in oops_variables_mod

   integer, parameter :: oops_int    = kind_int
   integer, parameter :: oops_long   = kind_long
   integer, parameter :: oops_single = kind_single
   integer, parameter :: oops_real   = kind_real
   
   integer, parameter :: aq_int = kind_int  
   integer, parameter :: aq_long = kind_long
   integer, parameter :: aq_single = kind_single
   integer, parameter :: aq_real = kind_real
   integer, parameter :: aq_double = kind_real
   
end module aq_constants_mod

