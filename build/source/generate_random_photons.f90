include 'mkl_vsl.f90'
module vsl_stream
use MKL_VSL_TYPE
use MKL_VSL
implicit none
integer, parameter ::  brng       = VSL_BRNG_MT19937,                          &
                       method_exp = VSL_RNG_METHOD_EXPONENTIAL_ICDF_ACCURATE,  &
                       method_uni = VSL_RNG_METHOD_UNIFORM_STD_ACCURATE
integer :: seed(1),seed_size
integer(kind=4) errcode
type (vsl_stream_state) :: stream
contains
subroutine initialize_stream()
implicit none
 call random_seed
 call random_seed(size=seed_size)
 call random_seed(get=seed)
 errcode=vslnewstream( stream,brng,seed(1) )
end subroutine initialize_stream
subroutine deinitialize_stream()
implicit none
 errcode=vsldeletestream( stream )
end subroutine
end module vsl_stream

subroutine generate_random_photons( k_lower_limit, k_upper_limit,              &
         aperture_half_angle, length_k, polar_angle, azimuth_angle, nn, k_degen)
use vsl_stream
implicit none
 real, intent(in)          :: k_lower_limit, k_upper_limit, k_degen,           &
                              aperture_half_angle
 integer, intent(in)       :: nn
 real(kind=8), intent(out) :: length_k(nn+1), polar_angle(nn+1), azimuth_angle(nn+1)
 real(kind=8)    :: r(nn+1) !buffer for random numbers
 real(kind=8)    :: a, b  !limits of uniform distribution
 integer(kind=4) :: i,j
 integer :: n

 n=nn+1; r(:)=0.0
 a = k_lower_limit
 b = k_upper_limit
 errcode= vdrnguniform( method_uni, stream, n, r, a, b )
 length_k(:) = r(:)
 r(:)=0.0; a = 0.0; b=8.0*atan(1.0)
 errcode= vdrnguniform( method_uni, stream, n, r, a, b )
 azimuth_angle(:) = r(:)
 r(:)=0.0; a = 0.0; b=1.0;
 errcode= vdrnguniform( method_uni, stream, n, r, a, b )
 polar_angle(:) = acos(1.0+r(:)*(cos(aperture_half_angle)-1.0))

 !fix the last element to match the central wavelength and aperture center
 polar_angle(nn+1) = 0.0
 length_k(nn+1)    = k_degen

end subroutine generate_random_photons
