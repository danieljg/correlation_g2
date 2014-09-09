include 'mkl_vsl.f90'
module vsl_stream
use MKL_VSL_TYPE
use MKL_VSL
implicit none
integer, parameter ::  brng       = VSL_BRNG_MT19937,                          &
                       method_exp = VSL_RNG_METHOD_EXPONENTIAL_ICDF_ACCURATE,  &
                       method_uni = VSL_RNG_METHOD_UNIFORM_STD_ACCURATE,       &
                       method_gauss = VSL_RNG_METHOD_GAUSSIAN_ICDF
integer :: seed(1),seed_size
integer(kind=4) errcode
type (vsl_stream_state) :: stream,stream_gauss
contains
subroutine initialize_stream()
implicit none
 call random_seed
 call random_seed(size=seed_size)
 call random_seed(get=seed)
 errcode=vslnewstream( stream,brng,seed(1) )
 call random_seed
 call random_seed(size=seed_size)
 call random_seed(get=seed)
 errcode=vslnewstream( stream_gauss,brng,seed(1) )
end subroutine initialize_stream
subroutine deinitialize_stream()
implicit none
 errcode=vsldeletestream( stream )
 errcode=vsldeletestream( stream_gauss )
end subroutine
end module vsl_stream

subroutine generate_random_photons_uniform( k_lower_limit, k_upper_limit,              &
         aperture_half_angle, length_k, polar_angle, azimuth_angle, ll, k_degen)
use vsl_stream
implicit none
 real, intent(in)          :: k_lower_limit, k_upper_limit, k_degen,           &
                              aperture_half_angle
 integer, intent(in)       :: ll
 real(kind=8), intent(out) :: length_k(ll+1), polar_angle(ll+1), azimuth_angle(ll+1)
 real(kind=8)    :: r(ll+1) !buffer for random numbers
 real(kind=8)    :: a, b  !limits of uniform distribution
 integer(kind=4) :: i,j
 integer :: n

 n=ll+1; r(:)=0.0
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
 polar_angle(ll+1) = 0.0
 length_k(ll+1)    = k_degen

end subroutine generate_random_photons_uniform

subroutine generate_random_photons_gaussian( k_lower_limit, k_upper_limit,      &
         aperture_half_angle, length_k, polar_angle, azimuth_angle, erf_factor, &
         ll, k_degen, spectral_bandwidth, length_k_uni, temp)
use vsl_stream
use vars_and_funcs
implicit none
 real, intent(in)          :: k_lower_limit, k_upper_limit, k_degen,           &
                              aperture_half_angle, spectral_bandwidth, temp
 integer, intent(in)       :: ll
 real(kind=8), intent(in)  :: length_k_uni(ll+1)
 real(kind=8), intent(out) :: length_k(ll+1), polar_angle(ll+1), azimuth_angle(ll+1),&
                              erf_factor(ll+1)
 real(kind=8)    :: r(ll+1) !buffer for random numbers
 real(kind=8)    :: a, b  !limits of uniform distribution
 integer(kind=4) :: i,j
 integer :: n
 real :: omega_gauss, lambda_gauss

 n=ll+1; r(:)=0.0
 a = 0.0
 b = spectral_bandwidth/sqrt(2.0)
 errcode= vdrnggaussian( method_gauss, stream_gauss, n, r, a, b )
 do j=1, ll
  omega_gauss  = 1.0*omega_pump-omega(length_k_uni(j),temp)+r(j)
  lambda_gauss = 2.0*pi*c/omega_gauss
  length_k(j)  = ktp_index(lambda_gauss,temp)*omega_gauss/c
  if( (length_k(j).gt.k_upper_limit) .or. (length_k(j).lt.k_lower_limit) )then
   omega_gauss  = 1.0*omega_pump-omega(length_k_uni(j),temp)-r(j)
   lambda_gauss = 2.0*pi*c/omega_gauss
   length_k(j)  = ktp_index(lambda_gauss,temp)*omega_gauss/c
  endif
 end do
 r(:)=0.0; a = 0.0; b=8.0*atan(1.0)
 errcode= vdrnguniform( method_uni, stream, n, r, a, b )
 azimuth_angle(:) = r(:)
 r(:)=0.0; a = 0.0; b=1.0;
 errcode= vdrnguniform( method_uni, stream, n, r, a, b )
 polar_angle(:) = acos(1.0+r(:)*(cos(aperture_half_angle)-1.0))

 !fix the last element to match the central wavelength and aperture center
 polar_angle(ll+1) = 0.0
 length_k(ll+1)    = k_degen

 !build error-function factor for importance sampling
 do j=1,ll+1
  erf_factor(j) = erf((omega(length_k_uni(j),temp)-omega(k_lower_limit,temp))/spectral_bandwidth)&
                + erf((omega(k_upper_limit,temp)-omega(length_k_uni(j),temp)/spectral_bandwidth))
 end do

end subroutine generate_random_photons_gaussian
