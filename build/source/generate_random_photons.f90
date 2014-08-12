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

end subroutine generate_random_photons_uniform

subroutine generate_random_photons_gaussian( k_lower_limit, k_upper_limit,      &
         aperture_half_angle, length_k, polar_angle, azimuth_angle, erf_factor, &
         nn, k_degen, spectral_bandwidth, length_k_uni, temp)
use vsl_stream
use vars_and_funcs, only: omega_pump, pi, c
implicit none
 real, intent(in)          :: k_lower_limit, k_upper_limit, k_degen,           &
                              aperture_half_angle, spectral_bandwidth, temp
 integer, intent(in)       :: nn
 real(kind=8), intent(in)  :: length_k_uni(nn+1)
 real(kind=8), intent(out) :: length_k(nn+1), polar_angle(nn+1), azimuth_angle(nn+1),&
                              erf_factor(nn+1)
 real(kind=8)    :: r(nn+1) !buffer for random numbers
 real(kind=8)    :: a, b  !limits of uniform distribution
 integer(kind=4) :: i,j
 integer :: n
 real :: omega_gauss, lambda_gauss

 n=nn+1; r(:)=0.0
 a = 0.0
 b = spectral_bandwidth/sqrt(2.0)
 errcode= vdrnggaussian( method_gauss, stream_gauss, n, r, a, b )
 do j=1, nn
  omega_gauss  = 1.0*omega_pump-omega(length_k_uni(j))+r(j)
  lambda_gauss = 2.0*pi*c/omega_gauss
  length_k(j)  = ktp_index(lambda_gauss,temp)*omega_gauss/c
  if( (length_k(j).gt.k_upper_limit) .or. (length_k(j).lt.k_lower_limit) )then
   omega_gauss  = 1.0*omega_pump-omega(length_k_uni(j))-r(j)
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
 polar_angle(nn+1) = 0.0
 length_k(nn+1)    = k_degen

 !build error-function factor for importance sampling
 do j=1,nn+1
  erf_factor(j) = erf((omega(length_k_uni(j))-omega(k_lower_limit))/spectral_bandwidth)&
                + erf((omega(k_upper_limit)-omega(length_k_uni(j))/spectral_bandwidth))
 end do

 contains

 real function ktp_index(x,temp)
 real :: x, temp
 real :: lx
 real, parameter :: na=2.12725, nb=1.18431, nc=5.14852e-2,&
                    nd=0.6603, ne=100.00506, nf=9.68956e-3
 real, parameter :: a0=9.9587e-6, a1=9.9228e-6, a2=-8.9603e-6, a3=4.1010e-6
 real, parameter :: b0=-1.1882e-8, b1=10.459e-8, b2=-9.8136e-8, b3=3.1481e-8
 lx = x*1.0e6
 ktp_index = sqrt( na + nb/(1.0-nc/(lx**2))+ nd/(1.0-ne/(lx**2)) - nf*(lx**2.0))&
       + (a0 + a1/lx + a2/(lx**2) + a3/(lx**3))*(temp-25.0)&
       + (b0 + b1/lx + b2/(lx**2) + b3/(lx**3))*(temp-25.0)**2.0
 end function ktp_index

 real function omega(k_len)
 use vars_and_funcs, only: low_wvln_signal,high_wvln_signal,low_wvln_idler,high_wvln_idler
 real :: k_len, y_left, y_right, lambda_left, lambda_right
 real :: m, b, lambda_save, lambda_new, y_new
 integer i
 lambda_left  = min(low_wvln_signal, low_wvln_idler)
 lambda_right = max(high_wvln_signal, high_wvln_idler)
 lambda_save  = lambda_left
 do i=1,16
  y_left  = ktp_index(lambda_left,temp)*2.0*pi/lambda_left - k_len
  y_right = ktp_index(lambda_right,temp)*2.0*pi/lambda_right - k_len
  m = (y_right-y_left)/(lambda_right-lambda_left)
  b = y_left-(m*lambda_left)
  lambda_new = -b/m
  y_new = ktp_index(lambda_new,temp)*2.0*pi/lambda_new - k_len
  if(y_new.gt.0.0)then
   lambda_save  = lambda_left
   lambda_left  = lambda_new
  else
   lambda_save  = lambda_right
   lambda_right = lambda_new
   lambda_left  = lambda_right-0.1*abs(lambda_save-lambda_right)
  endif
  if(abs(lambda_save-lambda_new)/lambda_new.lt.10*lambda_right*epsilon(k_len))exit
  if(lambda_left.eq.lambda_right)exit
 end do
 omega=2.0*pi*c/lambda_new
 end function omega
end subroutine generate_random_photons_gaussian
