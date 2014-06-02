subroutine evaluate_wavefunction(signal_k, signal_kx, signal_ky, signal_kz,    &
                            idler_k,  idler_kx,  idler_ky,  idler_kz,          &
                            wavefunction, temp, k_a, k_b)
use vars_and_funcs
implicit none
 real, intent(in) :: signal_k(nn), signal_kx(nn), signal_ky(nn), signal_kz(nn),&
                     idler_k(nn),  idler_kx(nn),  idler_ky(nn),  idler_kz(nn), &
                     temp, k_a, k_b
 real, intent(out) :: wavefunction(nn)
 real :: delta_kz(nn), delta_ksq(nn), sinc_vec(nn), omega_sum(nn)
 integer i
 delta_ksq(:) = (signal_kx(:) + idler_kx(:))**2.0                              &
              + (signal_ky(:) + idler_ky(:))**2.0
 do i=1,nn
  delta_kz(i) = k_p(omega(signal_k(i)) + omega(idler_k(i)))                    &
              - signal_kz(i) - idler_kz(i) - 2.0*pi/poling_period
  omega_sum(i) = omega(signal_k(i))+omega(idler_k(i))
  sinc_vec(i)  = (sinc(0.5*crystal_length*delta_kz(i)))**2
 end do
 wavefunction(:) = d_eff**2*sqrt(pump_power/2.0)*beam_waist/(pi*spectral_width)  &
                   *exp( -(beam_waist**2*delta_ksq(:)/2.0) )                   &!seems ok
                   *exp( -2.0*(omega_sum(:)-omega_pump)**2                     &
                          /spectral_width**2 )                                 &
                   *sinc_vec(:)
 contains
 real function omega(k_len)
 real :: k_len, y_left, y_right, lambda_left, lambda_right
 real :: m, b, lambda_save, 	lambda_new
 integer i
 lambda_left  = low_wavelength_limit
 lambda_right = high_wavelength_limit
 lambda_save  = lambda_left
 do i=1,16
  y_left  = ktp_index(lambda_left,temp)*2.0*pi/lambda_left - k_len
  y_right = ktp_index(lambda_right,temp)*2.0*pi/lambda_right - k_len
  m = (y_right-y_left)/(lambda_right-lambda_left)
  b = y_left-(m*lambda_left)
  lambda_new = -b/m
  if((ktp_index(lambda_new,temp)*lambda_new/c-k_len).gt.0.0)then
   lambda_save  = lambda_left
   lambda_left  = lambda_new
  else
   lambda_save  = lambda_right
   lambda_right = lambda_new
  endif
  if(abs(lambda_save-lambda_new)/lambda_new.lt.10*lambda_right*epsilon(k_len))exit
 end do
 omega=2.0*pi*c/lambda_new
 end function omega
 real function k_p(omega)
 real omega
  k_p = ktp_index(2.0*pi*c/omega,temp)*omega/c
 end function k_p
 real function sinc(x)
 real x
 real, parameter :: tn=tiny(x)
 sinc=sin(x+tn)/(x+tn)
 end function sinc
end subroutine evaluate_wavefunction
