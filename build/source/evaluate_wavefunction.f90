subroutine evaluate_wavefunction(signal_k, signal_kx, signal_ky, signal_kz,    &
                            idler_k,  idler_kx,  idler_ky,  idler_kz,          &
                            wavefunction, temp, k_a_signal, k_b_signal,        &
                            k_a_idler, k_b_idler,                              &
                            point_value, phase_mismatch )
use vars_and_funcs
implicit none
 real, intent(in) :: signal_k(nn+1), signal_kx(nn+1), signal_ky(nn+1), signal_kz(nn+1),&
                     idler_k(nn+1),  idler_kx(nn+1),  idler_ky(nn+1),  idler_kz(nn+1), &
                     temp, k_a_signal, k_b_signal, k_a_idler, k_b_idler
 real, intent(out) :: wavefunction(nn+1), point_value, phase_mismatch
 real :: delta_ksq(nn+1), sinc_vec(nn+1), omega_sum(nn+1), pump_profile(nn+1)
 integer i
 ! build vectors for wavefunction evaluation
 delta_ksq(:) = (signal_kx(:) + idler_kx(:))**2.0                              &
              + (signal_ky(:) + idler_ky(:))**2.0
 do i=1,nn+1
  omega_sum(i)    = omega(signal_k(i))+omega(idler_k(i))
  sinc_vec(i)     = calc_phase_mismatch(omega_sum(i),signal_kz(i),idler_kz(i))
  sinc_vec(i)     = (sinc(sinc_vec(i)))**2
  pump_profile(i) = exp( -(omega_sum(i)-omega_pump)**2/spectral_width**2 )
 end do
 !asign point-value and phase-mismatch
 point_value = d_eff*2*sqrt(pump_power/2.0)*beam_waist/(pi*spectral_width)*    &
               exp( - (beam_waist**2*delta_ksq(nn+1)/2.0) )                    &
               *pump_profile(nn+1)*sinc_vec(nn+1)
 phase_mismatch = calc_phase_mismatch(omega_pump,signal_kz(nn+1),idler_kz(nn+1))
 !assign the whole pack
 wavefunction(:) = d_eff**2*sqrt(pump_power/2.0)*beam_waist/(pi*spectral_width)*&
                   exp( -(beam_waist**2*delta_ksq(:)/2.0) )                    &!seems ok
                   *pump_profile(:)*sinc_vec(:)
 contains
 real function calc_phase_mismatch(pump_omega, signal_kz, idler_kz )
 real :: pump_omega, signal_kz, idler_kz
 calc_phase_mismatch = k_p(pump_omega) - signal_kz - idler_kz                  &
                     - 2.0*pi/poling_period_temp()
 calc_phase_mismatch = 0.5*crystal_length_temp()*calc_phase_mismatch
 end function
 real function omega(k_len)
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
 real function k_p(omega)
 real omega, lambda
  lambda = 2.0*pi*c/omega
  k_p = 2.0*pi*ktp_index(lambda,temp)/lambda
 end function k_p
 real function sinc(x)
 real x
 real, parameter :: tn=tiny(x)
 sinc=sin(x+tn)/(x+tn)
 end function sinc
end subroutine evaluate_wavefunction
