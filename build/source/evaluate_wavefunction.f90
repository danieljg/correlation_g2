subroutine evaluate_wavefunction(signal_k, signal_kx, signal_ky, signal_kz,    &
                            idler_k,  idler_kx,  idler_ky,  idler_kz,          &
                            wavefunction, temp, erf_factor,                    &
                            point_value, phase_mismatch )
use vars_and_funcs
implicit none
 real, intent(in) :: signal_k(nn+1), signal_kx(nn+1), signal_ky(nn+1), signal_kz(nn+1),&
                     idler_k(nn+1),  idler_kx(nn+1),  idler_ky(nn+1),  idler_kz(nn+1), &
                     erf_factor(nn+1), temp
 real, intent(out) :: wavefunction(nn+1), point_value, phase_mismatch
 real :: delta_ksq(nn+1), sinc_vec(nn+1), omega_sum(nn+1)!, pump_profile(nn+1)
 integer i
 ! build vectors for wavefunction evaluation
 delta_ksq(:) = (signal_kx(:) + idler_kx(:))**2.0                              &
              + (signal_ky(:) + idler_ky(:))**2.0
 do i=1,nn+1
  omega_sum(i)    = omega(signal_k(i),temp)+omega(idler_k(i),temp)
  sinc_vec(i)     = calc_phase_mismatch(omega_sum(i),signal_kz(i),idler_kz(i),temp)
  sinc_vec(i)     = (sinc(sinc_vec(i)))**2
!  pump_profile(i) = exp( -(omega_sum(i)-omega_pump)**2/spectral_width**2 )
 end do
 !asign point-value and phase-mismatch
 point_value = erf_factor(nn+1)*exp( - (beam_waist**2*delta_ksq(nn+1)/2.0) )*  &
               sinc_vec(nn+1)!&*pump_profile(nn+1)! not needed with importance sampling
!               d_eff*2*sqrt(pump_power/2.0)*                                   &
!               beam_waist/(pi*spectral_width)*                                 &
 phase_mismatch = calc_phase_mismatch(omega_pump,signal_kz(nn+1),idler_kz(nn+1),temp)
 !assign the whole pack
 wavefunction(:) = erf_factor(:)*exp( - (beam_waist**2*delta_ksq(:)/2.0) )* &
                   sinc_vec(:)!&*pump_profile(:)! not needed with importance sampling
!               d_eff*2*sqrt(pump_power/2.0)*                                   &
!               beam_waist/(pi*spectral_width)*                                 &

 contains
 real function calc_phase_mismatch(pump_omega, signal_kz, idler_kz, temp )
 real :: pump_omega, signal_kz, idler_kz, temp
 calc_phase_mismatch = k_p(pump_omega) - signal_kz - idler_kz                  &
                     - 2.0*pi/poling_period_temp(temp)
 calc_phase_mismatch = 0.5*crystal_length_temp(temp)*calc_phase_mismatch
 end function
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
