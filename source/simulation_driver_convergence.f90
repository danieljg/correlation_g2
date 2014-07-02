module vars_and_funcs
 integer, parameter :: nn = 1e3,                                               &
                       repetitions = 50
 real,    parameter ::                                                         &
   work_temp              = 51.5,                                              &
   axial_pm_temp          = 65.0,                                              &
   crystal_dist           = 0.5,                                               &
   pump_power             = 0.030,                                             &
   pump_wavelength        = 406.118e-9,                                        &
   low_wvln_signal        = 805.0e-9,                                          &
   high_wvln_signal       = 820.0e-9,                                          &
   beam_waist             = 50.0e-6,                   &!!!!GET THESE NUMBERS RIGHT
   spectral_width_nm      = 0.035e-9,                  &!!!!GET THESE NUMBERS RIGHT
   crystal_length         = 0.005,                                             &
   d_eff                  = 1.0,                       &!!!!GET THESE NUMBERS RIGHT
   pi                     = 4.0*atan(1.0),                                     &
   c                      = 3.0e8
 real, parameter :: omega_pump     = 2.0*pi*c/pump_wavelength,                 &
                    spectral_width = 2.0*pi*c*spectral_width_nm/pump_wavelength**2
 real :: poling_period, low_wvln_idler, high_wvln_idler
 contains
real function temp_factor(temp)
real :: temp
real, parameter :: alpha=6.7e-6, beta=11.0e-9
 temp_factor=1.0+alpha*(temp-25.0)+beta*(temp-25)*(temp-25)
end function
real function crystal_length_temp()
 crystal_length_temp=crystal_length*temp_factor(temp)
end function
real function poling_period_temp()
 poling_period_temp=poling_period*temp_factor(temp)
end function

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
end function
 real function ktp_index_new(lambda,temp)
 implicit none
 real :: lambda,ll,n,n1,n2,temp
 real :: A=2.12725, B=1.18431, CC=5.14852e-2, D=0.6603
 real :: E= 100.00506, F=9.68956e-3
  ll  = lambda*1.0e6
  n1  = n_one(lambda)
  n2  = n_two(lambda)
  n   = sqrt(A + B/(1.0- CC/(ll**2)) + D/(1.0-E/(ll**2)) -F*ll**2)
  ktp_index_new = n + n1*(temp-25.0) + n2*(temp-25.0)**2.0
 end function
 real function n_one(lambda)
 real ::lambda,ll
  ll  = lambda*1.0e6
  n_one = 9.9587e-6 + 9.9228e-6/(ll) - 8.9603e-6/((ll)**2) + 4.1010e-6/((ll)**3)
 end function n_one
 real function n_two(lambda)
 implicit none
 real ::lambda,ll
  ll  = lambda*1.0e6
  n_two = -1.1882e-8 + 10.459e-8/(ll) - 9.8136e-8/((ll)**2) + 3.1481e-8/((ll)**3)
 end function n_two
end module vars_and_funcs


program pair_detection_simulator_convergence_test
use vars_and_funcs
use vsl_stream
implicit none
 real, parameter :: number_of_temps = 1
 real :: k_a_signal, k_b_signal, k_a_idler, k_b_idler, k_degen,                 &
         signal_gamma, idler_gamma, signal_aperture_angle, idler_aperture_angle
 real, dimension(nn+1) :: signal_k, signal_polar, signal_azimuth,              &
                        idler_k,  idler_polar,  idler_azimuth,                 &
                        signal_kx, signal_ky, signal_kz,                       &
                        idler_kx,  idler_ky,  idler_kz,                        &
                        wavefunction
 real, dimension(number_of_temps) :: k_hypervolume, average,                   &
                                     point_value, phase_mismatch
 integer(kind=4) :: i = 1, j = 1
 real    :: temp, idler_bw_wvln,                                   &
            signal_angle, signal_aperture, idler_angle, idler_aperture

 call read_parameters()
 call initialize_stream()

 average(:)=0.0

 do j = 1, repetitions

 !write(*,*) 'repetition ', j, ' of ',repetitions

 call set_variables()

 call determine_k_limits_and_hypervolume(k_a_signal, k_b_signal,               &
                                         k_a_idler,  k_b_idler,                &
                                         k_hypervolume(i), k_degen)

 call generate_random_photons(k_a_signal, k_b_signal, signal_aperture_angle,   &
                          signal_k, signal_polar, signal_azimuth, nn, k_degen)
 call generate_random_photons(k_a_idler, k_b_idler, idler_aperture_angle,      &
                          idler_k,  idler_polar,  idler_azimuth,  nn, k_degen)

 call rotate_to_aperture_angle(signal_gamma,signal_k,signal_polar,signal_azimuth,&
                  signal_kx,signal_ky,signal_kz,nn+1)
 call rotate_to_aperture_angle(idler_gamma, idler_k, idler_polar, idler_azimuth, &
                  idler_kx, idler_ky, idler_kz, nn+1)

 call evaluate_wavefunction(signal_k, signal_kx, signal_ky, signal_kz,         &
                            idler_k,  idler_kx,  idler_ky,  idler_kz,          &
                            wavefunction, temp, k_a_signal, k_b_signal,        &
                            k_a_idler, k_b_idler,                              &
                            point_value(i), phase_mismatch(i))!!!

 average(i) = average(i)+sum(wavefunction(1:nn))/nn; wavefunction(:)=0.

 end do

 average(i) = average(i)/repetitions

 call write_data_to_file( idler_aperture, idler_bw_wvln, average(i),           &
                          phase_mismatch(i), k_hypervolume(i), point_value(i) )
 call deinitialize_stream()

contains

 subroutine read_parameters()
 implicit none
 character(len=30) :: value
  if(command_argument_count().eq.0)then
   write(*,*)'working with defaults, gotta give me 3 or 4 arguments like this:'
   write(*,*)'./coherence signal_aperture[mm] idler_aperture[mm] external_signal_angle[degrees] idler_angle[mm] idler_bw[nm] temp[C]'
   write(*,*)'[defaults]./coherence 1.0 1.0 2.4 2.4 20 51.5'
   signal_aperture = 0.001
   idler_aperture  = 0.001
   signal_angle  = 2.4
   idler_angle   = 2.4
   idler_bw_wvln = 20.0e-9
   temp = work_temp
  elseif(command_argument_count().eq.6)then
   call get_command_argument(1,value)
   read(value,*) signal_aperture; signal_aperture=signal_aperture/1e3
   call get_command_argument(2,value)
   read(value,*) idler_aperture; idler_aperture=idler_aperture/1e3
   call get_command_argument(3,value)
   read(value,*) signal_angle;
   call get_command_argument(4,value)
   read(value,*) idler_angle
   call get_command_argument(5,value)
   read(value,*) idler_bw_wvln; idler_bw_wvln=idler_bw_wvln/1e9
   call get_command_argument(6,value)
   read(value,*) temp
  else
   write(*,*) 'wrong number of arguments, try no arguments to see some instructions'
   stop
  endif
  poling_period = pump_wavelength/(ktp_index(pump_wavelength,axial_pm_temp)    &
                                  -ktp_index(2.0*pump_wavelength,axial_pm_temp))
  poling_period = poling_period/temp_factor(axial_pm_temp)
write(*,*)'poling_period',poling_period
 end subroutine read_parameters

 subroutine set_variables()
 implicit none
  low_wvln_idler  = 2.0*pump_wavelength-0.5*idler_bw_wvln
  high_wvln_idler = 2.0*pump_wavelength+0.5*idler_bw_wvln
  signal_aperture_angle = &
   tan(0.5*signal_aperture*cos(signal_angle*pi/180.0)/crystal_dist)            &
    /ktp_index(pump_wavelength*2.0,temp)
  idler_aperture_angle  = &
   tan(0.5*idler_aperture*cos(idler_angle*pi/180.0)/crystal_dist)              &
    /ktp_index(pump_wavelength*2.0,temp)
  signal_gamma  = asin(sin(signal_angle*(pi/180.0))                            &
                      /ktp_index(pump_wavelength*2.0,temp))
  idler_gamma   =-asin(sin(idler_angle*(pi/180.0))                             &
                      /ktp_index(pump_wavelength*2.0,temp))
 end subroutine set_variables

 subroutine determine_k_limits_and_hypervolume(k_low_signal, k_high_signal,    &
                                               k_low_idler,  k_high_idler,     &
                                               hypervolume,  k_degen)
 implicit none
 real, intent(out) :: k_low_signal, k_high_signal,                             &
                      k_low_idler, k_high_idler, hypervolume, k_degen
  k_degen = 2.0*pi*ktp_index(2.0*pump_wavelength,temp)/(2.0*pump_wavelength)
  k_low_signal  = min( &
           ktp_index(high_wvln_signal,temp)*2.0*pi/high_wvln_signal,&
           ktp_index(low_wvln_signal,temp)*2.0*pi/low_wvln_signal)
  k_high_signal = max( &
           ktp_index(high_wvln_signal,temp)*2.0*pi/high_wvln_signal,&
           ktp_index(low_wvln_signal,temp)*2.0*pi/low_wvln_signal)
  k_low_idler   = min( &
           ktp_index(high_wvln_idler,temp)*2.0*pi/high_wvln_idler,&
           ktp_index(low_wvln_idler,temp)*2.0*pi/low_wvln_idler)
  k_high_idler  = max( &
           ktp_index(high_wvln_idler,temp)*2.0*pi/high_wvln_idler,&
           ktp_index(low_wvln_idler,temp)*2.0*pi/low_wvln_idler)
  hypervolume = (2.0*pi/3.0)**2*(k_b_signal**3-k_a_signal**3)                  &
              * (k_b_idler**3-k_a_idler**3)                                    &
              * (1.0-cos(signal_aperture_angle))*(1.0-cos(idler_aperture_angle))
 end subroutine determine_k_limits_and_hypervolume

 subroutine write_data_to_file(idler_aperture, idler_bw, average, phase_mismatch, k_hypervolume, point_value)
 implicit none
 real, intent(in)    :: idler_aperture, idler_bw, average, phase_mismatch, k_hypervolume, point_value
 character(80) fmt_list
 open(unit=10,file='coherence.dat')
 write(10,*) "#idler_aprtr[mm] idlr_bw[nm] 2nd_order_coherence phase_msmtch[pi] pnt_value"
 fmt_list = '(E12.4, E12.4, E12.4, E12.4, E12.4)'
 write(10,fmt_list) idler_aperture*1e3, idler_bw*1e9, average*k_hypervolume, phase_mismatch/pi, point_value
 close(10)
 end subroutine
end program pair_detection_simulator_convergence_test
