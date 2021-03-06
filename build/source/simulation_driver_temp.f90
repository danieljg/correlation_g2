module vars_and_funcs
implicit none
 integer, parameter :: nn = 1e3,                                               &
                       repetitions = 1
 real,    parameter ::                                                         &
   min_temp          = 40.0, max_temp   = 60.0, temp_step      = 0.5,          &
   axial_pm_temp     = 65.0,                                                   &
   crystal_dist      = 0.5,                                                    &
   pump_power        = 0.030,                                                  &
   pump_wavelength   = 406.118e-9,                                             &
   low_wvln_signal   = 805.0e-9,                                               &
   high_wvln_signal  = 820.0e-9,                                               &
   low_wvln_idler    = 805.0e-9,                                               &
   high_wvln_idler   = 820.0e-9,                                               &
   beam_waist        = 70.0e-6,                   &!!!!GOT IT RIGHT 70.0e-6
   spectral_width_nm = 0.042e-9,                  &!!!!GOT IT RIGHT 0.042e-9
   crystal_length    = 0.005,                                                  &
   d_eff             = 1.0e-12,                   &!!!!GET THESE NUMBERS RIGHT
   pi                = 4.0*atan(1.0),                                          &
   c                 = 3.0e8
 real, parameter :: omega_pump     = 2.0*pi*c/pump_wavelength,                 &
                    spectral_width = 2.0*pi*c*spectral_width_nm/pump_wavelength**2
 real :: poling_period
 contains
 real function omega(k_len,temp)
 real :: k_len, temp
 omega = 2.0*pi*c/lambda_of_k(k_len,temp)
 end function
 real function lambda_of_k(k_len,temp)
 real :: k_len, y_left, y_right, lambda_left, lambda_right
 real :: m, b, lambda_save, lambda_new, y_new, temp
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
 lambda_of_k=lambda_new
 end function lambda_of_k

real function temp_factor(temp)
real temp
real, parameter :: alpha=6.7e-6, beta=11.0e-9
 temp_factor=1.0+alpha*(temp-25.0)+beta*(temp-25)*(temp-25)
end function
real function crystal_length_temp(temp)
real temp
 crystal_length_temp=crystal_length*temp_factor(temp)
end function
real function poling_period_temp(temp)
real temp
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
end module vars_and_funcs


program pair_detection_simulator_temp_dependence
use vars_and_funcs
use vsl_stream
implicit none
 integer, parameter ::                                                         &
   number_of_temps = 1+ceiling(abs(max_temp-min_temp)/temp_step)
 real :: k_low_limit, k_high_limit,                                            &
         signal_gamma, idler_gamma, signal_aperture_angle, idler_aperture_angle
 real :: signal_aperture, idler_aperture,  k_degen
 real, dimension(nn+1) :: signal_k, signal_polar, signal_azimuth,              &
                        idler_k,  idler_polar,  idler_azimuth,                 &
                        signal_kx, signal_ky, signal_kz,                       &
                        idler_kx,  idler_ky,  idler_kz,                        &
                        wavefunction, erf_factor
 real, dimension(number_of_temps) :: k_hypervolume, average,                   &
                                     point_value, phase_mismatch

 integer(kind=4) :: i = 1, j = 1
 real    :: temp = min_temp, signal_angle, idler_angle

 call read_parameters()
 call initialize_stream()

 average(:)=0.0

 do i = 1, number_of_temps

 write(*,*) 'temp ', i, ' of ',number_of_temps

 do j = 1, repetitions

! write(*,*) 'repetition ', j, ' of ',repetitions

 call set_variables()

 call determine_k_limits_and_hypervolume(k_low_limit, k_high_limit,            &
                                         k_hypervolume(i), k_degen, temp)

 call generate_random_photons_uniform(k_low_limit, k_high_limit, signal_aperture_angle,&
                          signal_k, signal_polar, signal_azimuth, nn, k_degen)
 call generate_random_photons_gaussian(k_low_limit, k_high_limit, idler_aperture_angle,&
                          idler_k,  idler_polar,  idler_azimuth, erf_factor,        &
                          nn, k_degen, spectral_width, signal_k, temp)

 call rotate_to_aperture_angle(signal_gamma,signal_k,signal_polar,signal_azimuth,&
                  signal_kx,signal_ky,signal_kz,nn+1)
 call rotate_to_aperture_angle(idler_gamma, idler_k, idler_polar, idler_azimuth, &
                  idler_kx, idler_ky, idler_kz, nn+1)

 call evaluate_wavefunction(signal_k, signal_kx, signal_ky, signal_kz,         &
                            idler_k,  idler_kx,  idler_ky,  idler_kz,          &
                            wavefunction, temp, erf_factor,                    &
                            point_value(i), phase_mismatch(i))!!!
 average(i) = average(i)+sum(wavefunction(1:nn))/nn; wavefunction(:)=0.

 end do

 average(i) = average(i)/repetitions
 temp = min_temp+(i)*temp_step

 end do

 call write_data_to_file(min_temp, temp_step, number_of_temps,                 &
                         average , phase_mismatch, k_hypervolume, point_value)
 call deinitialize_stream()

contains

 subroutine read_parameters()
 implicit none
 character(len=30) :: value
  if(command_argument_count().eq.0)then
   write(*,*)'working with defaults, gotta give me 3 or 4 arguments like this:'
   write(*,*)'./coherence signal_aperture[mm] idler_aperture[mm] external_signal_angle[degrees] {idler_angle[mm]}'
   write(*,*)'[defaults]./coherence 1.0 1.0 2.4 2.4'
   signal_aperture = 0.001
   idler_aperture  = 0.001
   signal_angle = 2.4
   idler_angle  = 2.4
  elseif(command_argument_count().eq.3)then
   call get_command_argument(1,value)
   read(value,*) signal_aperture; signal_aperture=signal_aperture/1e3
   call get_command_argument(2,value)
   read(value,*) idler_aperture; idler_aperture=idler_aperture/1e3
   call get_command_argument(3,value)
   read(value,*) signal_angle;   idler_angle = signal_angle
  elseif(command_argument_count().eq.4)then
   call get_command_argument(1,value)
   read(value,*) signal_aperture; signal_aperture=signal_aperture/1e3
   call get_command_argument(2,value)
   read(value,*) idler_aperture; idler_aperture=idler_aperture/1e3
   call get_command_argument(3,value)
   read(value,*) signal_angle;
   call get_command_argument(4,value)
   read(value,*) idler_angle
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

 subroutine determine_k_limits_and_hypervolume(k_low, k_high, hypervolume,     &
                                               k_degen, temp )
 implicit none
 real, intent(out) :: k_low, k_high, hypervolume, k_degen
 real, intent(in)  :: temp
 real :: k_low_signal, k_high_signal, k_low_idler, k_high_idler, omega_new
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
  k_low  = max( k_low_signal,  k_low_idler  )
  k_high = min( k_high_signal, k_high_idler )
  if( (omega_pump-omega(k_low,temp)).gt.omega(k_high,temp) )  then
   omega_new = omega_pump-omega(k_high,temp)
   k_low     = ktp_index(2.0*pi*c/omega_new,temp)*omega_new/c
  else
   omega_new = omega_pump-omega(k_low,temp)
   k_high    = ktp_index(2.0*pi*c/omega_new,temp)*omega_new/c
  endif
  hypervolume = (2.0*pi/3.0)**2*(k_high**3-k_low**3)**2                        &
              * (1.0-cos(signal_aperture_angle))*(1.0-cos(idler_aperture_angle))
 end subroutine determine_k_limits_and_hypervolume

 subroutine write_data_to_file(min_temp, temp_step, number_of_temps,           &
                            average, phase_mismatch, k_hypervolume, point_value)
 implicit none
 integer, intent(in) :: number_of_temps
 real, intent(in)    :: min_temp, temp_step,                                   &
                   average(number_of_temps), phase_mismatch(number_of_temps),  &
                   k_hypervolume(number_of_temps), point_value(number_of_temps)
 integer i
 character(80) fmt_list
 open(unit=10,file='coherence.dat')
 write(10,'(120a)') "#temperature phase_mismatch[pi] 2nd_order_coherence point_value"
 fmt_list = '(E12.4, E12.4, E12.4, E12.4)'
 do i=1,number_of_temps
 write(10,fmt_list) min_temp+(i-1)*temp_step, phase_mismatch(i)/pi, &
                    average(i)*k_hypervolume(i), point_value(i)
 end do
 close(10)
 end subroutine
end program pair_detection_simulator_temp_dependence
