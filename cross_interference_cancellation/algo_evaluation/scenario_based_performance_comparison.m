close all;
clear all;
clc;


addpath('victim_radar_hardware_specifications\')
addpath('aggressor_radar1_hardware_specification\')
addpath('aggressor_radar2_hardware_specification\')
addpath('baseband_signal_model\')
addpath('adaptive_noise_canceller\')
addpath('peak_detection\')
addpath('performance_metrics\')
addpath('proximal_sub_gradient_descent\')
addpath('sparse_henkel_matrix_decomposition\')
addpath('plots\')


sinr_int     = 0;
corr_int     = 0;

sinr_grad    = 0;
corr_grad    = 0;
time_grad    = 0;

sinr_anc    = 0;
corr_anc    = 0;
time_anc    = 0;

sinr_sparkle    = 0;
corr_sparkle    = 0;
time_sparkle    = 0;

itr = 1;
itr2 = 0;

disp('Running...');
for j = 1:itr


    ADC       = adcParam;      % victim radar's ADC specifications
    txchirp   = txchirpParam;  % victim radar's transmitted signal specifications
    intchirp  = intchirpParam1; % aggressor radar's transmitted signal specifications
    lowpaas   = filterParam;   % victim radar's lowpaas filter specifications

    %     target_range_list = 25;
    %     target_velocity_list = 3;
    %
    %     RCS_target_list = 0.2;
    %
    %     txAntennaGain_victim = 0.9;
    %     rxAntennaGain_victim = 0.9;
    %     txPower_victim = 1;
    %
    %     aggressor_range_list = 2;
    %     aggressor_velocity_list = 1;
    %     txAntennaGain_aggressor = 0.9;
    %     txPower_aggressor = 1;
    %
    %     [baseband_sig_reflections,baseband_sig_interferences] = scenario_based_adc_output...
    %     (ADC,txchirp,intchirp,lowpaas,...
    %     target_range_list,target_velocity_list,RCS_target_list,txAntennaGain_victim,rxAntennaGain_victim,txPower_victim,...
    %     aggressor_range_list,aggressor_velocity_list,txAntennaGain_aggressor,txPower_aggressor);

    %% data generation
    % reading parameter files


    % target
    target_range = 25;   % range of object
    target_velocity = 3; % velocity of object (+ve for velocity away from victim radar, -ve for velocity towards victim )
    % signal amplitudes
    tx_sig_amp   =   1;  % victim's transmitted signal amplitude
    rx_sig_amp   =   1;  % victim's received signal amplitude

    snr_db = -10; % SNR in dB

%     int_sig_amp  =  200; % amplitude of baseband interference signal after ADC in victim radar
%            
%     lambda = 5;
%     initial_step_size = 4;

%      int_sig_amp  =  600; % amplitude of baseband interference signal after ADC in victim radar
%            
%     lambda = 3;
%     initial_step_size = 2.4;

%% 1.
     int_sig_amp  =  1000; % amplitude of baseband interference signal after ADC in victim radar
           
     lambda = 5;
     initial_step_size = 4;

%% 2.
%       int_sig_amp  =  100; % amplitude of baseband interference signal after ADC in victim radar
%            
%      lambda = 0.5;
%      initial_step_size = 0.4;



    maxIter = 100;
    lowpaas.cutoff = 9e6;

    sir_db = 20*log10(tx_sig_amp*rx_sig_amp*(1/int_sig_amp))

    % baseband signal generation due to a singal object
    reflection_baseband_sig = tx_sig_amp*rx_sig_amp*baseband_reflection(target_range,target_velocity,txchirp,ADC,lowpaas);

    % aggressor location
    aggressor_range = 0;    % aggressor radar range from the victim (assumed to be at same(near) location)
    aggressor_velocity = 0; % aggressor's velocity

    % baseband signal generation due a single aggressor's transmitted signal in victim radar
    interference_baseband_sig =  int_sig_amp*baseband_interference(aggressor_range,aggressor_velocity,intchirp,txchirp,ADC,lowpaas);

%     figure;
%     imagesc(abs(reshape(interference_baseband_sig,[ ADC.count_sample ADC.count_chirp])))
%     xlabel('Chirp Index','FontSize',18)
%     ylabel('Sample Index','FontSize',18)


    % interference quantity
    int_samples = double(abs(interference_baseband_sig) > 0);% number of non-zero samples in interference baseband signal
    int_percentage = sum(sum(int_samples))*(1/(ADC.count_sample*ADC.count_chirp))*100; % percentange of interference

    % beat signals
    beat_wo_int = reflection_baseband_sig;  % beat signal without interference
    beat_wi_int = reflection_baseband_sig + interference_baseband_sig; % beat signal with interference

    % Thermal noise generation
    %     snr_db = -60; % SNR in dB
    beat_power_wo_int  = (tx_sig_amp*rx_sig_amp)^2; %signal power
    noise_power_wo_int = beat_power_wo_int*(10^(-1*(snr_db/10))); % noise power
    cmplx_noise_wo_int = sqrt(noise_power_wo_int /2)*complex(randn(size(beat_wo_int )),randn(size(beat_wo_int ))); % noise signal

    beat_wo_int = beat_wo_int + cmplx_noise_wo_int;
   
    %     hold on
    %title('noise + interference')
    % beat_wo_int = beat_wo_int + cmplx_noise_wo_int; % beat signal without interference and with noise
    beat_wi_int = beat_wi_int + cmplx_noise_wo_int; % beat signal with interference and with noise

    %     beat_wi_int = beat_wi_int ;%%%%%%%%%%%

    % peak in Range-doppler spectrum grid
    grid_size = [ADC.count_sample, ADC.count_chirp]; % size of range-doppler grid
    peakIdx = peak_bin(target_range,target_velocity,grid_size,txchirp,ADC); % peak bin in Range-Doppler spectrum

    data_sim_wi_int             =   beat_wi_int; % data for processing
    data_sim_wi_int_mat         =   reshape(data_sim_wi_int,grid_size); % N(samples per chirp)xL(number of chirps)


%     freq_vs_time_plot(intchirp,txchirp,ADC,lowpaas,2,0,0)
    xlabel('time ')

    %     int = zeros(512,64);
    %     int(56:60,:) = 0.1*complex(randn(5,64),randn(5,64));
    %     data_sim_wi_int_mat =  data_sim_wi_int_mat +int; %%%%%%%%%

    RFFT_data_sim_wi_int        =   fft(data_sim_wi_int_mat,ADC.count_sample,1); % range FFT
    RDFFT_data_sim_wi_int       =   fft( RFFT_data_sim_wi_int,ADC.count_chirp,2); % Range-Doppler FFT

    RD_mat         = RDFFT_data_sim_wi_int(1:ADC.count_sample/2,:) ; % positive half of Range-Doppler spectrum
    %     figure;
    %     surf(10*log10(abs(RD_mat)),EdgeColor="none")

    n_sig_bin      = [2 ,2];  % number of signal bins at left and right of peak bin for range and doppler respectively
    n_noise_bin    = [6 ,6];  % number of noise bins at left and right of peak bin for range and doppler respectively
    % local SINR calculation at peak bin
    [loc_sinr_int,loc_sinr_range_int,loc_sinr_doppler_int ] = local_snr(RD_mat,peakIdx(1,:),n_sig_bin, n_noise_bin);
    sinr_int = sinr_int + loc_sinr_int;

    % correlation coefficient calculation
    corr_coef_int = correlation_coeficient(beat_wo_int,beat_wi_int);
    corr_int = corr_int + corr_coef_int;


%     lambda = 10;
%     mu = 0.1;
%     addpath('SALSA\')
%     
%      [a_r,a_i] = SALSA( data_sim_wi_int_mat(:,5), lambda, mu,maxIter);
% 
%      figure;
%      plot(real(a_i))


    %% Interference Mitigation Alorithm
    %Sub-Gradient Descent
    %disp('Gradient Descent...')
    % Hyper parameters
    %     lambda = 0.4;
    %     initial_step_size = 0.5;
    limit = 1;
         maxIter = 100;
    %
    tic;
    [data_corr_mat_grad] = prox_sub_grad_des(data_sim_wi_int_mat,lambda,initial_step_size,limit,maxIter); %Interference matrix
    exec_time_grad_des = toc;

%      figure;
%     plot(real(cmplx_noise_wo_int(1:512) + real(interference_baseband_sig(1:512))),'b')
%     hold on
%     plot(real(data_corr_mat_grad(:,1)),'r')

    %     plot(real(data_corr_mat_grad(:,1)) ,'r');
    %     title('estimated intereference')
    %     figure;
    %     plot(real(interference_baseband_sig(1:512))-real(data_corr_mat_grad(:,1)))
    %     title('interference - estimated interference')
    %     hold on
    %      plot(real(int(:,1)),'b');
    %     figure;
    %     plot(real(data_sim_wi_int_mat(:,1))-real(data_corr_mat_grad(:,1)));
    %     figure;
    %     plot(real(data_sim_wi_int_mat(:,1)));
    %fprintf('*** data set 1 gradient descent time : %d ***\n',exec_time_grad_des);
    data_miti_mat_grad       =   data_sim_wi_int_mat - data_corr_mat_grad; % interefernce cancelled matrix
    RFFTGrad_data_miti       =   fft(data_miti_mat_grad,ADC.count_sample,1); % Range fft
    RDFFTGrad_data_miti      =   fft(RFFTGrad_data_miti,ADC.count_chirp,2); % Range-Doopler fft

    RD_mat         = RDFFTGrad_data_miti(1:ADC.count_sample/2,:); % positive half of RD FFT
    [loc_sinr_grad,loc_sinr_range_grad,loc_sinr_doppler_grad ] = local_snr(RD_mat,peakIdx(1,:),n_sig_bin, n_noise_bin); % SINR computation after interference cancellation
    sinr_grad = sinr_grad + loc_sinr_grad;

    % correlation coefficient
    corr_coef_grad = correlation_coeficient(beat_wo_int,reshape(data_miti_mat_grad,[numel(data_miti_mat_grad) 1]));
    corr_grad = corr_grad + corr_coef_grad;

    % execution time
    time_grad = time_grad + exec_time_grad_des;





    %Adaptive Noise canceller
    addpath('./adaptive_noise_canceller/')
    %disp('Adaptive Noise Canceller...')
    flen = 8; % lenght of adaptive filter 50 for 256
    tic;
    % anc_RDFFT_miti = anc_RDFFT(RFFT_data,NumSamplePerChirp,BlockSize,flen);
    eout = anc_RDFFT(RFFT_data_sim_wi_int,ADC.count_sample,ADC.count_chirp,flen);
    exec_time_anc  = toc;

    RDAfterANC_postive      = fft(eout,ADC.count_chirp, 2); %fftshift(fft(eout, BlockSize, 2), 2);
    %fprintf('*** RDFFT data set 1 adaptive noise cenceller time : %d ***\n',exec_time_anc);

    RD_mat         = RDAfterANC_postive;
    [loc_sinr_anc,loc_sinr_range_anc,loc_sinr_doppler_anc ] = local_snr(RD_mat,peakIdx(1,:),n_sig_bin, n_noise_bin);
    sinr_anc = sinr_anc + loc_sinr_anc;

    RDAfterANC_full = [RDAfterANC_postive;RDFFT_data_sim_wi_int((ADC.count_sample /2)+1:end,:)];
    data_miti_mat_anc = ifft(ifft(RDAfterANC_full,ADC.count_chirp,2),ADC.count_sample ,1);
    corr_coef_anc = correlation_coeficient(beat_wo_int,reshape(data_miti_mat_anc,[numel(data_miti_mat_anc) 1]));
    corr_anc      = corr_anc + corr_coef_anc;

    time_anc = time_anc + exec_time_anc;


  




    if (j <= itr2)
        %Sparse and Henkel Matrix Decomposition
        addpath('./sparse_henkel_matrix_decomposition/lowRaS_MC/')
        disp('IM-SPARKLE...')
        beta_1 = 0.1;           % 0.5
        mu     = 0.02;          % 0.05
        tau    = 0.02;          % 0.02
        k_beta = 1.6;
        k_mu   = 1.2;
        R      = 10;
        tic;
        [data_miti_vec_sparkle, i_lowRaS, rerr] = lowRaS_Hankel(data_sim_wi_int, R, beta_1, mu, tau, k_beta, k_mu);

        exec_time_sparkle = toc;

        % fprintf('*** data set 1 sparse and hankel matrix decomposition time : %d ***\n',exec_time_sparkle);
        data_miti_mat_sparkle       = reshape(data_miti_vec_sparkle,[ADC.count_sample , ADC.count_chirp]);
        RangeFFTsparkle_data_miti   = fft(data_miti_mat_sparkle,ADC.count_sample,1);
        RDFFTsparkle_data_miti      = fft(RangeFFTsparkle_data_miti,ADC.count_chirp,2);

        RD_mat         = RDFFTsparkle_data_miti(1:ADC.count_sample/2,:) ;
        [loc_sinr_sparkle,loc_sinr_range_sparkle,loc_sinr_doppler_sparkle ] = local_snr(RD_mat,peakIdx(1,:),n_sig_bin, n_noise_bin);
        sinr_sparkle = sinr_sparkle + loc_sinr_sparkle;

        corr_coef_sparkle = correlation_coeficient(beat_wo_int,data_miti_vec_sparkle);
        corr_sparkle      = corr_sparkle + corr_coef_sparkle;

        time_sparkle = time_sparkle + exec_time_sparkle;

        %             figure;
        %             surf((20*log10(abs(RDFFTsparkle_data_miti))),'EdgeColor','none');%(1:2:end,1:2:end)
        %             hold on
        %             colormap default
    end






end

sinr_grad    = sinr_grad/itr;
sinr_int     = sinr_int/itr;
corr_grad    = corr_grad/itr ;
corr_int     = corr_int/itr;
time_grad    = time_grad/itr;

sinr_anc    = sinr_anc/itr;
sinr_sparkle     = sinr_sparkle/itr;
corr_anc    = corr_anc/itr ;
corr_sparkle     = corr_sparkle/itr;
time_anc    = time_anc/itr;
time_sparkle    = time_sparkle/itr;

disp('*********************************************');
fprintf('Interference Properties\n');
fprintf('Samples percentage: %s \n', num2str(int_percentage));
fprintf('SINR : %s \n', num2str(sinr_int));
fprintf('Correlation Coef. : %s \n', num2str(corr_int));
fprintf('Correlation Coef. in magn and angle : %s %s\n', num2str(abs(corr_int)),num2str(angle(corr_int)));

disp('*********************************************');
fprintf('Sub-gradient Descent\n');
fprintf('SINR : %s \n', num2str(sinr_grad));
fprintf('Correlation Coef. : %s\n', num2str(corr_grad));
fprintf('Correlation Coef. in magn and angle : %s %s\n', num2str(abs(corr_grad)),num2str(angle(corr_grad)));
fprintf('Execution Time : %s\n', num2str(time_grad));

disp('*********************************************');
fprintf('ANC\n');
fprintf('SINR : %s \n', num2str(sinr_anc));
fprintf('Correlation Coef. : %s\n', num2str(corr_anc));
fprintf('Correlation Coef. in magn and angle : %s %s\n', num2str(abs(corr_anc)),num2str(angle(corr_anc)));
fprintf('Execution Time : %s\n', num2str(time_anc));

disp('*********************************************');
fprintf('SPARKLE\n');
fprintf('SINR : %s \n', num2str(sinr_sparkle));
fprintf('Correlation Coef. : %s\n', num2str(corr_sparkle));
fprintf('Correlation Coef. in magn and angle : %s %s\n', num2str(abs(corr_sparkle)),num2str(angle(corr_sparkle)));
fprintf('Execution Time : %s\n', num2str(time_sparkle));


clear all;

sinr_int     = 0;
corr_int     = 0;

sinr_grad    = 0;
corr_grad    = 0;
time_grad    = 0;

sinr_anc    = 0;
corr_anc    = 0;
time_anc    = 0;

sinr_sparkle    = 0;
corr_sparkle    = 0;
time_sparkle    = 0;

itr = 1;
itr2 = 0;


for j = 1:itr


    ADC       = adcParam;      % victim radar's ADC specifications
    txchirp   = txchirpParam;  % victim radar's transmitted signal specifications
    intchirp  = intchirpParam1; % aggressor radar's transmitted signal specifications
    lowpaas   = filterParam;   % victim radar's lowpaas filter specifications

    %     target_range_list = 25;
    %     target_velocity_list = 3;
    %
    %     RCS_target_list = 0.2;
    %
    %     txAntennaGain_victim = 0.9;
    %     rxAntennaGain_victim = 0.9;
    %     txPower_victim = 1;
    %
    %     aggressor_range_list = 2;
    %     aggressor_velocity_list = 1;
    %     txAntennaGain_aggressor = 0.9;
    %     txPower_aggressor = 1;
    %
    %     [baseband_sig_reflections,baseband_sig_interferences] = scenario_based_adc_output...
    %     (ADC,txchirp,intchirp,lowpaas,...
    %     target_range_list,target_velocity_list,RCS_target_list,txAntennaGain_victim,rxAntennaGain_victim,txPower_victim,...
    %     aggressor_range_list,aggressor_velocity_list,txAntennaGain_aggressor,txPower_aggressor);

    %% data generation
    % reading parameter files


    % target
    target_range = 25;   % range of object
    target_velocity = 3; % velocity of object (+ve for velocity away from victim radar, -ve for velocity towards victim )
    % signal amplitudes
    tx_sig_amp   =   1;  % victim's transmitted signal amplitude
    rx_sig_amp   =   1;  % victim's received signal amplitude

    snr_db = -10; % SNR in dB

%     int_sig_amp  =  200; % amplitude of baseband interference signal after ADC in victim radar
%            
%     lambda = 5;
%     initial_step_size = 4;

%      int_sig_amp  =  600; % amplitude of baseband interference signal after ADC in victim radar
%            
%     lambda = 3;
%     initial_step_size = 2.4;

% %% 1.
%      int_sig_amp  =  1000; % amplitude of baseband interference signal after ADC in victim radar
%            
%      lambda = 5;
%      initial_step_size = 4;

%% 2.
      int_sig_amp  =  500; % amplitude of baseband interference signal after ADC in victim radar
           
     lambda = 2.5;
     initial_step_size = 2;


    maxIter = 100;
    lowpaas.cutoff = 9e6;

    sir_db = 20*log10(tx_sig_amp*rx_sig_amp*(1/int_sig_amp))

    % baseband signal generation due to a singal object
    reflection_baseband_sig = tx_sig_amp*rx_sig_amp*baseband_reflection(target_range,target_velocity,txchirp,ADC,lowpaas);

    % aggressor location
    aggressor_range = 0;    % aggressor radar range from the victim (assumed to be at same(near) location)
    aggressor_velocity = 0; % aggressor's velocity

    % baseband signal generation due a single aggressor's transmitted signal in victim radar
    interference_baseband_sig =  int_sig_amp*baseband_interference(aggressor_range,aggressor_velocity,intchirp,txchirp,ADC,lowpaas);
% 
%     figure;
%     imagesc(abs(reshape(interference_baseband_sig,[ ADC.count_sample ADC.count_chirp])))
%     


    % interference quantity
    int_samples = double(abs(interference_baseband_sig) > 0);% number of non-zero samples in interference baseband signal
    int_percentage = sum(sum(int_samples))*(1/(ADC.count_sample*ADC.count_chirp))*100; % percentange of interference

    % beat signals
    beat_wo_int = reflection_baseband_sig;  % beat signal without interference
    beat_wi_int = reflection_baseband_sig + interference_baseband_sig; % beat signal with interference

    % Thermal noise generation
    %     snr_db = -60; % SNR in dB
    beat_power_wo_int  = (tx_sig_amp*rx_sig_amp)^2; %signal power
    noise_power_wo_int = beat_power_wo_int*(10^(-1*(snr_db/10))); % noise power
    cmplx_noise_wo_int = sqrt(noise_power_wo_int /2)*complex(randn(size(beat_wo_int )),randn(size(beat_wo_int ))); % noise signal

    beat_wo_int = beat_wo_int + cmplx_noise_wo_int;
   
    %     hold on
    %title('noise + interference')
    % beat_wo_int = beat_wo_int + cmplx_noise_wo_int; % beat signal without interference and with noise
    beat_wi_int = beat_wi_int + cmplx_noise_wo_int; % beat signal with interference and with noise

    %     beat_wi_int = beat_wi_int ;%%%%%%%%%%%

    % peak in Range-doppler spectrum grid
    grid_size = [ADC.count_sample, ADC.count_chirp]; % size of range-doppler grid
    peakIdx = peak_bin(target_range,target_velocity,grid_size,txchirp,ADC); % peak bin in Range-Doppler spectrum

    data_sim_wi_int             =   beat_wi_int; % data for processing
    data_sim_wi_int_mat         =   reshape(data_sim_wi_int,grid_size); % N(samples per chirp)xL(number of chirps)

%      freq_vs_time_plot(intchirp,txchirp,ADC,lowpaas,2,0,0)

    %     int = zeros(512,64);
    %     int(56:60,:) = 0.1*complex(randn(5,64),randn(5,64));
    %     data_sim_wi_int_mat =  data_sim_wi_int_mat +int; %%%%%%%%%

    RFFT_data_sim_wi_int        =   fft(data_sim_wi_int_mat,ADC.count_sample,1); % range FFT
    RDFFT_data_sim_wi_int       =   fft( RFFT_data_sim_wi_int,ADC.count_chirp,2); % Range-Doppler FFT

    RD_mat         = RDFFT_data_sim_wi_int(1:ADC.count_sample/2,:) ; % positive half of Range-Doppler spectrum
    %     figure;
    %     surf(10*log10(abs(RD_mat)),EdgeColor="none")

    n_sig_bin      = [2 ,2];  % number of signal bins at left and right of peak bin for range and doppler respectively
    n_noise_bin    = [6 ,6];  % number of noise bins at left and right of peak bin for range and doppler respectively
    % local SINR calculation at peak bin
    [loc_sinr_int,loc_sinr_range_int,loc_sinr_doppler_int ] = local_snr(RD_mat,peakIdx(1,:),n_sig_bin, n_noise_bin);
    sinr_int = sinr_int + loc_sinr_int;

    % correlation coefficient calculation
    corr_coef_int = correlation_coeficient(beat_wo_int,beat_wi_int);
    corr_int = corr_int + corr_coef_int;


    %% Interference Mitigation Alorithm
    %Sub-Gradient Descent
    %disp('Gradient Descent...')
    % Hyper parameters
    %     lambda = 0.4;
    %     initial_step_size = 0.5;
    limit = 0.01;
    %     maxIter = 100;
    %
    tic;
    [data_corr_mat_grad] = prox_sub_grad_des(data_sim_wi_int_mat,lambda,initial_step_size,limit,maxIter); %Interference matrix
    exec_time_grad_des = toc;

%      figure;
%     plot(real(cmplx_noise_wo_int(1:512) + real(interference_baseband_sig(1:512))),'b')
%     hold on
%     plot(real(data_corr_mat_grad(:,1)),'r')

    %     plot(real(data_corr_mat_grad(:,1)) ,'r');
    %     title('estimated intereference')
    %     figure;
    %     plot(real(interference_baseband_sig(1:512))-real(data_corr_mat_grad(:,1)))
    %     title('interference - estimated interference')
    %     hold on
    %      plot(real(int(:,1)),'b');
    %     figure;
    %     plot(real(data_sim_wi_int_mat(:,1))-real(data_corr_mat_grad(:,1)));
    %     figure;
    %     plot(real(data_sim_wi_int_mat(:,1)));
    %fprintf('*** data set 1 gradient descent time : %d ***\n',exec_time_grad_des);
    data_miti_mat_grad       =   data_sim_wi_int_mat - data_corr_mat_grad; % interefernce cancelled matrix
    RFFTGrad_data_miti       =   fft(data_miti_mat_grad,ADC.count_sample,1); % Range fft
    RDFFTGrad_data_miti      =   fft(RFFTGrad_data_miti,ADC.count_chirp,2); % Range-Doopler fft

    RD_mat         = RDFFTGrad_data_miti(1:ADC.count_sample/2,:); % positive half of RD FFT
    [loc_sinr_grad,loc_sinr_range_grad,loc_sinr_doppler_grad ] = local_snr(RD_mat,peakIdx(1,:),n_sig_bin, n_noise_bin); % SINR computation after interference cancellation
    sinr_grad = sinr_grad + loc_sinr_grad;

    % correlation coefficient
    corr_coef_grad = correlation_coeficient(beat_wo_int,reshape(data_miti_mat_grad,[numel(data_miti_mat_grad) 1]));
    corr_grad = corr_grad + corr_coef_grad;

    % execution time
    time_grad = time_grad + exec_time_grad_des;



    %Adaptive Noise canceller
    addpath('./adaptive_noise_canceller/')
    %disp('Adaptive Noise Canceller...')
    flen = 8; % lenght of adaptive filter 50 for 256
    tic;
    % anc_RDFFT_miti = anc_RDFFT(RFFT_data,NumSamplePerChirp,BlockSize,flen);
    eout = anc_RDFFT(RFFT_data_sim_wi_int,ADC.count_sample,ADC.count_chirp,flen);
    exec_time_anc  = toc;

    RDAfterANC_postive      = fft(eout,ADC.count_chirp, 2); %fftshift(fft(eout, BlockSize, 2), 2);
    %fprintf('*** RDFFT data set 1 adaptive noise cenceller time : %d ***\n',exec_time_anc);

    RD_mat         = RDAfterANC_postive;
    [loc_sinr_anc,loc_sinr_range_anc,loc_sinr_doppler_anc ] = local_snr(RD_mat,peakIdx(1,:),n_sig_bin, n_noise_bin);
    sinr_anc = sinr_anc + loc_sinr_anc;

    RDAfterANC_full = [RDAfterANC_postive;RDFFT_data_sim_wi_int((ADC.count_sample /2)+1:end,:)];
    data_miti_mat_anc = ifft(ifft(RDAfterANC_full,ADC.count_chirp,2),ADC.count_sample ,1);
    corr_coef_anc = correlation_coeficient(beat_wo_int,reshape(data_miti_mat_anc,[numel(data_miti_mat_anc) 1]));
    corr_anc      = corr_anc + corr_coef_anc;

    time_anc = time_anc + exec_time_anc;

    if (j <= itr2)
        %Sparse and Henkel Matrix Decomposition
        addpath('./sparse_henkel_matrix_decomposition/lowRaS_MC/')
        disp('IM-SPARKLE...')
        beta_1 = 0.1;           % 0.5
        mu     = 0.02;          % 0.05
        tau    = 0.02;          % 0.02
        k_beta = 1.6;
        k_mu   = 1.2;
        R      = 10;
        tic;
        [data_miti_vec_sparkle, i_lowRaS, rerr] = lowRaS_Hankel(data_sim_wi_int, R, beta_1, mu, tau, k_beta, k_mu);

        exec_time_sparkle = toc;

        % fprintf('*** data set 1 sparse and hankel matrix decomposition time : %d ***\n',exec_time_sparkle);
        data_miti_mat_sparkle       = reshape(data_miti_vec_sparkle,[ADC.count_sample , ADC.count_chirp]);
        RangeFFTsparkle_data_miti   = fft(data_miti_mat_sparkle,ADC.count_sample,1);
        RDFFTsparkle_data_miti      = fft(RangeFFTsparkle_data_miti,ADC.count_chirp,2);

        RD_mat         = RDFFTsparkle_data_miti(1:ADC.count_sample/2,:) ;
        [loc_sinr_sparkle,loc_sinr_range_sparkle,loc_sinr_doppler_sparkle ] = local_snr(RD_mat,peakIdx(1,:),n_sig_bin, n_noise_bin);
        sinr_sparkle = sinr_sparkle + loc_sinr_sparkle;

        corr_coef_sparkle = correlation_coeficient(beat_wo_int,data_miti_vec_sparkle);
        corr_sparkle      = corr_sparkle + corr_coef_sparkle;

        time_sparkle = time_sparkle + exec_time_sparkle;

        %             figure;
        %             surf((20*log10(abs(RDFFTsparkle_data_miti))),'EdgeColor','none');%(1:2:end,1:2:end)
        %             hold on
        %             colormap default
    end






end

sinr_grad    = sinr_grad/itr;
sinr_int     = sinr_int/itr;
corr_grad    = corr_grad/itr ;
corr_int     = corr_int/itr;
time_grad    = time_grad/itr;

sinr_anc    = sinr_anc/itr;
sinr_sparkle     = sinr_sparkle/itr;
corr_anc    = corr_anc/itr ;
corr_sparkle     = corr_sparkle/itr;
time_anc    = time_anc/itr;
time_sparkle    = time_sparkle/itr;

disp('*********************************************');
fprintf('Interference Properties\n');
fprintf('Samples percentage: %s \n', num2str(int_percentage));
fprintf('SINR : %s \n', num2str(sinr_int));
fprintf('Correlation Coef. : %s \n', num2str(corr_int));
fprintf('Correlation Coef. in magn and angle : %s %s\n', num2str(abs(corr_int)),num2str(angle(corr_int)));

disp('*********************************************');
fprintf('Sub-gradient Descent\n');
fprintf('SINR : %s \n', num2str(sinr_grad));
fprintf('Correlation Coef. : %s\n', num2str(corr_grad));
fprintf('Correlation Coef. in magn and angle : %s %s\n', num2str(abs(corr_grad)),num2str(angle(corr_grad)));
fprintf('Execution Time : %s\n', num2str(time_grad));

disp('*********************************************');
fprintf('ANC\n');
fprintf('SINR : %s \n', num2str(sinr_anc));
fprintf('Correlation Coef. : %s\n', num2str(corr_anc));
fprintf('Correlation Coef. in magn and angle : %s %s\n', num2str(abs(corr_anc)),num2str(angle(corr_anc)));
fprintf('Execution Time : %s\n', num2str(time_anc));

disp('*********************************************');
fprintf('SPARKLE\n');
fprintf('SINR : %s \n', num2str(sinr_sparkle));
fprintf('Correlation Coef. : %s\n', num2str(corr_sparkle));
fprintf('Correlation Coef. in magn and angle : %s %s\n', num2str(abs(corr_sparkle)),num2str(angle(corr_sparkle)));
fprintf('Execution Time : %s\n', num2str(time_sparkle));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

sinr_int     = 0;
corr_int     = 0;

sinr_grad    = 0;
corr_grad    = 0;
time_grad    = 0;

sinr_anc    = 0;
corr_anc    = 0;
time_anc    = 0;

sinr_sparkle    = 0;
corr_sparkle    = 0;
time_sparkle    = 0;

itr = 1;
itr2 = 0;

disp('Running...');
for j = 1:itr


    ADC       = adcParam;      % victim radar's ADC specifications
    txchirp   = txchirpParam;  % victim radar's transmitted signal specifications
    intchirp  = intchirpParam2; % aggressor radar's transmitted signal specifications
    lowpaas   = filterParam;   % victim radar's lowpaas filter specifications

    %     target_range_list = 25;
    %     target_velocity_list = 3;
    %
    %     RCS_target_list = 0.2;
    %
    %     txAntennaGain_victim = 0.9;
    %     rxAntennaGain_victim = 0.9;
    %     txPower_victim = 1;
    %
    %     aggressor_range_list = 2;
    %     aggressor_velocity_list = 1;
    %     txAntennaGain_aggressor = 0.9;
    %     txPower_aggressor = 1;
    %
    %     [baseband_sig_reflections,baseband_sig_interferences] = scenario_based_adc_output...
    %     (ADC,txchirp,intchirp,lowpaas,...
    %     target_range_list,target_velocity_list,RCS_target_list,txAntennaGain_victim,rxAntennaGain_victim,txPower_victim,...
    %     aggressor_range_list,aggressor_velocity_list,txAntennaGain_aggressor,txPower_aggressor);

    %% data generation
    % reading parameter files


    % target
    target_range = 25;   % range of object
    target_velocity = 3; % velocity of object (+ve for velocity away from victim radar, -ve for velocity towards victim )
    % signal amplitudes
    tx_sig_amp   =   1;  % victim's transmitted signal amplitude
    rx_sig_amp   =   1;  % victim's received signal amplitude

    snr_db = -10; % SNR in dB

%     int_sig_amp  =  200; % amplitude of baseband interference signal after ADC in victim radar
%            
%     lambda = 5;
%     initial_step_size = 4;

%      int_sig_amp  =  600; % amplitude of baseband interference signal after ADC in victim radar
%            
%     lambda = 3;
%     initial_step_size = 2.4;

%% 1.
     int_sig_amp  =  1000; % amplitude of baseband interference signal after ADC in victim radar
           
     lambda = 5;
     initial_step_size = 4;

%% 2.
%       int_sig_amp  =  100; % amplitude of baseband interference signal after ADC in victim radar
%            
%      lambda = 0.5;
%      initial_step_size = 0.4;



    maxIter = 100;
    lowpaas.cutoff = 9e6;

    sir_db = 20*log10(tx_sig_amp*rx_sig_amp*(1/int_sig_amp))

    % baseband signal generation due to a singal object
    reflection_baseband_sig = tx_sig_amp*rx_sig_amp*baseband_reflection(target_range,target_velocity,txchirp,ADC,lowpaas);

    % aggressor location
    aggressor_range = 0;    % aggressor radar range from the victim (assumed to be at same(near) location)
    aggressor_velocity = 0; % aggressor's velocity

    % baseband signal generation due a single aggressor's transmitted signal in victim radar
    interference_baseband_sig =  int_sig_amp*baseband_interference(aggressor_range,aggressor_velocity,intchirp,txchirp,ADC,lowpaas);

%     figure;
%     imagesc(abs(reshape(interference_baseband_sig,[ ADC.count_sample ADC.count_chirp])))
%     xlabel('Chirp Index','FontSize',18)
%     ylabel('Sample Index','FontSize',18)


    % interference quantity
    int_samples = double(abs(interference_baseband_sig) > 0);% number of non-zero samples in interference baseband signal
    int_percentage = sum(sum(int_samples))*(1/(ADC.count_sample*ADC.count_chirp))*100; % percentange of interference

    % beat signals
    beat_wo_int = reflection_baseband_sig;  % beat signal without interference
    beat_wi_int = reflection_baseband_sig + interference_baseband_sig; % beat signal with interference

    % Thermal noise generation
    %     snr_db = -60; % SNR in dB
    beat_power_wo_int  = (tx_sig_amp*rx_sig_amp)^2; %signal power
    noise_power_wo_int = beat_power_wo_int*(10^(-1*(snr_db/10))); % noise power
    cmplx_noise_wo_int = sqrt(noise_power_wo_int /2)*complex(randn(size(beat_wo_int )),randn(size(beat_wo_int ))); % noise signal

    beat_wo_int = beat_wo_int + cmplx_noise_wo_int;
   
    %     hold on
    %title('noise + interference')
    % beat_wo_int = beat_wo_int + cmplx_noise_wo_int; % beat signal without interference and with noise
    beat_wi_int = beat_wi_int + cmplx_noise_wo_int; % beat signal with interference and with noise

    %     beat_wi_int = beat_wi_int ;%%%%%%%%%%%

    % peak in Range-doppler spectrum grid
    grid_size = [ADC.count_sample, ADC.count_chirp]; % size of range-doppler grid
    peakIdx = peak_bin(target_range,target_velocity,grid_size,txchirp,ADC); % peak bin in Range-Doppler spectrum

    data_sim_wi_int             =   beat_wi_int; % data for processing
    data_sim_wi_int_mat         =   reshape(data_sim_wi_int,grid_size); % N(samples per chirp)xL(number of chirps)


%       freq_vs_time_plot(intchirp,txchirp,ADC,lowpaas,2,0,0)

    %     int = zeros(512,64);
    %     int(56:60,:) = 0.1*complex(randn(5,64),randn(5,64));
    %     data_sim_wi_int_mat =  data_sim_wi_int_mat +int; %%%%%%%%%

    RFFT_data_sim_wi_int        =   fft(data_sim_wi_int_mat,ADC.count_sample,1); % range FFT
    RDFFT_data_sim_wi_int       =   fft( RFFT_data_sim_wi_int,ADC.count_chirp,2); % Range-Doppler FFT

    RD_mat         = RDFFT_data_sim_wi_int(1:ADC.count_sample/2,:) ; % positive half of Range-Doppler spectrum
    %     figure;
    %     surf(10*log10(abs(RD_mat)),EdgeColor="none")

    n_sig_bin      = [2 ,2];  % number of signal bins at left and right of peak bin for range and doppler respectively
    n_noise_bin    = [6 ,6];  % number of noise bins at left and right of peak bin for range and doppler respectively
    % local SINR calculation at peak bin
    [loc_sinr_int,loc_sinr_range_int,loc_sinr_doppler_int ] = local_snr(RD_mat,peakIdx(1,:),n_sig_bin, n_noise_bin);
    sinr_int = sinr_int + loc_sinr_int;

    % correlation coefficient calculation
    corr_coef_int = correlation_coeficient(beat_wo_int,beat_wi_int);
    corr_int = corr_int + corr_coef_int;


    %% Interference Mitigation Alorithm
    %Sub-Gradient Descent
    %disp('Gradient Descent...')
    % Hyper parameters
    %     lambda = 0.4;
    %     initial_step_size = 0.5;
    limit = 0.01;
    %     maxIter = 100;
    %
    tic;
    [data_corr_mat_grad] = prox_sub_grad_des(data_sim_wi_int_mat,lambda,initial_step_size,limit,maxIter); %Interference matrix
    exec_time_grad_des = toc;

%      figure;
%     plot(real(cmplx_noise_wo_int(1:512) + real(interference_baseband_sig(1:512))),'b')
%     hold on
%     plot(real(data_corr_mat_grad(:,1)),'r')

    %     plot(real(data_corr_mat_grad(:,1)) ,'r');
    %     title('estimated intereference')
    %     figure;
    %     plot(real(interference_baseband_sig(1:512))-real(data_corr_mat_grad(:,1)))
    %     title('interference - estimated interference')
    %     hold on
    %      plot(real(int(:,1)),'b');
    %     figure;
    %     plot(real(data_sim_wi_int_mat(:,1))-real(data_corr_mat_grad(:,1)));
    %     figure;
    %     plot(real(data_sim_wi_int_mat(:,1)));
    %fprintf('*** data set 1 gradient descent time : %d ***\n',exec_time_grad_des);
    data_miti_mat_grad       =   data_sim_wi_int_mat - data_corr_mat_grad; % interefernce cancelled matrix
    RFFTGrad_data_miti       =   fft(data_miti_mat_grad,ADC.count_sample,1); % Range fft
    RDFFTGrad_data_miti      =   fft(RFFTGrad_data_miti,ADC.count_chirp,2); % Range-Doopler fft

    RD_mat         = RDFFTGrad_data_miti(1:ADC.count_sample/2,:); % positive half of RD FFT
    [loc_sinr_grad,loc_sinr_range_grad,loc_sinr_doppler_grad ] = local_snr(RD_mat,peakIdx(1,:),n_sig_bin, n_noise_bin); % SINR computation after interference cancellation
    sinr_grad = sinr_grad + loc_sinr_grad;

    % correlation coefficient
    corr_coef_grad = correlation_coeficient(beat_wo_int,reshape(data_miti_mat_grad,[numel(data_miti_mat_grad) 1]));
    corr_grad = corr_grad + corr_coef_grad;

    % execution time
    time_grad = time_grad + exec_time_grad_des;



    %Adaptive Noise canceller
    addpath('./adaptive_noise_canceller/')
    %disp('Adaptive Noise Canceller...')
    flen = 8; % lenght of adaptive filter 50 for 256
    tic;
    % anc_RDFFT_miti = anc_RDFFT(RFFT_data,NumSamplePerChirp,BlockSize,flen);
    eout = anc_RDFFT(RFFT_data_sim_wi_int,ADC.count_sample,ADC.count_chirp,flen);
    exec_time_anc  = toc;

    RDAfterANC_postive      = fft(eout,ADC.count_chirp, 2); %fftshift(fft(eout, BlockSize, 2), 2);
    %fprintf('*** RDFFT data set 1 adaptive noise cenceller time : %d ***\n',exec_time_anc);

    RD_mat         = RDAfterANC_postive;
    [loc_sinr_anc,loc_sinr_range_anc,loc_sinr_doppler_anc ] = local_snr(RD_mat,peakIdx(1,:),n_sig_bin, n_noise_bin);
    sinr_anc = sinr_anc + loc_sinr_anc;

    RDAfterANC_full = [RDAfterANC_postive;RDFFT_data_sim_wi_int((ADC.count_sample /2)+1:end,:)];
    data_miti_mat_anc = ifft(ifft(RDAfterANC_full,ADC.count_chirp,2),ADC.count_sample ,1);
    corr_coef_anc = correlation_coeficient(beat_wo_int,reshape(data_miti_mat_anc,[numel(data_miti_mat_anc) 1]));
    corr_anc      = corr_anc + corr_coef_anc;

    time_anc = time_anc + exec_time_anc;

    if (j <= itr2)
        %Sparse and Henkel Matrix Decomposition
        addpath('./sparse_henkel_matrix_decomposition/lowRaS_MC/')
        disp('IM-SPARKLE...')
        beta_1 = 0.1;           % 0.5
        mu     = 0.02;          % 0.05
        tau    = 0.02;          % 0.02
        k_beta = 1.6;
        k_mu   = 1.2;
        R      = 10;
        tic;
        [data_miti_vec_sparkle, i_lowRaS, rerr] = lowRaS_Hankel(data_sim_wi_int, R, beta_1, mu, tau, k_beta, k_mu);

        exec_time_sparkle = toc;

        % fprintf('*** data set 1 sparse and hankel matrix decomposition time : %d ***\n',exec_time_sparkle);
        data_miti_mat_sparkle       = reshape(data_miti_vec_sparkle,[ADC.count_sample , ADC.count_chirp]);
        RangeFFTsparkle_data_miti   = fft(data_miti_mat_sparkle,ADC.count_sample,1);
        RDFFTsparkle_data_miti      = fft(RangeFFTsparkle_data_miti,ADC.count_chirp,2);

        RD_mat         = RDFFTsparkle_data_miti(1:ADC.count_sample/2,:) ;
        [loc_sinr_sparkle,loc_sinr_range_sparkle,loc_sinr_doppler_sparkle ] = local_snr(RD_mat,peakIdx(1,:),n_sig_bin, n_noise_bin);
        sinr_sparkle = sinr_sparkle + loc_sinr_sparkle;

        corr_coef_sparkle = correlation_coeficient(beat_wo_int,data_miti_vec_sparkle);
        corr_sparkle      = corr_sparkle + corr_coef_sparkle;

        time_sparkle = time_sparkle + exec_time_sparkle;

        %             figure;
        %             surf((20*log10(abs(RDFFTsparkle_data_miti))),'EdgeColor','none');%(1:2:end,1:2:end)
        %             hold on
        %             colormap default
    end






end

sinr_grad    = sinr_grad/itr;
sinr_int     = sinr_int/itr;
corr_grad    = corr_grad/itr ;
corr_int     = corr_int/itr;
time_grad    = time_grad/itr;

sinr_anc    = sinr_anc/itr;
sinr_sparkle     = sinr_sparkle/itr;
corr_anc    = corr_anc/itr ;
corr_sparkle     = corr_sparkle/itr;
time_anc    = time_anc/itr;
time_sparkle    = time_sparkle/itr;

disp('*********************************************');
fprintf('Interference Properties\n');
fprintf('Samples percentage: %s \n', num2str(int_percentage));
fprintf('SINR : %s \n', num2str(sinr_int));
fprintf('Correlation Coef. : %s \n', num2str(corr_int));
fprintf('Correlation Coef. in magn and angle : %s %s\n', num2str(abs(corr_int)),num2str(angle(corr_int)));

disp('*********************************************');
fprintf('Sub-gradient Descent\n');
fprintf('SINR : %s \n', num2str(sinr_grad));
fprintf('Correlation Coef. : %s\n', num2str(corr_grad));
fprintf('Correlation Coef. in magn and angle : %s %s\n', num2str(abs(corr_grad)),num2str(angle(corr_grad)));
fprintf('Execution Time : %s\n', num2str(time_grad));

disp('*********************************************');
fprintf('ANC\n');
fprintf('SINR : %s \n', num2str(sinr_anc));
fprintf('Correlation Coef. : %s\n', num2str(corr_anc));
fprintf('Correlation Coef. in magn and angle : %s %s\n', num2str(abs(corr_anc)),num2str(angle(corr_anc)));
fprintf('Execution Time : %s\n', num2str(time_anc));

disp('*********************************************');
fprintf('SPARKLE\n');
fprintf('SINR : %s \n', num2str(sinr_sparkle));
fprintf('Correlation Coef. : %s\n', num2str(corr_sparkle));
fprintf('Correlation Coef. in magn and angle : %s %s\n', num2str(abs(corr_sparkle)),num2str(angle(corr_sparkle)));
fprintf('Execution Time : %s\n', num2str(time_sparkle));


clear all;

sinr_int     = 0;
corr_int     = 0;

sinr_grad    = 0;
corr_grad    = 0;
time_grad    = 0;

sinr_anc    = 0;
corr_anc    = 0;
time_anc    = 0;

sinr_sparkle    = 0;
corr_sparkle    = 0;
time_sparkle    = 0;

itr = 1;
itr2 = 0;


for j = 1:itr


    ADC       = adcParam;      % victim radar's ADC specifications
    txchirp   = txchirpParam;  % victim radar's transmitted signal specifications
    intchirp  = intchirpParam2; % aggressor radar's transmitted signal specifications
    lowpaas   = filterParam;   % victim radar's lowpaas filter specifications

    %     target_range_list = 25;
    %     target_velocity_list = 3;
    %
    %     RCS_target_list = 0.2;
    %
    %     txAntennaGain_victim = 0.9;
    %     rxAntennaGain_victim = 0.9;
    %     txPower_victim = 1;
    %
    %     aggressor_range_list = 2;
    %     aggressor_velocity_list = 1;
    %     txAntennaGain_aggressor = 0.9;
    %     txPower_aggressor = 1;
    %
    %     [baseband_sig_reflections,baseband_sig_interferences] = scenario_based_adc_output...
    %     (ADC,txchirp,intchirp,lowpaas,...
    %     target_range_list,target_velocity_list,RCS_target_list,txAntennaGain_victim,rxAntennaGain_victim,txPower_victim,...
    %     aggressor_range_list,aggressor_velocity_list,txAntennaGain_aggressor,txPower_aggressor);

    %% data generation
    % reading parameter files


    % target
    target_range = 25;   % range of object
    target_velocity = 3; % velocity of object (+ve for velocity away from victim radar, -ve for velocity towards victim )
    % signal amplitudes
    tx_sig_amp   =   1;  % victim's transmitted signal amplitude
    rx_sig_amp   =   1;  % victim's received signal amplitude

    snr_db = -10; % SNR in dB

%     int_sig_amp  =  200; % amplitude of baseband interference signal after ADC in victim radar
%            
%     lambda = 5;
%     initial_step_size = 4;

%      int_sig_amp  =  600; % amplitude of baseband interference signal after ADC in victim radar
%            
%     lambda = 3;
%     initial_step_size = 2.4;

% %% 1.
%      int_sig_amp  =  1000; % amplitude of baseband interference signal after ADC in victim radar
%            
%      lambda = 5;
%      initial_step_size = 4;

%% 2.
      int_sig_amp  =  500; % amplitude of baseband interference signal after ADC in victim radar
           
     lambda = 2.5;
     initial_step_size =2;


    maxIter = 100;
    lowpaas.cutoff = 9e6;

    sir_db = 20*log10(tx_sig_amp*rx_sig_amp*(1/int_sig_amp))

    % baseband signal generation due to a singal object
    reflection_baseband_sig = tx_sig_amp*rx_sig_amp*baseband_reflection(target_range,target_velocity,txchirp,ADC,lowpaas);

    % aggressor location
    aggressor_range = 0;    % aggressor radar range from the victim (assumed to be at same(near) location)
    aggressor_velocity = 0; % aggressor's velocity

    % baseband signal generation due a single aggressor's transmitted signal in victim radar
    interference_baseband_sig =  int_sig_amp*baseband_interference(aggressor_range,aggressor_velocity,intchirp,txchirp,ADC,lowpaas);

%     figure;
%     imagesc(abs(reshape(interference_baseband_sig,[ ADC.count_sample ADC.count_chirp])))



    % interference quantity
    int_samples = double(abs(interference_baseband_sig) > 0);% number of non-zero samples in interference baseband signal
    int_percentage = sum(sum(int_samples))*(1/(ADC.count_sample*ADC.count_chirp))*100; % percentange of interference

    % beat signals
    beat_wo_int = reflection_baseband_sig;  % beat signal without interference
    beat_wi_int = reflection_baseband_sig + interference_baseband_sig; % beat signal with interference

    % Thermal noise generation
    %     snr_db = -60; % SNR in dB
    beat_power_wo_int  = (tx_sig_amp*rx_sig_amp)^2; %signal power
    noise_power_wo_int = beat_power_wo_int*(10^(-1*(snr_db/10))); % noise power
    cmplx_noise_wo_int = sqrt(noise_power_wo_int /2)*complex(randn(size(beat_wo_int )),randn(size(beat_wo_int ))); % noise signal

    beat_wo_int = beat_wo_int + cmplx_noise_wo_int;
   
    %     hold on
    %title('noise + interference')
    % beat_wo_int = beat_wo_int + cmplx_noise_wo_int; % beat signal without interference and with noise
    beat_wi_int = beat_wi_int + cmplx_noise_wo_int; % beat signal with interference and with noise

    %     beat_wi_int = beat_wi_int ;%%%%%%%%%%%

    % peak in Range-doppler spectrum grid
    grid_size = [ADC.count_sample, ADC.count_chirp]; % size of range-doppler grid
    peakIdx = peak_bin(target_range,target_velocity,grid_size,txchirp,ADC); % peak bin in Range-Doppler spectrum

    data_sim_wi_int             =   beat_wi_int; % data for processing
    data_sim_wi_int_mat         =   reshape(data_sim_wi_int,grid_size); % N(samples per chirp)xL(number of chirps)

%       freq_vs_time_plot(intchirp,txchirp,ADC,lowpaas,2,0,0)

    %     int = zeros(512,64);
    %     int(56:60,:) = 0.1*complex(randn(5,64),randn(5,64));
    %     data_sim_wi_int_mat =  data_sim_wi_int_mat +int; %%%%%%%%%

    RFFT_data_sim_wi_int        =   fft(data_sim_wi_int_mat,ADC.count_sample,1); % range FFT
    RDFFT_data_sim_wi_int       =   fft( RFFT_data_sim_wi_int,ADC.count_chirp,2); % Range-Doppler FFT

    RD_mat         = RDFFT_data_sim_wi_int(1:ADC.count_sample/2,:) ; % positive half of Range-Doppler spectrum
    %     figure;
    %     surf(10*log10(abs(RD_mat)),EdgeColor="none")

    n_sig_bin      = [2 ,2];  % number of signal bins at left and right of peak bin for range and doppler respectively
    n_noise_bin    = [6 ,6];  % number of noise bins at left and right of peak bin for range and doppler respectively
    % local SINR calculation at peak bin
    [loc_sinr_int,loc_sinr_range_int,loc_sinr_doppler_int ] = local_snr(RD_mat,peakIdx(1,:),n_sig_bin, n_noise_bin);
    sinr_int = sinr_int + loc_sinr_int;

    % correlation coefficient calculation
    corr_coef_int = correlation_coeficient(beat_wo_int,beat_wi_int);
    corr_int = corr_int + corr_coef_int;


    %% Interference Mitigation Alorithm
    %Sub-Gradient Descent
    %disp('Gradient Descent...')
    % Hyper parameters
    %     lambda = 0.4;
    %     initial_step_size = 0.5;
     limit = 0.01;
    %     maxIter = 100;
    %
    tic;
    [data_corr_mat_grad] = prox_sub_grad_des(data_sim_wi_int_mat,lambda,initial_step_size,limit,maxIter); %Interference matrix
    exec_time_grad_des = toc;

%      figure;
%     plot(real(cmplx_noise_wo_int(1:512) + real(interference_baseband_sig(1:512))),'b')
%     hold on
%     plot(real(data_corr_mat_grad(:,1)),'r')

    %     plot(real(data_corr_mat_grad(:,1)) ,'r');
    %     title('estimated intereference')
    %     figure;
    %     plot(real(interference_baseband_sig(1:512))-real(data_corr_mat_grad(:,1)))
    %     title('interference - estimated interference')
    %     hold on
    %      plot(real(int(:,1)),'b');
    %     figure;
    %     plot(real(data_sim_wi_int_mat(:,1))-real(data_corr_mat_grad(:,1)));
    %     figure;
    %     plot(real(data_sim_wi_int_mat(:,1)));
    %fprintf('*** data set 1 gradient descent time : %d ***\n',exec_time_grad_des);
    data_miti_mat_grad       =   data_sim_wi_int_mat - data_corr_mat_grad; % interefernce cancelled matrix
    RFFTGrad_data_miti       =   fft(data_miti_mat_grad,ADC.count_sample,1); % Range fft
    RDFFTGrad_data_miti      =   fft(RFFTGrad_data_miti,ADC.count_chirp,2); % Range-Doopler fft

    RD_mat         = RDFFTGrad_data_miti(1:ADC.count_sample/2,:); % positive half of RD FFT
    [loc_sinr_grad,loc_sinr_range_grad,loc_sinr_doppler_grad ] = local_snr(RD_mat,peakIdx(1,:),n_sig_bin, n_noise_bin); % SINR computation after interference cancellation
    sinr_grad = sinr_grad + loc_sinr_grad;

    % correlation coefficient
    corr_coef_grad = correlation_coeficient(beat_wo_int,reshape(data_miti_mat_grad,[numel(data_miti_mat_grad) 1]));
    corr_grad = corr_grad + corr_coef_grad;

    % execution time
    time_grad = time_grad + exec_time_grad_des;



    %Adaptive Noise canceller
    addpath('./adaptive_noise_canceller/')
    %disp('Adaptive Noise Canceller...')
    flen = 8; % lenght of adaptive filter 50 for 256
    tic;
    % anc_RDFFT_miti = anc_RDFFT(RFFT_data,NumSamplePerChirp,BlockSize,flen);
    eout = anc_RDFFT(RFFT_data_sim_wi_int,ADC.count_sample,ADC.count_chirp,flen);
    exec_time_anc  = toc;

    RDAfterANC_postive      = fft(eout,ADC.count_chirp, 2); %fftshift(fft(eout, BlockSize, 2), 2);
    %fprintf('*** RDFFT data set 1 adaptive noise cenceller time : %d ***\n',exec_time_anc);

    RD_mat         = RDAfterANC_postive;
    [loc_sinr_anc,loc_sinr_range_anc,loc_sinr_doppler_anc ] = local_snr(RD_mat,peakIdx(1,:),n_sig_bin, n_noise_bin);
    sinr_anc = sinr_anc + loc_sinr_anc;

    RDAfterANC_full = [RDAfterANC_postive;RDFFT_data_sim_wi_int((ADC.count_sample /2)+1:end,:)];
    data_miti_mat_anc = ifft(ifft(RDAfterANC_full,ADC.count_chirp,2),ADC.count_sample ,1);
    corr_coef_anc = correlation_coeficient(beat_wo_int,reshape(data_miti_mat_anc,[numel(data_miti_mat_anc) 1]));
    corr_anc      = corr_anc + corr_coef_anc;

    time_anc = time_anc + exec_time_anc;

    if (j <= itr2)
        %Sparse and Henkel Matrix Decomposition
        addpath('./sparse_henkel_matrix_decomposition/lowRaS_MC/')
        disp('IM-SPARKLE...')
        beta_1 = 0.1;           % 0.5
        mu     = 0.02;          % 0.05
        tau    = 0.02;          % 0.02
        k_beta = 1.6;
        k_mu   = 1.2;
        R      = 10;
        tic;
        [data_miti_vec_sparkle, i_lowRaS, rerr] = lowRaS_Hankel(data_sim_wi_int, R, beta_1, mu, tau, k_beta, k_mu);

        exec_time_sparkle = toc;

        % fprintf('*** data set 1 sparse and hankel matrix decomposition time : %d ***\n',exec_time_sparkle);
        data_miti_mat_sparkle       = reshape(data_miti_vec_sparkle,[ADC.count_sample , ADC.count_chirp]);
        RangeFFTsparkle_data_miti   = fft(data_miti_mat_sparkle,ADC.count_sample,1);
        RDFFTsparkle_data_miti      = fft(RangeFFTsparkle_data_miti,ADC.count_chirp,2);

        RD_mat         = RDFFTsparkle_data_miti(1:ADC.count_sample/2,:) ;
        [loc_sinr_sparkle,loc_sinr_range_sparkle,loc_sinr_doppler_sparkle ] = local_snr(RD_mat,peakIdx(1,:),n_sig_bin, n_noise_bin);
        sinr_sparkle = sinr_sparkle + loc_sinr_sparkle;

        corr_coef_sparkle = correlation_coeficient(beat_wo_int,data_miti_vec_sparkle);
        corr_sparkle      = corr_sparkle + corr_coef_sparkle;

        time_sparkle = time_sparkle + exec_time_sparkle;

        %             figure;
        %             surf((20*log10(abs(RDFFTsparkle_data_miti))),'EdgeColor','none');%(1:2:end,1:2:end)
        %             hold on
        %             colormap default
    end






end

sinr_grad    = sinr_grad/itr;
sinr_int     = sinr_int/itr;
corr_grad    = corr_grad/itr ;
corr_int     = corr_int/itr;
time_grad    = time_grad/itr;

sinr_anc    = sinr_anc/itr;
sinr_sparkle     = sinr_sparkle/itr;
corr_anc    = corr_anc/itr ;
corr_sparkle     = corr_sparkle/itr;
time_anc    = time_anc/itr;
time_sparkle    = time_sparkle/itr;

disp('*********************************************');
fprintf('Interference Properties\n');
fprintf('Samples percentage: %s \n', num2str(int_percentage));
fprintf('SINR : %s \n', num2str(sinr_int));
fprintf('Correlation Coef. : %s \n', num2str(corr_int));
fprintf('Correlation Coef. in magn and angle : %s %s\n', num2str(abs(corr_int)),num2str(angle(corr_int)));

disp('*********************************************');
fprintf('Sub-gradient Descent\n');
fprintf('SINR : %s \n', num2str(sinr_grad));
fprintf('Correlation Coef. : %s\n', num2str(corr_grad));
fprintf('Correlation Coef. in magn and angle : %s %s\n', num2str(abs(corr_grad)),num2str(angle(corr_grad)));
fprintf('Execution Time : %s\n', num2str(time_grad));

disp('*********************************************');
fprintf('ANC\n');
fprintf('SINR : %s \n', num2str(sinr_anc));
fprintf('Correlation Coef. : %s\n', num2str(corr_anc));
fprintf('Correlation Coef. in magn and angle : %s %s\n', num2str(abs(corr_anc)),num2str(angle(corr_anc)));
fprintf('Execution Time : %s\n', num2str(time_anc));

disp('*********************************************');
fprintf('SPARKLE\n');
fprintf('SINR : %s \n', num2str(sinr_sparkle));
fprintf('Correlation Coef. : %s\n', num2str(corr_sparkle));
fprintf('Correlation Coef. in magn and angle : %s %s\n', num2str(abs(corr_sparkle)),num2str(angle(corr_sparkle)));
fprintf('Execution Time : %s\n', num2str(time_sparkle));


