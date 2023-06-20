close all;
clear all;
clc;

% NumChirpsPerFrame     = 128;RDFFT_data_sim_wo_intt
NumSamplePerChirp     = 512;
%% DATA
StartChirpIdx         = 0;      % for processing, start from StartChirpIdx + 1
EndChirpIdx           = 64;     % for processing
BlockSize             = EndChirpIdx - StartChirpIdx;

Int_percentage = [];

Sinr_grad      = [];
Sinr_anc       = [];
Sinr_sparkle   = [];
Loc_sinr_int   = [];

sinr_anc     = 0;
sinr_grad    = 0;

sinr_sparkle = 0;
sinr_int     = 0;

Corr_grad      = [];
Corr_anc       = [];
Corr_sparkle   = [];
Corr_int   = [];

corr_anc     = 0;
corr_grad    = 0;

corr_sparkle = 0;
corr_int     = 0;

Time_grad      = [];
Time_anc       = [];
Time_sparkle   = [];

time_anc     = 0;
time_grad    = 0;
time_sparkle = 0;

itr  = 1;
itr1 = 100;
itr2 = 0;
for i= 1:itr
    for j = 1:itr1

        % SIMULATED
        % with interference
        ADC = adcParam;
        txchirp = txchirpParam;
        intchirp = intchirpParam;
        lowpaas = filterParam;

        target_range = 32;
        target_velocity = 4;

        tx_sig_amp   =   1;
        rx_sig_amp   =   1;

        reflection_baseband_sig = baseband_reflection(target_range,target_velocity,txchirp,ADC,lowpaas);

        aggressor_range = 0;
        aggressor_velocity = 0;

        int_sig_amp   =   50;
        interference_baseband_sig =  int_sig_amp*baseband_interference(aggressor_range,aggressor_velocity,intchirp,txchirp,ADC,lowpaas);

        int_per = double(abs(interference_baseband_sig) > 0);%int_beat1 +int_beat2 +
        int_percentage = sum(sum(int_per))*(1/(ADC.count_sample*ADC.count_chirp))*100;

        beat_wo_int = reflection_baseband_sig;
        beat_wi_int = reflection_baseband_sig + interference_baseband_sig;
        snr_db = 6;

        beat_power_wo_int  = sum(abs(beat_wo_int (:)).^2)/numel(beat_wo_int ); %signal power
        noise_power_wo_int = beat_power_wo_int*(10^(-1*(snr_db/10)));
        cmplx_noise_wo_int = sqrt(noise_power_wo_int /2)*complex(randn(size(beat_wo_int )),randn(size(beat_wo_int )));

        beat_wo_int = beat_wo_int + cmplx_noise_wo_int;
        beat_wi_int = beat_wi_int + cmplx_noise_wo_int;

        grid_size = [512, 64];
        peakIdx = peak_bin(target_range,target_velocity,grid_size,txchirp,ADC);

        path4figure = './data/';
        save( [path4figure 'data' num2str(i)  'wo_int_.mat'],'beat_wo_int');
        save( [path4figure 'data' num2str(i)  'wi_int_.mat'],'beat_wi_int');

        sig_bin = peakIdx;
        data_sim_wi_int             =   beat_wi_int(1+NumSamplePerChirp*StartChirpIdx:NumSamplePerChirp*EndChirpIdx);
        data_sim_wi_int_mat         =   reshape(data_sim_wi_int,[NumSamplePerChirp,(EndChirpIdx-StartChirpIdx)]); % N(samples per chirp)xL(number of chirps)
        RFFT_data_sim_wi_int        =   fft(data_sim_wi_int_mat,NumSamplePerChirp,1);
        RDFFT_data_sim_wi_int       =   fft(RFFT_data_sim_wi_int,BlockSize,2);

        RD_mat         = RDFFT_data_sim_wi_int(1:NumSamplePerChirp/2,:) ;

        n_sig_bin      = [2 ,2];  % 1st entry for range direction and 2nd entry is for doppler direction
        n_noise_bin    = [6 ,6]; % 1st entry for range direction and 2nd entry is for doppler direction
        [loc_sinr_int,loc_sinr_range_int,loc_sinr_doppler_int ] = local_snr(RD_mat,sig_bin(1,:),n_sig_bin, n_noise_bin);
        sinr_int = sinr_int + loc_sinr_int;

        corr_coef_int = correlation_coeficient(beat_wo_int,beat_wi_int);
        corr_int = corr_int + corr_coef_int;


        %% Interference Mitigation Alorithms


        %Sub-Gradient Descent
        %disp('Gradient Descent...')
        addpath('./gradient_descent/')
        % Hyper parameters
        lambda = 0.4;
        initial_step_size = 0.5;
        limit = 0.1;
        maxIter = 128;
        %
        tic;
        [data_corr_mat_grad] = sub_grad_des_cmplx_mat_updates_v3(data_sim_wi_int_mat,lambda,initial_step_size,limit,maxIter);
        exec_time_grad_des = toc;
        %fprintf('*** data set 1 gradient descent time : %d ***\n',exec_time_grad_des);
        data_miti_mat_grad       =   data_sim_wi_int_mat + data_corr_mat_grad;
        RFFTGrad_data_miti       =   fft(data_miti_mat_grad,NumSamplePerChirp,1);
        RDFFTGrad_data_miti      =   fft(RFFTGrad_data_miti,BlockSize,2);

        RD_mat         = RDFFTGrad_data_miti(1:NumSamplePerChirp/2,:);
        [loc_sinr_grad,loc_sinr_range_grad,loc_sinr_doppler_grad ] = local_snr(RD_mat,sig_bin(1,:),n_sig_bin, n_noise_bin);
        sinr_grad = sinr_grad + loc_sinr_grad;

        corr_coef_grad = correlation_coeficient(beat_wo_int,reshape(data_miti_mat_grad,[numel(data_miti_mat_grad) 1]));
        corr_grad = corr_grad + corr_coef_grad;

        time_grad = time_grad + exec_time_grad_des;


% %         %Adaptive Noise canceller
% %         addpath('./adaptive_noise_canceller/')
% %         %disp('Adaptive Noise Canceller...')
% %         flen = 4; % lenght of adaptive filter 50 for 256
% %         tic;
% %         % anc_RDFFT_miti = anc_RDFFT(RFFT_data,NumSamplePerChirp,BlockSize,flen);
% %         eout = anc_RDFFT(RFFT_data_sim_wi_int,NumSamplePerChirp,BlockSize,flen);
% %         exec_time_anc  = toc;
% %         path4figure = './data/';
% %         save( [path4figure 'anc_output' num2str(i)  'mit_.mat'],'eout');
% %         RDAfterANC_postive      = fft(eout, BlockSize, 2); %fftshift(fft(eout, BlockSize, 2), 2);
% %         %fprintf('*** RDFFT data set 1 adaptive noise cenceller time : %d ***\n',exec_time_anc);
% % 
% %         RD_mat         = RDAfterANC_postive;
% %         [loc_sinr_anc,loc_sinr_range_anc,loc_sinr_doppler_anc ] = local_snr(RD_mat,sig_bin(1,:),n_sig_bin, n_noise_bin);
% %         sinr_anc = sinr_anc + loc_sinr_anc;
% % 
% %         RDAfterANC_full = [RDAfterANC_postive;RDFFT_data_sim_wi_int((NumSamplePerChirp/2)+1:end,:)];
% %         data_miti_mat_anc = ifft(ifft(RDAfterANC_full,BlockSize,2),NumSamplePerChirp,1);
% %         corr_coef_anc = correlation_coeficient(beat_wo_int,reshape(data_miti_mat_anc,[numel(data_miti_mat_anc) 1]));
% %         corr_anc      = corr_anc + corr_coef_anc;
% % 
% %         time_anc = time_anc + exec_time_anc;
% % % 
% % %                 figure;
% % %                 surf((20*log10(abs(RDAfterANC_postive))),'EdgeColor','none');%(1:2:end,1:2:end)
% % %                 hold on
% % %                 colormap default
% % %                 figure;
% % %                 surf((20*log10(abs(RDFFT_data_sim_wi_int))),'EdgeColor','none');%(1:2:end,1:2:end)
% % %                 hold on
% % %                 colormap default
% % % 
% % % 
% % %         RangeFFTANC_data_miti = RFFT_data;
% % %         RangeFFTANC_data_miti(1:NumSamplePerChirp/2,:) = anc_RFFT_miti_positive;
% % %         RDFFTANC_data_miti      = fft(RangeFFTANC_data_miti,BlockSize,2);

% %         if (j <= itr2)
% %             %Sparse and Henkel Matrix Decomposition
% %             addpath('./sparse_henkel_matrix_decomposition/lowRaS_MC/')
% %             disp('IM-SPARKLE...')
% %             beta_1 = 0.1;           % 0.5
% %             mu     = 0.02;          % 0.05
% %             tau    = 0.02;          % 0.02
% %             k_beta = 1.6;
% %             k_mu   = 1.2;
% %             R = 10;
% %             tic;
% %             [data_miti_vec_sparkle, i_lowRaS, rerr] = lowRaS_Hankel(data_sim_wi_int, R, beta_1, mu, tau, k_beta, k_mu);
% % 
% %             exec_time_sparkle = toc;
% % 
% %             path4figure = './data/';
% %             save( [path4figure 'sparkle_output' num2str(i)  'mit_.mat'],'data_miti_vec_sparkle');
% %             save( [path4figure 'sparkle_output' num2str(i)  'ilow_.mat'],'i_lowRaS');
% %             save( [path4figure 'sparkle_output' num2str(i)  'rerr_.mat'],'rerr');
% %             % fprintf('*** data set 1 sparse and hankel matrix decomposition time : %d ***\n',exec_time_sparkle);
% %             data_miti_mat_sparkle       = reshape(data_miti_vec_sparkle,[NumSamplePerChirp , BlockSize]);
% %             RangeFFTsparkle_data_miti   = fft(data_miti_mat_sparkle,NumSamplePerChirp,1);
% %             RDFFTsparkle_data_miti      = fft(RangeFFTsparkle_data_miti,BlockSize,2);
% % 
% %             RD_mat         = RDFFTsparkle_data_miti(1:NumSamplePerChirp/2,:) ;
% %             [loc_sinr_sparkle,loc_sinr_range_sparkle,loc_sinr_doppler_sparkle ] = local_snr(RD_mat,sig_bin(1,:),n_sig_bin, n_noise_bin);;
% %             sinr_sparkle = sinr_sparkle + loc_sinr_sparkle;
% % 
% %             corr_coef_sparkle = correlation_coeficient(beat_wo_int,data_miti_vec_sparkle);
% %             corr_sparkle      = corr_sparkle + corr_coef_sparkle;
% % 
% %             time_sparkle = time_sparkle + exec_time_sparkle;
% % 
% %             %             figure;
% %             %             surf((20*log10(abs(RDFFTsparkle_data_miti))),'EdgeColor','none');%(1:2:end,1:2:end)
% %             %             hold on
% %             %             colormap default
% %         end
        %% Performance evaluation and  comparison

        % figure();
        % % surf(range_axis,doppler_axis,(20*log10(abs(fftshift(Range_doppler_DFT)))),'EdgeColor','none');%(1:2:end,1:2:end)
        % surf((20*log10(abs(RDFFT_data_sim_wo_int ))),'EdgeColor','none');%(1:2:end,1:2:end)
        % hold on
        % colormap default
        % figure;
        % % surf(range_axis,doppler_axis,(20*log10(abs(fftshift(Range_doppler_DFT)))),'EdgeColor','none');%(1:2:end,1:2:end)
        % surf((20*log10(abs(RDFFT_data_sim_wi_int ))),'EdgeColor','none');%(1:2:end,1:2:end)
        % hold on
        % colormap default


        %sig_bin = peak_idx_r_d (RDFFT_data_sim_wo_int(1:NumSamplePerChirp/2,:));

    end

    Sinr_grad = [Sinr_grad , sinr_grad/itr1 ];
    Sinr_anc = [Sinr_anc, sinr_anc/itr1];
    Sinr_sparkle = [Sinr_sparkle sinr_sparkle/itr2];
    Loc_sinr_int = [Loc_sinr_int sinr_int/itr1];

    sinr_anc     = 0;
    sinr_grad    = 0;

    sinr_sparkle = 0;
    sinr_int     = 0;


    Corr_grad = [Corr_grad , corr_grad/itr1 ];
    Corr_anc = [Corr_anc, corr_anc/itr1];
    Corr_sparkle = [Corr_sparkle corr_sparkle/itr2];
    Corr_int = [Corr_int corr_int/itr1];

    corr_anc     = 0;
    corr_grad    = 0;

    corr_sparkle = 0;
    corr_int     = 0;


    Time_grad      = [Time_grad    time_grad/itr1];
    Time_anc       = [Time_anc    time_anc/itr1];
    Time_sparkle   = [Time_sparkle  time_sparkle/itr2];

    time_anc     = 0;
    time_grad    = 0;

    time_sparkle = 0;

    Int_percentage = [Int_percentage, int_percentage];


end

disp('*********************************************');
fprintf('Interference Properties\n');
fprintf('Samples percentage: %s \n', num2str(Int_percentage));
fprintf('SINR : %s \n', num2str(Loc_sinr_int));
fprintf('Correlation Coef. : %s \n', num2str(Corr_int));
fprintf('Correlation Coef. in magn and angle : %s %s\n', num2str(abs(Corr_int)),num2str(angle(Corr_int)));

disp('*********************************************');
fprintf('Sub-gradient Descent\n');
fprintf('SINR : %s \n', num2str(Sinr_grad));
fprintf('Correlation Coef. : %s\n', num2str(Corr_grad));
fprintf('Correlation Coef. in magn and angle : %s %s\n', num2str(abs(Corr_grad)),num2str(angle(Corr_grad)));
fprintf('Execution Time : %s\n', num2str(Time_grad));

disp('*********************************************');
fprintf('Adaptive Noise Canceller\n');
fprintf('SINR : %s \n', num2str(Sinr_anc));
fprintf('Correlation Coef. : %s \n', num2str(Corr_anc));
fprintf('Correlation Coef. in magn and angle : %s %s\n', num2str(abs(Corr_anc)),num2str(angle(Corr_anc)));
fprintf('Execution Time : %s \n', num2str(Time_anc));
disp('*********************************************');

fprintf('Sparse and Low rank Hankel Matrix Decomposition\n');
fprintf('SINR : %s \n', num2str(Sinr_sparkle));
fprintf('Correlation Coef.: %s \n', num2str(Corr_sparkle));
fprintf('Correlation Coef. in magn and angle : %s %s\n', num2str(abs(Corr_sparkle)),num2str(angle(corr_sparkle)));
fprintf('Execution Time : %s\n', num2str(Time_sparkle));
disp('*********************************************');



