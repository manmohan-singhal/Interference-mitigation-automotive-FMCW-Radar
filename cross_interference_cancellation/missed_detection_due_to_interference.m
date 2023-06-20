close all; clear all; clc;

addpath('peak_detection\')
addpath('plots\')
addpath('hardware_specification_read_files\')
addpath('baseband_signal_model\')
ADC       = adcParam;      % victim radar's ADC specifications
lowpaas   = filterParam;   % victim radar's lowpaas filter specificatio

%% sweeping interference
txchirp   = txchirpParam;  % victim radar's transmitted signal specifications
intchirp  = intchirpParam; % aggressor radar's transmitted signal specifications


% signal amplitude
tx_sig_amp   =   1;  % victim's transmitted signal amplitude
rx_sig_amp   =   1;  % victim's received signal amplitude
int_sig_amp  =   30;

target_range_1 = 10;
target_range_2 = 20;
target_velocity_1 = 3;
target_velocity_2 = 4;

% aggressor location
aggressor_range    = 0;    % aggressor radar range from the victim (assumed to be at same(near) location)
aggressor_velocity = 0; % aggressor's velocity


% baseband signal generation due to a single object
reflection_baseband_sig1 = tx_sig_amp*rx_sig_amp*baseband_reflection(target_range_1,target_velocity_1,txchirp,ADC,lowpaas);
reflection_baseband_sig2 = (1/4)*tx_sig_amp*rx_sig_amp*baseband_reflection(target_range_2,target_velocity_2,txchirp,ADC,lowpaas);
reflection_baseband_sig  = reflection_baseband_sig1+ reflection_baseband_sig2;


% baseband signal generation due a single aggressor's transmitted signal in victim radar
interference_baseband_sig =  int_sig_amp*baseband_interference(aggressor_range,aggressor_velocity,intchirp,txchirp,ADC,lowpaas);

% interference quantity
int_samples = double(abs(interference_baseband_sig) > 0);% number of non-zero samples in interference baseband signal
int_percentage = sum(sum(int_samples))*(1/(ADC.count_sample*ADC.count_chirp))*100 % percentange of interference
%
% beat signals
beat_wo_int = reflection_baseband_sig;  % beat signal without interference
beat_wi_int = reflection_baseband_sig + interference_baseband_sig; % beat signal with interference

% Thermal noise generation
snr_db = -20; % SNR in dB
beat_power_wo_int  = ((tx_sig_amp*rx_sig_amp)^2 + ((1/2)*tx_sig_amp*rx_sig_amp)^2)/2; %signal power
noise_power_wo_int = beat_power_wo_int*(10^(-1*(snr_db/10))); % noise power
cmplx_noise_wo_int = sqrt(noise_power_wo_int /2)*complex(randn(size(beat_wo_int )),randn(size(beat_wo_int ))); % noise signal

beat_wo_int = beat_wo_int + cmplx_noise_wo_int; % beat signal without interference and with noise
% beat_wi_int = beat_wi_int + cmplx_noise_wo_int; % beat signal with interference and with noise

% peak in Range-doppler spectrum grid
grid_size = [ADC.count_sample, ADC.count_chirp]; % size of range-doppler grid
%peakIdx = peak_bin(target_range,target_velocity,grid_size,txchirp,ADC); % peak bin in Range-Doppler spectrum

data_sim             =   beat_wo_int + beat_wi_int; % data for processing
data_sim_mat         =   reshape(data_sim,grid_size); % N(samples per chirp)xL(number of chirps)
RFFT_data_sim        =   fft(data_sim_mat,ADC.count_sample,1); % range FFT
RDFFT_data_sim       =   fft( RFFT_data_sim,ADC.count_chirp,2); % Range-Doppler FFT

cell.NumTrainingCellRange = 5;
cell.NumTrainingCellDoppler = 5;
cell.NumGuardCellRange = 2;
cell.NumGuardCellDoppler = 2;

Threshold_dB = 16;
Threshold = 10^(Threshold_dB/10);

peak_bins_rdm = constant_false_alarm_rate(cell,fftshift(transpose(RDFFT_data_sim)),Threshold);

%[~,index] = sort(peak_bins_rdm(:,3),'descend');

range_doppler_map(transpose(RDFFT_data_sim),peak_bins_rdm,ADC,txchirp)
