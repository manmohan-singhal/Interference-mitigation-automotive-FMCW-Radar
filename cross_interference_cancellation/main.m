close all;
clear all;
clc;

ADC = adcParam;
txchirp = txchirpParam;
intchirp = intchirpParam;
lowpaas = filterParam;
target_range = 32;
target_velocity = 4;
aggressor_range = 0;
aggressor_velocity = 0;

txAntennaGain_victim = 0.9;
rxAntennaGain_victim = 0.9;
txPower_victim = 1;
centerfreq_victim = txchirp.basefreq + txchirp.slope*(txchirp.duration/2);
RCS_object = 0.01;

txAntennaGain_aggressor= 0.9;
txPower_aggressor = 1;
centerfreq_aggressor = intchirp.basefreq + intchirp.slope*(intchirp.duration/2);

rx_power_reflection = reflection_power(target_range,txAntennaGain_victim,rxAntennaGain_victim,txPower_victim,centerfreq_victim,RCS_object);
reflection_baseband_sig = rx_power_reflection*baseband_reflection(target_range,target_velocity,txchirp,ADC,lowpaas);

rx_Power_interference = interference_power(aggressor_range,txAntennaGain_aggressor,rxAntennaGain_victim,txPower_aggressor,centerfreq_aggressor);
interference_baseband_sig =  rx_Power_interference*baseband_interference(aggressor_range,aggressor_velocity,intchirp,txchirp,ADC,lowpaas);

x = reflection_baseband_sig + interference_baseband_sig;

X = reshape(x,[256 64]);
Y = fft(X);
Z = fft(transpose(Y));
figure;
surf(abs(Z),'EdgeColor','none');

nchirps = 5;
freqVStime(intchirp,txchirp,ADC,lowpaas,nchirps,aggressor_range,aggressor_velocity);
