function [baseband_sig_reflections,baseband_sig_interferences] = scenario_based_adc_output...
    (ADC,txchirp,intchirp,lowpaas,...
    target_range_list,target_velocity_list,RCS_target_list,txAntennaGain_victim,rxAntennaGain_victim,txPower_victim,...
    aggressor_range_list,aggressor_velocity_list,txAntennaGain_aggressor,txPower_aggressor)


baseband_sig_reflections = zeros(ADC.count_sample*ADC.count_chirp,1);
baseband_sig_interferences = zeros(ADC.count_sample*ADC.count_chirp,1);

centerfreq_victim = txchirp.bandwidth/2 + txchirp.basefreq;
tx_sig_amp_victim = sqrt(txPower_victim);

[m,n] = size(target_range_list);
k = max(m,n);
if k ~=0
    for i = 1:k
        rx_Power_reflection = reflection_power(target_range_list(i),txAntennaGain_victim,rxAntennaGain_victim,...
            txPower_victim,centerfreq_victim,RCS_target_list(i));
        rx_sig_amp_reflection = sqrt(rx_Power_reflection);
        % baseband signal generation due to a singal object
        reflection_baseband_sig = tx_sig_amp_victim*rx_sig_amp_reflection*baseband_reflection...
            (target_range_list(i),target_velocity_list(i),txchirp,ADC,lowpaas);
        baseband_sig_reflections = baseband_sig_reflections + reflection_baseband_sig;
    end
end

centerfreq_aggressor = intchirp.bandwidth/2 + intchirp.basefreq;
tx_sig_amp_aggressor = sqrt(txPower_aggressor);

[m,n] = size(aggressor_range_list);
k = max(m,n);
if k~=0
    for i = 1:k
        rx_Power_interference = interference_power(aggressor_range_list(i),txAntennaGain_aggressor,...
            rxAntennaGain_victim,txPower_aggressor,centerfreq_aggressor);
        rx_sig_amp_interference = sqrt(rx_Power_interference);

        interference_baseband_sig =  tx_sig_amp_aggressor*rx_sig_amp_interference*baseband_interference...
            (aggressor_range_list(i),aggressor_velocity_list,intchirp,txchirp,ADC,lowpaas);
        baseband_sig_interferences = baseband_sig_interferences + interference_baseband_sig;
   
    end


end



