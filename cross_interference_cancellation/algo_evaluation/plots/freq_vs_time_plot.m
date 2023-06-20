function freq_vs_time_plot(intchirp,txchirp,ADC,lowpaas,nchirps,aggressor_distance,aggressor_velocity)

c = 3e8;
% generate time stamps vector

time_stamps = ADC.start_time:ADC.period:nchirps*(txchirp.repetetion_period);
[~,m] = size(time_stamps);

% delay vector for each sample
delay_vector = (aggressor_distance - aggressor_velocity*time_stamps)*(1/c);

txfreq = zeros(m,1);
for idx = 1:m
    %    chirp_idx = floor((time_stamps(idx)-ADC.start_time)*(1/txchirp.repetetion_time));
    chirp_instant = rem((time_stamps(idx)-txchirp.start_time),txchirp.repetetion_period);
    if chirp_instant <= txchirp.duration
        txfreq(idx) = txchirp.basefreq + (txchirp.bandwidth*(1/txchirp.duration)*(chirp_instant)) ;
    else
        txfreq(idx) = txchirp.basefreq;
    end

end


intfreq = zeros(m,1);
for idx = 1:m
    %    chirp_idx = floor((time_stamps(idx)-ADC.start_time)*(1/intchirp.repetetion_time));
    chirp_instant = rem((time_stamps(idx)-intchirp.start_time),intchirp.repetetion_period);
    if chirp_instant <= intchirp.duration
        intfreq(idx) = intchirp.basefreq + (intchirp.bandwidth*(1/intchirp.duration)*(chirp_instant - delay_vector(idx))) ;
    else
        intfreq(idx) = intchirp.basefreq;
    end

end

upper_lowpaas_boundary = txfreq + lowpaas.cutoff;
lower_lowpaas_boundary = txfreq - lowpaas.cutoff;

time_stamps = time_stamps*1e6;
txfreq = txfreq*1e-9;
intfreq = intfreq*1e-9;

figure;
plot(time_stamps,txfreq,LineWidth=1.5);
hold on 
plot(time_stamps,intfreq,LineWidth=1.5);
hold on 

xlabel('Time (\mu s)','FontSize',18)
ylabel('Frequency (GHz)','FontSize',18)
% % plot(time_stamps,upper_lowpaas_boundary,'Color','r')
% % hold on 
% % plot(time_stamps,lower_lowpaas_boundary,'Color','r')
% % hold on
% for idx = 0:nchirps-1
%     xline(ADC.start_time + ADC.repetition_time*idx,'Color','g','LineWidth',2)
%     hold on
%    xline(ADC.start_time+ ADC.period*ADC.count_sample + ADC.repetition_time*idx,'Color','g','LineWidth',2)
%    hold on
% end



