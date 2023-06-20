function range_doppler_map(Range_doppler_DFT,peak_list,ADC,txchirp)

 % reading parameter files
%  ADC       = adcParam;      % victim radar's ADC specifications
%  txchirp   = txchirpParam;  % victim radar's transmitted signal specificationS

chirp_duration = txchirp.duration  ;
B              = txchirp.bandwidth ;
fc             = txchirp.basefreq ;
TRRI           = txchirp.repetetion_period ;
chirps         = ADC.count_chirp ;

Ts             = ADC.period;
N_samples      = ADC.count_sample  ;


range_axis =  1*((-N_samples/2:N_samples/2-1)*(chirp_duration*3e8/(2*B*Ts*N_samples))); % mapping x axis

doppler_axis =  1*(-chirps/2:chirps/2-1)*(3e8/(2*TRRI*fc*chirps)); % mapping y axis


axes1 = axes('Parent',figure);
hold(axes1,'on');
surf(range_axis,doppler_axis,(20*log10(abs(fftshift(Range_doppler_DFT)))),'EdgeColor','none');%(1:2:end,1:2:end)
%surf(range_axis,doppler_axis,(20*log10(abs(Range_doppler_DFT))),'EdgeColor','none');%(1:2:end,1:2:end)
xlabel('range axis (m)','FontSize',18)
ylabel('doppler axis (m/s)','FontSize',18)
zlabel('magnitude axis (dB)','FontSize',18)
% title('Range-Doppler Map','FontSize',18)

view(axes1,[153.847798683667 39.9236569097532]);
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'GridAlpha',1,'GridColor',...
    [0.650980392156863 0.650980392156863 0.650980392156863],'LineWidth',1,...
    'MinorGridColor',[1 1 1]);



xh = get(gca,'XLabel'); % Handle of the x label
set(xh, 'Units', 'Normalized')
pos = get(xh, 'Position');
set(xh, 'Position',pos.*[1,1,1],'Rotation',11)
yh = get(gca,'YLabel'); % Handle of the y label
set(yh, 'Units', 'Normalized')
pos = get(yh, 'Position');
set(yh, 'Position',pos.*[1,1,1],'Rotation',-11)
zh = get(gca,'ZLabel'); % Handle of the z label
set(zh, 'Units', 'Normalized')
pos = get(zh, 'Position');
set(zh, 'Position',pos.*[1,1,0],'Rotation',90)
hold on

colormap default


if size(peak_list) ~= [0,0]
    h = scatter3(range_axis(peak_list(:,1)),doppler_axis(peak_list(:,2)),20*log10(peak_list(:,3)),'filled');
    % h = scatter3(-1*Range_indices*(3e8/(2*TRRI*fc*chirps)),-1*Doppler_indices*(chirp_duration*3e8/(2*B*Ts*N_samples)),ones(n),'filled');
    h.MarkerFaceColor = [1 0 0];
end


end