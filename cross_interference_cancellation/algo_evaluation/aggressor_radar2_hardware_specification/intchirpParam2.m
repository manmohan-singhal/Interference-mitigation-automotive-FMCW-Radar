function intchirp = intchirpParam2

intchirp.count             = 64;
intchirp.basefreq          = 76.98e9;%%
intchirp.slope             = 45e12;%%
intchirp.duration          = 50e-6;%% 
intchirp.bandwidth         = intchirp.slope*intchirp.duration;
intchirp.repetetion_period = 159.95e-6;
intchirp.start_time        = 3e-6;


