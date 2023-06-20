function intchirp = intchirpParam1

intchirp.count             = 64;
intchirp.basefreq          = 76.6e9;
intchirp.slope             = 60e12;
intchirp.duration          = 40e-6;
intchirp.bandwidth         = intchirp.slope*intchirp.duration;
intchirp.repetetion_period = 159.95e-6 ;%- 5e-9;
intchirp.start_time        = 3e-6;