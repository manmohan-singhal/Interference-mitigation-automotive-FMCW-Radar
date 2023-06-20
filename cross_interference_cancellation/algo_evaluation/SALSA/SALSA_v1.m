function [a_r,a_i] = SALSA_v1(baseband_vec_signal, lambda, mu,max_itr)

[m,~] = size(baseband_vec_signal);
Fr = dftmtx(m)/(sqrt(m));
window_length = 128;
overlapping_length =122;
stft_size = 128;

d_r = rand(m,1);
a_r = Fr*baseband_vec_signal;%zeros(512,1);%

g = ones(window_length,1);%hann(window_length,'periodic');
a_i = stft(baseband_vec_signal,window = g,OverlapLength=overlapping_length,FFTLength=stft_size); %zeros(128,65);%

[a,b] = size(a_i);
d_i = rand(size(a_i));
convergence_vec = zeros(max_itr,1);

for i = 1:max_itr
    a_r_previous = a_r;
    a_i_previous = a_i;
    v_r = (a_r + d_r).*max(0, 1 - (lambda/2*mu)./abs(a_r + d_r)) - d_r;
   % keyboard
    v_i = (a_i + d_i).*max(0, 1 - ((1-lambda)/(2*mu))./abs(a_i + d_i)) - d_i;
    x = baseband_vec_signal - Fr'*a_r - istft(a_i,window = g,OverlapLength=overlapping_length,FFTLength=stft_size);
    d_r = (1/(2))*Fr*x;
    d_i = (1/(2))*stft(x,window = g,OverlapLength=overlapping_length,FFTLength=stft_size);
    a_r = d_r + v_r;
    a_i = d_i + v_i;

   

    relative_change_ar = norm(a_r - a_r_previous)/norm(a_r_previous);
    convergence_vec_ar(i) = relative_change_ar;
    
    relative_change_ai = norm(a_i - a_i_previous)/norm(a_i_previous);
    convergence_vec_ai(i) = relative_change_ai;
    disp([relative_change_ar relative_change_ai]);

end

% figure;
% plot(abs(a_r))
% title('a_r')
% 
% figure;
% imagesc(abs(a_i))
% title('a_i')
figure;
plot(convergence_vec_ar)
title('normalised relative change ar')
figure;
plot(convergence_vec_ai)
title('normalised relative change ai')



% fs = 4096;
% t = 0:1/fs:2-1/fs;
% x = chirp(t,250,1,500,'q');
% size(x)
% g = ones(128,1);
% y = stft(x,Window = g,OverlapLength=122,FFTLength=128);
% size(y)
% z = istft(y,Window = g,OverlapLength=122,FFTLength=128);
% size(z)
% z-y

