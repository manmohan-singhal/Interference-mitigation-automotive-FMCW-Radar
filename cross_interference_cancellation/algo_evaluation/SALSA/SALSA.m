function [a_r,a_i] = SALSA(baseband_vec_signal, lambda, mu,max_itr)

[m,~] = size(baseband_vec_signal);
Fr = dftmtx(m)/(sqrt(m));

delta = 0.1;
a_r = Fr*baseband_vec_signal;
a_i = zeros(m,1);
d = zeros(m,1);

convergence_vec_ar = zeros(max_itr,1);
convergence_vec_ai = zeros(max_itr,1);

for i = 1:max_itr
    a_r_previous = a_r;
    a_i_previous = a_i;
   
    a_r = sign(a_r - delta*mu*(Fr'*a_r + a_i - d)).*max(abs(a_r - delta*mu*(Fr'*a_r + a_i - d)) - lambda*delta,0);

    a_i = sign(a_i - delta*mu*(Fr'*a_r + a_i - d)).*max(abs(a_i - delta*mu*(Fr'*a_r + a_i - d)) - (1-lambda)*delta,0);

     d = d - (Fr'*a_r + a_i - baseband_vec_signal);

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




