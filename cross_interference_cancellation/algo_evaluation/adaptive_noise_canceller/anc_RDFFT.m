function output = anc_RDFFT(RangeFFT,NFast,NSlow,flen)
% flen = 8;
% NFast = 256;
% NSlow = 64;
% threshold = 0.5; % If the reference input power is greater than the threshold, there must be interferences.
for i = 1 : NSlow
 
    priCH                     = RangeFFT(1:NFast/2, i);
    refCH                     = conj(flipud( RangeFFT(NFast/2+1 : NFast , i)));  
    
    priCH               = flip(priCH);
    % refCH               = circshift(refCH,1);
    refCH               = flip(refCH);
    
    RefPower            = sum(abs(refCH).^2)/(NFast/2);
%     if (RefPower > threshold)
        % With adaptive filter
        M   = flen;      % Length of adaptive filter
        mu  = 2/RefPower/50;   % Weight update constant
        [e, wo] = af_with_ref_cx(priCH, refCH, M, mu);
        e                  =flip(e);
        eout(:, i) = e;
%     else
%         % Without Apaptive filter
%         priCH               =flip(priCH);
%         eout(:, i) = priCH;
%     end
  
end
% RDAfterANC          = fftshift(fft(eout, NSlow, 2), 2);
% output = RDAfterANC;

output = eout;

end