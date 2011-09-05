% G729B VAD:
% Input: frame xi [Nw,1] in the time domain
% Output: speech/non-speech decision
function s = g729b(xi)
    N = length(xi);
    r = autocorr(xi);
    
    % features
    LSF = lsf(r);
    Ef  = 10*log10(r(1)/N);
    El  = lowenergy(xi);
    ZC  = zcr(xi);
end

% LSF: 
%   Line spectral frequencies, calculated from the first
%   11 terms of the autocorrelation using G729A procedures
% Input:
%   autocorrelation coefficients r [n,1]
% Output:
%   LSF coefficients l
function l = lsf(r)
    a = ac2poly(r);     % LPC coefficients from AutoCorr
    l = poly2lsf(a);    % LSF coefficients from LPC
end

% Full-band energy in the log10 scale
% Input:
%   e = 10 log10 (r_0/10)
% 
function e = fullenergy(r,N)
    e = 10*log10(r(1)/N);
end

% Low-band energy (<1kH)
function e = lowenergy(xi)
    fs = 16000; % TODO
    Wn = 1000 * 2/fs;
    [b,a] = butter(13,Wn,'low');
    xi_low = filtfilt(b,a,xi);
    e = sum(xi_low.^2);
end

% Zero-crossing rate
function z = zcr(xi)
    M = length(xi)-1;
    
    z = sum(abs(sgn(xi(2:end))-sgn(xi(1:end-1)))) / (2*M);
end

% the sign function
function s = sgn(x)
    s = x>=0;
end