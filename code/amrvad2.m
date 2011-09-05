% AMR VAD option 2
%
% @author: Pham Chau Khoa
% @date:   3 Sept 2011
%
% Implemented based on ETSI TS 126.094 v10.0.0 2011-04
%
function flag = amrvad2(s_m,params)

if nargin<2,
    params.lastsample = 0;
end

% Filter + downscaling
% speech = Pre_Process( y2, y1, x0, x1, new_speech );


% 4.3.1 Frequency domain conversion
% L = fs*0.02;                % subframe length
zeta = -0.8;                % pre-amphasis coefficient
% d = filter([1 -zeta],1,s);  % pre-emphasis filter, zeta can be .95-.98
                            % d(n) = s(n) + zeta*s(n-1)
d = s + zeta*[params.lastsample,s(1:end-1)];

M = 1024;                   % DFT length                           
G = fft(d,1024);

% 4.3.2 Channel energy estimator

% The channel table is defined below.  In this table, the
% lower and higher frequency coefficients for each of the 16
% channels are specified.  The table excludes the coefficients
% with numbers 0 (DC), 1, and 64 (Foldover frequency).  For
% these coefficients, the gain is always set at 1.0 (0 dB). */

ch_tbl = [
     2,  3
     4,  5
     6,  7
     8,  9
    10, 11
    12, 13
    14, 16
    17, 19
    20, 22
    23, 26
    27, 30
    31, 35
    36, 41
    42, 48
    49, 55
    56, 63]; % looks like can use polyfit n=4 to fit the curve??
  
% 4.3.3 Channel SNR estimator

% 4.3.4 Voice metric calculation

% 4.3.5 Frame SNR and long-term peak SNR calculation

% 4.3.6 Negative SNR sesitivity bias

% 4.3.7 VAD decision

% 4.3.8 spectral deviation estimator

% 4.3.9 Sinewave detection

% 4.3.10 Background noise update decision

% 4.3.11 Background noise estimate update

end
