%
% ADDNOISE
% @author: Chng Eng Siong
% return: y     noisy signal
%         n1    noise at the desired power
%
function [y,n1] = addnoise(x,n,SNR,fs)

P_n = sum(n.^2)/length(n);

frameDuration = 0.01; % 10msec per frame
frameLen = frameDuration*fs;
nFrame = floor(length(x)/frameLen);

xStart = 1; 
xEnd = frameLen;
for i=1:nFrame
    x_i = x(xStart:xEnd);
    frameEnergy(i) = sum(x_i.^2)/frameLen;     
    xStart = xEnd+1;
    xEnd = xEnd+frameLen;
end
sortFrameEnergy = sort(frameEnergy);
bottom20 = sortFrameEnergy(ceil(0.2*nFrame));
top20 = sortFrameEnergy(ceil(0.8*nFrame));
threshold2 = 0.1*(top20-bottom20)+bottom20;

sigMask = frameEnergy > threshold2;
% plot(sigMask*max(frameEnergy),'g');
P_x = sum(sigMask.*frameEnergy)/sum(sigMask);

%SNR = 10*log10(sigPower/noisePower);
%10^(SNR/10) = sigPower/noisePower;
%desired_noisePower = sigPower/(10^(SNR/10))
%

P_n1 = P_x/(10^(SNR/10));
gain = sqrt(P_n1/P_n);
n1 = gain * n;
% genNoisePower = sum(n1.^2)/length(n);

y = x + n1(1:length(x));

end