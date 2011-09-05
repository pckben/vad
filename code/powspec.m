function S = powspec(x,Nw,Nsh,nfft,fs)
    
    N = floor((length(x)-Nw)/Nsh);
    W = floor(nfft/2+1);
    S = zeros(W,N);
    
%     window = hann(Nw);
    Hs1 = spectrum.mtm(2);
    
    for n = 1:N
        x_i = x((1:Nw) + (n-1)*Nsh);
%         x_i = x_i.*window;
%         X_i = fft(x_i,nfft);
%         S(:,n) = X_i(1:W);
        hpsd = psd(Hs1,x_i,'Fs',fs,'NFFT',nfft);
        S(:,n) = hpsd.Data;
%         X_i = fft(x_i.*window,nfft);
    end
end