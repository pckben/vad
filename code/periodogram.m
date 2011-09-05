function X = periodogram(x,Nw,Nsh,nfft)
    
    window = hann(Nw);
    T = floor((length(x)-Nw)/Nsh);
    W = floor(nfft/2);
    X = zeros(W,T);
    
%     if nargout>=2
%         A = zeros(W,N);
%     end

    for t = 1:T
        x_t = x((1:Nw) + (t-1)*Nsh);
        x_t = x_t.*window;
        X_t = fft(x_t,nfft);
        X(:,t) = X_t(1:W);
        
%         a = lpc(x_i,18);
%         h = freqz(1,a,W,16000);
%         h = 20*log10(abs(h));
%         S(:,n) = S(:,n) + max(h)-h;
        
%         if nargout>=2
%             a = lpc(x_i,18);
%             h = freqz(1,a,W,16000);
%             A(:,n) = abs(h);
%         end
    end
end