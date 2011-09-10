classdef vadghosh11
    properties
        nfft
        win
        Nw
        Nsh
        R=30            % number of last frames used to compute LTSV
        mode            % mode flags. 
                        %   'b': bartlett-welch method, if not: periodogram
        M=20            % Bartlett-Welch spectrum smoothing window size
    end
    
    methods
        % constructor
        function obj = vadghosh11(nfft,Nw,Nsh,R,mode,M)
            obj.nfft = nfft;
            obj.Nw = Nw;
            obj.Nsh = Nsh;
            obj.win = hanning(Nw);
            obj.R   = R;
            obj.mode = mode;
            obj.M    = M;
        end
        
        % compute the Long-Term Signal Variability measure
        function [flag,f] = ltsv(obj,x)
            % compute the speech spectrum
            x = enframe(x,obj.win,obj.Nsh);
            X = fft(x,obj.nfft,2);
            X = X(:,1:end/2)';
            S = X.*conj(X);

            if any(obj.mode=='b')
                % Bartlett-Welch spectrum smoothing
                S1 = S;
                for p=2:obj.M
                    S1(p:end) = S1(p:end)+S(1:end-p+1);
                end
                S = S1;
            end
            
            L = size(X,2);                   % number of frames
            f = zeros(L,1);                  % R-order LTSV
            flag = zeros(L,1);
            R = obj.R;
            for m=R:L
                sumS = sum(S(:,m-R+1:m),2);  % sum of recent subband values
                scaledS = S(:,m-R+1:m)./repmat(sumS,1,R);
                entropy = -sum(scaledS.*log(scaledS),2);
                f(m) = var(entropy);
            end
            
            f(1:R-1) = min(f(R:end));
            f = log10(f);
        end
    end
    
end