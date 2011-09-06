% implementation of VAD proposed in
% Ramirez et al. Efficient voice activity detection algorithms using
% long-term speech information. Speech Communication. 2004
%
classdef vadramirez04
    
    properties
        NFFT
        fs          % sampling rate
        Nw          % window size
        Nsh         % window shift size
        win         % the hamming window
        
        M           % number of neighbour frames each side used to calculate LTSE
        
        % noise updating:
        N           % average noise spectrum
        E0          % energy of the cleanest condition
        E1          % energy of the noisiest condition
        gamma0      % low energy threshold
        gamma1      % high energy threshold
        K           % number of neighbor frames (each side) for updating noise spectrum
        alpha       % noise spectrum updating rate
    end
    
    methods
        % constructor
        function obj = vadramirez04(NFFT,fs,Nw,Nsh,n,M)
            obj.NFFT = NFFT;
            obj.fs   = fs;
            obj.Nw   = Nw;
            obj.Nsh  = Nsh;
            obj.win  = hamming(Nw);
            
            % compute the average noise spectrum
            n = enframe(n,obj.win,Nsh);
            N = fft(n,NFFT,2);
            N = abs(N(:,1:end/2)');
            obj.N = mean(N,2);
            
            enrgy = sum(n.^2,2);
            obj.E0 = min(enrgy);
            obj.E1 = max(enrgy);

            obj.M = M;
            obj.K = 3;
            obj.alpha = 0.95;
        end
        
        % computes the Long-Term Spectral Divergence feature vector
        function [flag,d] = ltsd(obj,s)
            
            % compute the speech spectrum
            x = enframe(s,obj.win,obj.Nsh);
            X = fft(x,obj.NFFT,2);
            X = abs(X(:,1:end/2)');

            L = size(X,2);                   % number of frames
            d = zeros(L,1);                  % M-order LTSD
            flag = zeros(L,1);
            for l=1+obj.M:L-obj.M
                ltse = max(X(:,l-obj.M:l+obj.M),[],2); % max by column. Equation (1)
                d(l) = 10*log10(ltse'.^2 * (1./obj.N.^2)/obj.NFFT); % just following the paper. Equation (2)
                
%                 e = X(:,l)' .* X(:,l);      % frame energy
%                 flag = d(l)<obj.threshold(e); % vad decision
%                 if ~flag                    % update noise spectrum
%                     obj.updateN(X,l);
%                 end
                flag(l) = NaN; % not implemented
            end
        end
        
        % threshold as a function of noise energy
        function gamma = threshold(obj,e)
            % Equation (3):
            if e<=obj.E0,     gamma = obj.gamma0;
            elseif e>=obj.E1, gamma = obj.gamma1;
            else
                gamma = obj.gamma0 - ...
                        (obj.gamma0-obj.gamma1)/(1-obj.E1/obj.E0) + ...
                        (obj.gamma0-obj.gamma1)/(obj.E0-obj.E1)*e;
            end
        end
        
        % update noise spectrum during non-speech period at frame l-th
        function updateN(obj,X,l)
            N = sum(X(:,l-obj.K:l+obj.K),2)/(2*obj.K+1); % Equation (5)
            obj.N = obj.alpha*obj.N + (1-obj.alpha)*N;   % Equation (4)
        end
        
        % use training data to estimate optimal threshold value for gammas
        function train(obj,x,ref)
        end
        
        % perform VAD on a speech instance
        function flag = vad(obj,x)
        end
    end
    
end