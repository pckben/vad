% implementation of (Huang and Yang, 2000)
% A Novel Approach to Robust Speech Endpoint Detection in Car Environments
classdef vadhuang00
    
    properties
        nfft
        Nw
        Nsh
        win
        E0      % estimated average noise energy
        H0      % estimated average noise entropy
    end
    
    methods
        function obj = vadhuang00(nfft,Nw,Nsh,n)
            obj.nfft = nfft;
            obj.Nw   = Nw;
            obj.Nsh  = Nsh;
            obj.win  = hamming(Nw);
            
            % estimate noise energy and entropy
            n        = enframe(n,obj.win,Nsh);
            obj.E0   = mean(sum(n.^2,2));
            obj.H0   = mean(obj.spectral_entropy(n));
        end
        
        function EE = eefeature(obj,x)
            H  = obj.spectral_entropy(x);   % entropy
            x  = enframe(x,obj.win,obj.Nsh);
            E  = sum(x.^2,2);               % energy
            EE = sqrt(1+abs((E-obj.E0).*(H-obj.H0))); % equation (7) (8)
        end
        
        function H = spectral_entropy(obj,x)
            x = enframe(x,obj.win,obj.Nsh); % enframed,             [nframes x Nw]
            X = fft(x,obj.nfft,2);
            X = X(:,1:end/2);              % frequency spectrum    [nframes x nbands]
            M = size(X,2);                  % number of bands
            S = X.*conj(X);                 % power spectrum        [nframes x nbands]
            P = S ./ repmat(sum(S,2),1,M);  % normalized spectrum,  eq.(3)
                                            % equals to subband pdfs [nframes x nbands] 
            H = sum(P.*log(P),2);           % spectral entropy      [nframes x 1]       eq.(6)
        end
    end
    
end