% comparing feature distribution of various VAD's
%

if matlabpool('size')==0
    matlabpool open;
end

load('./speakers.mat');
nSpeakers = length(SPEAKERS);

% constants
fs  = 16000;
Nw  = fs*.02;
Nsh = fs*.01;
nfft = 1024; %2^nextpow2(Nw*4);

% function to calculate short-term energy of a signal x
energy = @(x) sum(enframe(x,Nw,Nsh).^2,2);

% initialize noise signal
[n,n_fs] = read_audio('../NOISEX-92/white.wav','wav');
n = resample(n,fs,n_fs); % NOISEX-92 data was sampled at 19980Hz
n = n/max(abs(n));
En = energy(n);      
Pn = sum(En)/length(En);

% there are certain types of noise that need a mask to extract energy.
% such as: babble, machinegun, ...
% En_sorted = sort(En);
% bottom20  = En_sorted(round(0.2*length(En)));
% top20     = En_sorted(round(0.8*length(En)));
% En_thresh = bottom20 + 0.1*(top20-bottom20);
% n_mask    = En > En_thresh;
% Pn        = sum(En(n_mask))/sum(n_mask); % noise power

SNR = 20:-5:-10;
nSNR = length(SNR);
n_Gain = zeros(nSpeakers,nSNR); % debug purpose, vector of noise gains

fx = cell(nSNR,nSpeakers); % feature value for signal
fn = cell(nSNR,nSpeakers); % feature value for noise

nSpeakersUsing = 100; % 100 speakers, 1000 utterances

for snrcnt=1:nSNR
    snr = SNR(snrcnt);
    fprintf('SNR=%d dB--------------------------------------------\n',snr);
    
    parfor s = 1:nSpeakersUsing
        data = load([SPEAKERS{s},'/data.mat']);
        x = data.x;
        ref = data.ref;
        
        % add noise to the desired SNR
        assert(length(n)>length(x));
        Ex = energy(x);
        Px = sum(Ex(ref==1))/sum(ref);  % speech power
        Pn1 = Px/(10^(snr/10));         % expected noise power: SNR = 10*log10(Px/Pn1);
        n_gain = sqrt(Pn1/Pn);          % noise gain
                                        n_Gain(s,snrcnt)=n_gain; %debug
        n1 = n*n_gain;                  % adjust the signals
        y  = x + n1(1:length(x));

        n_sample = n(end-fs*20:end);    % sample noise for calibration: last 20s
%         rmr = vadramirez04(nfft,fs,Nw,Nsh,n_sample,5);
%         [~,d] = rmr.ltsd(y);      

        pp = struct;
        pp.ti = Nsh/fs;                 % frame increment in seconds
%         pp.xn = min(snr,0);             % minimum prior SNR
        [~,z] = vadsohn(n_sample,fs,'bn');
        d = vadsohn(y,z,'bn');
        d = d(:,3);

        fx{snrcnt,s} = [fx{snrcnt,s};d(ref==1)];
        fn{snrcnt,s} = [fn{snrcnt,s};d(ref==0)];
        
%         subplot(311);
%         plot(x,'Color',[.85,.85,.85]);
%         hold on;
%         plot(r*max(x)/2,'r');
%         hold off;
%         axis tight;
%         subplot(312);
%         plot(y,'Color',[.8,.8,.8]);axis tight;
%         subplot(313);
%         plot(d);axis tight;
        
%         fprintf('(%3d%%) %s\n',round(s/nSpeakers*100),SPEAKERS{s});
        fprintf('.'); if mod(s,10)==0, fprintf('\n'); end
        
    end
    
    FX = []; FN = [];
    for s=1:nSpeakersUsing
        FX = [FX;fx{snrcnt,s}];
        FN = [FN;fn{snrcnt,s}];
    end
    
    figure(snrcnt);
    [cx,xx] = hist(FX,100);
    [cn,nn] = hist(FN,100);
    plot(xx,cx/trapz(xx,cx),'b',nn,cn/trapz(nn,cn),'r');
    title(['SNR=',num2str(snr),'dB']);

end
save('experiment1_sohn_white.mat','fx','fn','SNR','fs','Nw','Nsh','nfft');

d = zeros(nSNR,1);
for snrcnt=1:nSNR
    FX = []; FN = [];
    for s=1:nSpeakersUsing
        FX = [FX;fx{snrcnt,s}];
        FN = [FN;fn{snrcnt,s}];
    end
    d(snrcnt) = bhattacharyya(FX,FN);
end
plot(SNR,d);
