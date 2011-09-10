% comparing feature distribution of various VAD's
%

clear;clc;

if matlabpool('size')==0
    matlabpool open;
end
fprintf('Running on %d parallel labs.\n',matlabpool('size'));

load('./speakers.mat');
nSpeakers = length(SPEAKERS);

% constants
fs  = 16000;
Nw  = fs*.02;
Nsh = fs*.01;
nfft = 1024; %2^nextpow2(Nw*4);

% function to calculate short-term energy of a signal x
energy = @(x) sum(enframe(x,Nw,Nsh).^2,2);
% zero crossing rate:
sgn    = @(f) f>=0;
zcr    = @(f) sum(abs(sgn(f(:,2:end))-sgn(f(:,1:end-1))),2);

% initialize noise signal
noisetypes = {'babble';'buccaneer1';'buccaneer2';'destroyerengine';...
              'destroyerops';'f16';'factory1';'factory2';'hfchannel';...
              'leopard';'m109';'pink';'volvo';'machinegun';'white'};
          
for nt = 1:size(noisetypes,1)
    noisetype = noisetypes{nt};
    
    fprintf('NOISE TYPE=%s--------------------------------------------\n',noisetype);
        
    if any(strfind(['babble','white'],noisetype)) continue; end
    
    [n,n_fs] = read_audio(sprintf('../NOISEX-92/%s.wav',noisetype),'wav');
    n = resample(n,fs,n_fs); % NOISEX-92 data was sampled at 19980Hz
    n = n/max(abs(n));
    En = energy(n);      
    Pn = sum(En)/length(En);

    % there are certain types of noise that need a mask to extract energy.
    % such as: babble, machinegun, ...
    
    if any(strfind(['machinegun'],noisetype))
        En_sorted = sort(En);
        bottom20  = En_sorted(round(0.2*length(En)));
        top20     = En_sorted(round(0.8*length(En)));
        En_thresh = bottom20 + 0.1*(top20-bottom20);
        n_mask    = En > En_thresh;
        Pn        = sum(En(n_mask))/sum(n_mask); % noise power
    end
    
    SNR = 20:-5:-10;
    nSNR = length(SNR);

    nSpeakersUsing = 100; % 100 speakers, 1000 utterances

    RAMIREZ04 = 1;
    SOHN99    = 2;
    GHOSH11   = 3;
    HUANG00   = 4;
    ENERGY    = 5;
    SPENTROPY = 6;
    ZCR       = 7;
    nFeatures = 7;

    for snrcnt=1:nSNR
        snr = SNR(snrcnt);
        fprintf('SNR=%d dB--------------------------------------------\n',snr);

        ramirez04 = cell(nSpeakers,1); % feature value output of RAMIREZ 2004
        sohn99    = cell(nSpeakers,1); % feature value output of SOHN 1999
        ghosh11   = cell(nSpeakers,1); % feature value output of GHOSH 2011
        huang00   = cell(nSpeakers,1); % EnergyEntropy feature (HUANG&YANG 2000)
        energyfeature    = cell(nSpeakers,1); % log energy feature
        spectralentropy  = cell(nSpeakers,1); % entropy feature
        zerocrossingrate = cell(nSpeakers,1); % zero crossing rate

        parfor s=1:nSpeakersUsing
            data = load([SPEAKERS{s},'/data.mat']);
            x = data.x;
            ref = data.ref;

            % add noise to the desired SNR
            assert(length(n)>length(x));
            Ex = energy(x);
            Px = sum(Ex(ref==1))/sum(ref);  % speech power
            Pn1 = Px/(10^(snr/10));         % expected noise power: SNR = 10*log10(Px/Pn1);
            n_gain = sqrt(Pn1/Pn);          % noise gain
                                            %n_Gain(s,snrcnt)=n_gain; %debug
            n1 = n*n_gain;                  % adjust the signals
            y  = x + n1(1:length(x));

            n_sample = n(end-fs*20:end);    % sample noise for calibration: last 20s

            for fi=1:nFeatures
                switch(fi)
                    case RAMIREZ04
                        % RAMIREZ 2004
                        rmr = vadramirez04(nfft,fs,Nw,Nsh,n_sample,5);
                        [~,ramirez04{s}] = rmr.ltsd(y);      
                    case SOHN99
                        % SOHN 1999
                        pp = struct;
                        pp.ti = Nsh/fs;                 % frame increment in seconds
                        [~,z] = vadsohn(n_sample,fs,'bn');
                        temp = vadsohn(y,z,'bn');
                        sohn99{s} = temp(:,3);
                    case GHOSH11
                        % GHOSH 2011
                        ghosh = vadghosh11(nfft,Nw,Nsh,30,'b',20);
                        [~,ghosh11{s}] = ghosh.ltsv(y);
                    case HUANG00
                        % HUANG 2000
                        huang = vadhuang00(nfft,Nw,Nsh,n_sample);
                        huang00{s} = huang.eefeature(y);
                    case ENERGY
                        % Energy
                        energyfeature{s} = log10(energy(y));
                    case SPENTROPY
                        % Spectral Entropy
                        spectralentropy{s} = huang.spectral_entropy(y);
                    case ZCR
                        % Zero crossing rate
                        zerocrossingrate{s} = zcr(enframe(y,Nw,Nsh));
                end
            end
            fprintf('s=%3d,(%3d%%) %s\n',s,round(s/nSpeakersUsing*100),SPEAKERS{s});
        end


        save(sprintf('%s_%s_%ddB.mat','ramirez04',noisetype,snr),'ramirez04','SNR','fs','Nw','Nsh','nfft');
        save(sprintf('%s_%s_%ddB.mat','sohn99',noisetype,snr),'sohn99');
        save(sprintf('%s_%s_%ddB.mat','ghosh11',noisetype,snr),'ghosh11');
        save(sprintf('%s_%s_%ddB.mat','huang00',noisetype,snr),'huang00');
        save(sprintf('%s_%s_%ddB.mat','energy',noisetype,snr),'energyfeature');
        save(sprintf('%s_%s_%ddB.mat','spectralEntropy',noisetype,snr),'spectralentropy');
        save(sprintf('%s_%s_%ddB.mat','zcr',noisetype,snr),'zerocrossingrate');

    end
end

% d = zeros(nSNR,1);
% for snrcnt=1:nSNR
%     FX = []; FN = [];
%     for s=1:nSpeakersUsing
%         FX = [FX;fx{snrcnt,s}];
%         FN = [FN;fn{snrcnt,s}];
%     end
%     d(snrcnt) = bhattacharyya(FX,FN);
% end
% figure(nSNR+1);
% plot(SNR,d);
