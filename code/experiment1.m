% comparing feature distribution of various VAD's
%

load('./speakers.mat');

% constants
fs  = 16000;
Nw  = fs*.02;
Nsh = fs*.01;
nfft = 1024; %2^nextpow2(Nw*4);
MIN_NS = 0.5; % minimum non-speech segment=50ms

fx = []; % feature value for signal
fn = []; % feature value for noise

% function to calculate short-term energy of a signal x
energy = @(x) sum(enframe(x,Nw,Nsh).^2,2);

n = read_audio('../NOISEX-92/volvo.wav','wav');
en = energy(n);      % noise energy


for s = 1:length(SPEAKERS)
    fprintf('(%3d%%) %s\n',round(s/length(SPEAKERS)*100),SPEAKERS{s});
    load([SPEAKERS{s},'/data.mat']);
    r = frame2signal(ref,Nw,Nsh);
    
    subplot(311);
    plot(x,'Color',[.85,.85,.85]);
    hold on;
    plot(r*max(x)/2,'r');
    hold off;
    axis tight;
    
    % add noise to the desired SNR
    assert(length(n)>length(x));
    ex = energy(x);
    [xw,n_gain] = addnoise1(x,nw(3000:end),10,e_vad,Nw,Nsh,ref);
    nw = nw * n_gain;
    subplot(312);
    plot(xw,'Color',[.8,.8,.8]);axis tight;
    
    rmr = vadramirez04(nfft,fs,Nw,Nsh,nw(3000:end),5);
    [~,d] = rmr.ltsd(xw);
    subplot(313);
    plot(d);axis tight;
    
    fx = [fx;d(ref==1)];
    fn = [fn;d(ref==0)];
    
    pause(.1);
    
    
end