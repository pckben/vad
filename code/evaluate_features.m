% Evaluation

% database
TRAINDIR = '../SPEECHDATA/TRAIN/';
SPEAKERS = {};
s = 0;
dir1 = dir(TRAINDIR);
for i=1:length(dir1)
    if dir1(i).isdir && dir1(i).name(1)~='.'
        dir2 = dir([TRAINDIR,dir1(i).name]);
        for j=1:length(dir2)
            if dir2(j).isdir && dir2(j).name(1)~='.'
                s = s+1;
                SPEAKERS{s} = [TRAINDIR,dir1(i).name,'/',dir2(j).name];
            end
        end
    end
end

% constants
fs  = 16000;
Nw  = fs*.02;
Nsh = fs*.01;
nfft = 1024; %2^nextpow2(Nw*4);
MIN_NS = 0.5; % minimum non-speech segment=50ms


for s = 1:length(SPEAKERS)
    fprintf('(%3d%%) %s\n',round(s/length(SPEAKERS)*100),SPEAKERS{s});
    x = [];
    l = [];
    e_vad = [];
    ref = [];
    
    dir3 = dir([SPEAKERS{s} '/*.WAV']);
    for i=1:length(dir3)
        fn = [SPEAKERS{s} '/' dir3(i).name];
    
        % read audio file
        xi = read_audio(fn,'sphere');
        xi = xi./max(abs(xi));

        % read and extract manual labels
        [path,name,~] = fileparts(fn);
        phn = read_phn(fullfile(path,[name,'.PHN']));
        li  = extractLabel(length(xi),phn);
        
        % get energy measure on the clean speech
        ei    = energy(xi,Nw,Nsh);
        ei_th  = max(xi(1:3));
        ei_vad = ei > ei_th;
        ei_vad = frame2signal(ei_vad,Nw,Nsh);
        ei_vad(end:length(xi)) = 0;

        % combine manual label and energy-based as referenced ground truth
        refi = zeros(length(xi),1);
        refi(1:length(ei_vad)) = li(1:length(ei_vad)) & ei_vad;
        
        % balancing the number of speech and non-speech frames such that
        % non-speech has 60% of the total 
        n_s = sum(refi==1);
        n_ns = sum(refi==0);
        n_ns_required = round(n_s/.4*.6);
        n_ns_adding   = n_ns_required-n_ns;
        
        l1  = floor(n_ns_adding/2);
        sil = 1e-20*ones(l1,1);      
        x   = [x;sil;xi;sil];
        ref = [ref;zeros(l1,1);refi;zeros(l1,1)] == 1;
        e_vad = [e_vad;zeros(l1,1);ei_vad;zeros(l1,1)];
        
        if i==1, 
            first = find(ref==1,1,'first');
            x = x(first:end);
            ref = ref(first:end);
            e_vad = e_vad(first:end);
        end
    end
    ref = vadhangover(ref,fs*MIN_NS); 
    
    subplot(211);
    plot(x,'Color',[.85,.85,.85]);
    hold on;
    plot(ref*max(x)/2,'r');%,'LineWidth',1.2);
    hold off;
    axis tight;
    
    nw = read_audio('../Noisefiles/pinknoise.wav','wav');
    nw = nw./max(abs(nw));
    [xw,n_gain] = addnoise1(x,nw(3000:end),10,e_vad,ref,Nw,Nsh);
    nw = nw * n_gain;
    subplot(212);
    plot(xw,'Color',[.8,.8,.8]);axis tight;
    
    
    g729.vad = G729(xw,fs,Nw,Nsh);
    g729.vad = vadhangover(g729.vad,fs*MIN_NS);
    g729.result = vadperformance(ref,g729.vad);
    subplot(212);hold on;plot(g729.vad*.5,'m');hold off;
    
    fprintf('G729:\nspeech:\t\trecall=%.2f\tprecision=%.2f\n',...
            g729.result.speech.recall,...
            g729.result.speech.precision);
%     fprintf('nonspeech:\trecall=%.2f\tprecision=%.2f\n',...
%             g729.result.nonspeech.recall,...
%             g729.result.nonspeech.precision);
        
%     pp.of = Nw/Nsh; % overlap factor
%     sohn.feature = vadsohn(xw,fs,'b');
    sohn.vad = vadsohn(xw,fs,'a');
    sohn.vad = vadhangover(sohn.vad,fs*MIN_NS);
    sohn.result = vadperformance(ref,sohn.vad);
    subplot(212);hold on;plot(sohn.vad*.25);axis tight;hold off;
    
    fprintf('Sohn:\nspeech:\t\trecall=%.2f\tprecision=%.2f\n',...
            sohn.result.speech.recall,...
            sohn.result.speech.precision);
%     fprintf('nonspeech:\trecall=%.2f\tprecision=%.2f\n',...
%             sohn.result.nonspeech.recall,...
%             sohn.result.nonspeech.precision);

    rmr = vadramirez04(nfft,fs,Nw,Nsh,nw(3000:end),5);
    ramirez04.feature = rmr.ltsd(xw);
        
    pause(.1);
    
    
end
