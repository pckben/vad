
% database
DATADIR = '../SPEECHDATA/';
SPEAKERS = {};
s = 0;

d0 = pwd;
cd(DATADIR);
d1 = dir;
checkdir = @(d) d.isdir && d.name(1)~='.';

for i=1:length(d1)
    if checkdir(d1(i))
        cd(d1(i).name);
        d2 = dir;
        for j=1:length(d2)
            if checkdir(d2(j))
                cd(d2(j).name);
                d3 = dir;
                for k=1:length(d3)
                    if checkdir(d3(k))
                        s = s+1;
                        SPEAKERS{s,1} = [pwd,'/',d3(k).name];
                    end
                end
                cd ..
            end
        end
        cd ..
    end
end
cd(d0);

% constants
fs  = 16000;
Nw  = fs*.02;
Nsh = fs*.01;
assert(Nw/Nsh==round(Nw/Nsh)); % only support Nw in multiples of Nsh
nfft = 1024; %2^nextpow2(Nw*4);

for s = 1:length(SPEAKERS)
    fprintf('(%3d%%) %s\n',round(s/length(SPEAKERS)*100),SPEAKERS{s});
    x   = [];
    ref = false(0);
    
    dir3 = dir([SPEAKERS{s} '/*.WAV']);
    
    for i=1:length(dir3)
        fn = [SPEAKERS{s} '/' dir3(i).name];
    
        % read audio file
        xi = read_audio(fn,'sphere');
        xi = xi./max(abs(xi));
        
        % get energy measure on the clean speech
        f      = enframe(xi,Nw,Nsh);
        nf     = size(f,1);
        ei     = sum(f.^2,2);
        ei_th  = max(f(1:3));
        ei_vad = ei > ei_th;        
%         ei_vad = frame2signal(ei_vad,Nw,Nsh);

        xi(Nw+(nf-1)*Nsh+1:end) = [];   % trimming xi to make it exact nf frames
        nx = length(xi);                % number of samples
        
        % read and extract manual labels
        [path,name,~] = fileparts(fn);
        phn = read_phn(fullfile(path,[name,'.PHN']));
        label = extractLabel(nx,phn);
        fl  = enframe(label,Nw,Nsh);  % frames of label vector
        li  = sum(fl,2)>=Nw/2;        % marked as 1 if no less than half are 1s
        
        % combine manual label and energy-based as referenced ground truth
        refi = li & ei_vad;
        
        % balancing the number of speech and non-speech frames such that
        % non-speech has 60% of the total 
        n_s = sum(refi==1);
        n_ns = sum(refi==0);
        n_ns_required = round(n_s/.4*.6);
        n_ns_adding   = n_ns_required-n_ns;
        
        len1 = floor(n_ns_adding/2);
        len2 = n_ns_adding-len1;
        sil1 = min(ei)*ones(len1*Nsh,1);  
        sil2 = min(ei)*ones(len2*Nsh,1);
        x   = [x;sil1;xi;sil2];
        
        if i>1, ref = [ref;false((Nw-Nsh)/Nsh,1)]; end
        
        ref = [ref;false(len1,1);refi;false(len2,1)];
        
    end
    
    assert((length(x)-Nw)/Nsh+1 == length(ref));
    save([SPEAKERS{s} '/data.mat'],'x','ref');
end

save('./speakers.mat','SPEAKERS');