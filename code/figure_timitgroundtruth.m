% Evaluation

% database
TRAINDIR = '../SPEECHDATA/TRAIN/';
FILELIST = [TRAINDIR 'allwavs.txt'];
fid = fopen(FILELIST);
Files = textscan(fid,'%s'); Files = Files{1};
fclose(fid);

% constants
fs  = 16000;
Nw  = fs*.03;
Nsh = fs*.01;
nfft = 1024; %2^nextpow2(Nw*4);


fn = 35;

fprintf('(%3d%%) %s\n',round(fn/length(Files)*100),Files{fn});

% read audio file
x = read_audio([TRAINDIR Files{fn}],'sphere');
x = x/max(abs(x));

% read word file and extract manual labels
[path,name,~] = fileparts(Files{fn});
phn = read_phn(fullfile([TRAINDIR path],[name,'.PHN']));
lbl = extractLabel(length(x),phn);

% get energy measure on the clean speech
e = log10(energy(x,Nw,Nsh));
e_max = max(e);
e_min = min(e);
e_th = e_min + (e_max-e_min)*.2;    
e_decision = e > e_th;
e_decision = frame2signal(e_decision,Nw,Nsh);

% combine manual label and energy-based as referenced ground truth
ref = lbl;
ref(1:length(e_decision)) = ref(1:length(e_decision)) & e_decision;

figure(1);clf;
t = (1:length(x))/fs;
plot(t,x,'Color',[.85,.85,.85]);
hold on;
plot(t,-lbl*max(x)/4,'m');%,'LineWidth',1.2);
plot(t(1:length(e_decision)),-e_decision*max(x)/2,'b--');%,'LineWidth',1.2);
plot(t,ref*max(x)/2,'r','LineWidth',2);
axis tight;
colormap gray;
h=legend('speech','manual','energy','combined');
set(h,'Orientation','Horizontal');
hold off;
xlabel('time (s)');
ylabel('amplitude');
text(.2,.3,'(a) breath');
text(.77,.3,'(b)');
xlim([0 1.5]);