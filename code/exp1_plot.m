
load('./speakers.mat');
nSpeakers = length(SPEAKERS);
nSpeakersUsing = 100;

SNR = 20:-5:-10;
nSNR = length(SNR);

nFeatures = 8;
D = zeros(nSNR,nFeatures);


% initialize noise signal
noisetypes = {'babble';'buccaneer1';'buccaneer2';'destroyerengine';...
              'destroyerops';'f16';'factory1';'factory2';'hfchannel';...
              'leopard';'m109';'pink';'volvo';'machinegun';'white'};
          
for nt = 1:size(noisetypes,1)
    noisetype = noisetypes{nt};
    
    for snrcnt=1:nSNR
        snr = SNR(snrcnt);

        f = cell(nFeatures,2);
        load(sprintf('%s_%s_%ddB.mat','ramirez04',noisetype,snr));
        load(sprintf('%s_%s_%ddB.mat','sohn99',noisetype,snr));
        load(sprintf('%s_%s_%ddB.mat','ghosh11',noisetype,snr));
        load(sprintf('%s_%s_%ddB.mat','huang00',noisetype,snr));
        load(sprintf('%s_%s_%ddB.mat','energy',noisetype,snr));
        load(sprintf('%s_%s_%ddB.mat','spectralEntropy',noisetype,snr));
        load(sprintf('%s_%s_%ddB.mat','zcr',noisetype,snr));

        for s = 1:nSpeakersUsing
            data = load([SPEAKERS{s},'/data.mat']);
            x = data.x;
            ref = data.ref;

            f1 = ramirez04{s};
            f2 = sohn99{s};
            f3 = ghosh11{s};
            f4 = huang00{s};
            f5 = energyfeature{s};
            f6 = spectralentropy{s};
            f7 = zerocrossingrate{s};
            f8 = round(rand(size(ref)));

            f{1,1} = [f{1,1};f1(ref==1)]; f{1,2} = [f{1,2};f1(ref==0)];
            f{2,1} = [f{2,1};f2(ref==1)]; f{2,2} = [f{2,2};f2(ref==0)];
            f{3,1} = [f{3,1};f3(ref==1)]; f{3,2} = [f{3,2};f3(ref==0)];
            f{4,1} = [f{4,1};f4(ref==1)]; f{4,2} = [f{4,2};f4(ref==0)];
            f{5,1} = [f{5,1};f5(ref==1)]; f{5,2} = [f{5,2};f5(ref==0)];
            f{6,1} = [f{6,1};f6(ref==1)]; f{6,2} = [f{6,2};f6(ref==0)];
            f{7,1} = [f{7,1};f7(ref==1)]; f{7,2} = [f{7,2};f7(ref==0)];
            f{8,1} = [f{8,1};f8(ref==1)]; f{8,2} = [f{8,2};f8(ref==0)];

        end

    %     figure(snrcnt);
    %     [cx,xx] = hist(d1{1},100);
    %     [cn,nn] = hist(d1{2},100);
    %     plot(xx,cx/trapz(xx,cx),'b',nn,cn/trapz(nn,cn),'r');
    %     title(['SNR=',num2str(snr),'dB']);

        parfor i=1:nFeatures
            D(snrcnt,i) = bhattacharyya(f{i,1},f{i,2});
        end
    end

    figure;
    markers = '+o*.xsd^v><ph';
    colors  = gray(nFeatures+5);
    getMarker = @(i) markers(mod(i,numel(markers))+1);
    clf;hold on;
    for i=1:nFeatures
        plot(SNR,D(:,i),[getMarker(i),'-'],'Color',colors(i,:));
    end
    hold off;
    axis tight;
    xlabel('SNR (dB)');
    ylabel('Bhattacharyya Distance');
    h=legend('Ramirez 2004','Sohn 1999','Ghosh 2011','Huang 2000','Energy','Spectral Entropy','ZCR','Random');
    set(h,'Location','EastOutside');
    title(sprintf('Bhattacharyya Distance - %s noise',noisetype));
end