function TransShiftDemo_Triphthong
%% Path configuration
addpath('../BIN/debug/');
toPlay=1;

%% Set default parameters
p=getDefaultParams('male');
MexIO('init',p);

%% Load the utterance data
utterFileName='./diao1_male.mat';

load(utterFileName);  % gives data
fs=data.params.sr;
sigIn = data.signalIn;
    
TransShiftMex(6);   % Reset;

sigIn = resample(sigIn, 48000, data.params.sr);     
sigInCell = makecell(sigIn, 64);

% Makeshift
data.params.avgLen=data.params.pitchLen; rmfield(data.params,'pitchLen');
rmfield(data.params,'bTriPertBaseF2');
data.params.F2Min=data.params.triF2Min; rmfield(data.params,'triF2Min');
data.params.F2Max=data.params.triF2Max; rmfield(data.params,'triF2Max');
data.params.F1Min=data.params.triF1Min; rmfield(data.params,'triF1Min');
data.params.F1Max=data.params.triF1Max; rmfield(data.params,'triF1Max');
data.params.pertF2=data.params.triPertF2; rmfield(data.params,'triPertF2');
data.params.pertAmp=data.params.triPertAmp; rmfield(data.params,'triPertAmp');
data.params.pertPhi=data.params.triPertPhi; rmfield(data.params,'triPertPhi');
data.params.LBk=data.params.triLBk; rmfield(data.params,'triLBk');
data.params.LBb=data.params.triLBb; rmfield(data.params,'triLBb');
rmfield(data.params,'triPertL');
rmfield(data.params,'lvNoise');
data.params.bWeight=data.params.bWeigth; rmfield(data.params,'bWeigth');
data.params.minVowelLen=60;
% ~Makeshift

MexIO('init', data.params); % Set speaker-specific parameters

%% Run TransShiftMex over the input signal
procTimes=nan(size(sigInCell));
for n = 1 : length(sigInCell)
    tic;
    TransShiftMex(5,sigInCell{n});
    procTimes(n)=toc;
end

figure;
plot(procTimes*1e3);
set(gca,'YLim',[0,2]);

[sig, b, c] = TransShiftMex(4);
dataOut=MexIO('getData');

if (toPlay)
    wavplay(2*resample(sig(:, 1),fs*4,fs), fs*4);
    wavplay(2*resample(sig(:, 2),fs*4,fs), fs*4);      
end

%% Spectrogram with estimated formants overlaid

fl=data.params.frameLen;
taxis2=0:(1/fs*fl):(1/fs*fl)*(size(b,1)-1);
% get the formant plot bounds
[i1,i2,f1,f2,iv1,iv2]=getFmtPlotBounds(b(:,5),b(:,6));

figure;    
[s,f,t]=spectrogram(sig(:,1),64,48,1024,fs);  
imagesc(t,f,10*log10(abs(s)));  hold on;
axis xy;    
plot(taxis2,b(:,5:7),'w');
set(gca,'XLim',[taxis2(i1),taxis2(i2)]);
set(gca,'YLim',[0,4000]);

%% MATLAB-Based detection of transition

[k1,k2]=detectTriphTrans(b(:,5),b(:,6),iv1,iv2);

%%
figure('Position',[200,200,700,320]);
subplot(1,2,1);
plot(b(:, 5:6), 'b','LineWidth',1);    hold on;
plot(data.fmts(:,1:2),'k','LineWidth',1);  hold on;
plot(b(:, 15:16), 'r','LineWidth',1);
plot(data.sfmts, 'g');
if (i2>i1)
    set(gca, 'XLim', [i1,i2]);
end
%     legend({'F1 original', 'F2 original', 'F1 shifted', 'F2 shifted'});
xlabel('Frame #');
ylabel('Formants (Hz)');
xs=get(gca,'XLim'); ys=get(gca,'YLim');
plot([data.transStart,data.transStart],[ys(1),ys(2)],'--','Color',[0.6,0.6,0.6]);
plot([data.transStop,data.transStop],[ys(1),ys(2)],'-','Color',[0.6,0.6,0.6]);    
plot([dataOut.transStart,dataOut.transStart],[ys(1),ys(2)],'--','Color',[0.125,0.125,0.125]);
plot([dataOut.transStop,dataOut.transStop],[ys(1),ys(2)],'-','Color',[0.125,0.125,0.125]);
plot([k1,k1],[ys(1),ys(2)],'--','Color',[0.125,0.125,0.125]);
plot([k2,k2],[ys(1),ys(2)],'--','Color',[0.125,0.125,0.125]);    

subplot(1,2,2);
%     plot(b(:,5),b(:,6),'b');    hold on;
%     plot(b(:,15),b(:,16),'r');
if (dataOut.transStop>dataOut.transStart)
    plot(data.fmts(dataOut.transStart:dataOut.transStop,1),...
        data.fmts(dataOut.transStart:dataOut.transStop,2),'b');

end
if (k2>k1)
    plot(dataOut.fmts(k1:k2,1),dataOut.fmts(k1:k2,2),'b','LineWidth',1);  hold on;
    plot(data.fmts(k1:k2,1),data.fmts(k1:k2,2),'k','LineWidth',1);  hold on;
    sf1=dataOut.sfmts(k1:k2,1);
    sf2=dataOut.sfmts(k1:k2,2);
    idxnz=find(sf1~=0);
    sf1=sf1(idxnz);
    sf2=sf2(idxnz);
    plot(sf1,sf2,'r','LineWidth',1);
end
set(gca,'XLim',[0,1200]);   set(gca,'YLim',[500,3000]);
xlabel('F1 (Hz)');
ylabel('F2 (Hz)');
return