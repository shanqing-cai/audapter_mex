function TransShiftDemo_Monophthong_Ratio
%% Config
gender='female';
toPlay=0;

%% Path configuration
addpath('../BIN/release/');
addpath('commonmcode');
addpath('commonmcode/graph');

%% Set default parameters
p=getDefaultParams(gender);
MexIO('init',p);    % Initialize

%% Load the utterance data
utterFileName=['./da1_',gender,'.mat'];

load(utterFileName);  % gives data
fs=data.params.sr;
sigIn = data.signalIn;
    
TransShiftMex(6);   % Reset;

sigIn = resample(sigIn, 48000, data.params.sr);     
sigInCell = makecell(sigIn, 64);

data.params.bMelShift=0;    % Define the perturbation field in Hz, not in mel.
data.params.bRatioShift=1;  % Let the data in pertAmp be ratio of shift, instead of absolute amount of shift
data.params.LBk=0;          
data.params.LBb=0;
ddata.params.F2Min=1000; data.params.F2Max=2000;    % Redefine the lower and upper boundaries of the perturbation field
data.params.F1Min=0;  data.params.F1Max=1000;       % Redefine the left ant right boundaries of the field
data.params.pertF2=linspace(data.params.F2Min,data.params.F2Max,257);   % Redefine the independent variable of the perturbation field, in Hz
data.params.pertAmp=0.2*ones(size(data.params.pertAmp));            % Uniform 20% increase in F1

MexIO('init',data.params); % Set speaker-specific parameters

%% Run TransShiftMex over the input signal
procTimes=nan(size(sigInCell));
for n = 1 : length(sigInCell)
    tic;
    TransShiftMex(5,sigInCell{n});
    procTimes(n)=toc;
end

% figure;
% plot(procTimes*1e3);
% set(gca,'YLim',[0,2]);

[sig, b] = TransShiftMex(4);

if (toPlay)
    wavplay(2*resample(sig(:, 1),fs*4,fs), fs*4);
    wavplay(2*resample(sig(:, 2),fs*4,fs), fs*4);      
end

%% Spectrogram with estimated formants overlaid
fl=data.params.frameLen;
taxis2=0:(1/fs*fl):(1/fs*fl)*(size(b,1)-1);
% get the formant plot bounds
[i1,i2,f1,f2,iv1,iv2]=getFmtPlotBounds(b(:,5),b(:,6));

figure('Position',[200,200,800,400]);
[s,f,t]=spectrogram(sig(:,1),64,48,1024,fs);  
[s2,f2,t2]=spectrogram(sig(:,2),64,48,1024,fs);
subplot(121);
imagesc(t,f,10*log10(abs(s)));  hold on;
axis xy;    
plot(taxis2,b(:,5:6),'k','LineWidth',2);
set(gca,'XLim',[taxis2(i1),taxis2(i2)]);
set(gca,'YLim',[0,3000]);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Original');

subplot(122);
imagesc(t2,f2,10*log10(abs(s2))); hold on;
axis xy;
plot(taxis2,b(:,15:16),'g','LineWidth',2);
set(gca,'XLim',[taxis2(i1),taxis2(i2)]);
set(gca,'YLim',[0,3000]);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Shifted');
%% MATLAB-Based detection of transition

[k1,k2]=detectTriphTrans(b(:,5),b(:,6),iv1,iv2);

%%
figure('Position',[200,200,700,320]);
subplot(1,2,1);
plot(b(:,5:6), 'k','LineWidth',2);    hold on;
plot(b(:,15:16), 'g','LineWidth',2);
if (i2>i1)
    set(gca, 'XLim', [i1,i2]);
end
xlabel('Frame #');
ylabel('Formants (Hz)');
xs=get(gca,'XLim'); ys=get(gca,'YLim');
plot([k1,k1],[ys(1),ys(2)],'--','Color',[0.125,0.125,0.125]);
plot([k2,k2],[ys(1),ys(2)],'--','Color',[0.125,0.125,0.125]);    

subplot(1,2,2);
if (k2>k1)
    plot(b(k1:k2,5),b(k1:k2,6),'k','LineWidth',2);  hold on;
    sf1=b(k1:k2,15);
    sf2=b(k1:k2,16);
    idxnz=find(sf1~=0);
    sf1=sf1(idxnz);
    sf2=sf2(idxnz);
    plot(sf1,sf2,'g','LineWidth',2);
end
legend({'Original F1-F2 traj.', 'Shifted F1-F2 traj.'});
set(gca,'XLim',[0,1200]);   set(gca,'YLim',[500,3000]);
xlabel('F1 (Hz)');
ylabel('F2 (Hz)');
return