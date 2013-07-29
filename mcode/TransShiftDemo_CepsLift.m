function varargout=TransShiftDemo_CepsLift(varargin)
%% Config
fileName='../testVowels/female_EH_F0Rise';
if (~isempty(findStringInCell(varargin,'fileName')))
    fileName=varargin{findStringInCell(varargin,'fileName')+1};
end

toPlay=0;
toPlot=1;
bCepsLift=0;

if (~isempty(findStringInCell(varargin,'bCepsLift')))
    bCepsLift=varargin{findStringInCell(varargin,'bCepsLift')+1};
end

if (~isempty(findStringInCell(varargin,'toPlot')))
    toPlot=varargin{findStringInCell(varargin,'toPlot')+1};
end

sex='female';
if (~isempty(findStringInCell(varargin,'sex')))
    sex=varargin{findStringInCell(varargin,'sex')+1};
end

%% Path configuration
addpath('./');
addpath('../speechsyn');
addpath('../BIN/debug/');

%% Set default parameters
p=getDefaultParams(sex);
MexIO('init',p);

%% Load the utterance data
fs=p.sr;
    
TransShiftMex(6);   % Reset;

[sigIn,fs0]=wavread(fileName);
vowelInfo=load(strrep(fileName,'.wav','.mat'));
meanF0=mean(vowelInfo.f0);

sigIn = resample(sigIn, 48000, fs0);
sigInCell = makecell(sigIn, 64);

% p.nLPC=11;
p.cepsWinWidth=round(0.54/meanF0*12000);
p.bCepsLift=bCepsLift;

MexIO('init',p); % Set speaker-specific parameters

%% Run TransShiftMex over the input signal
procTimes=nan(size(sigInCell));
for n = 1 : length(sigInCell)
    tic;
    TransShiftMex(5,sigInCell{n});
    procTimes(n)=toc;
end

[sig, b] = TransShiftMex(4);

if (toPlay)
    wavplay(sig(:,1),p.sr);
    wavplay(sig(:,2),p.sr);
end

%% Calculate the absolute and relative errors
fl=p.frameLen;
taxis2=0:(1/fs*fl):(1/fs*fl)*(size(b,1)-1);
% get the formant plot bounds
% [i1,i2,f1,f2,iv1,iv2]=getFmtPlotBounds(b(:,5),b(:,6));

ti=vowelInfo.t(find(vowelInfo.t>=0.05 & vowelInfo.t<0.45));
f1real=vowelInfo.f1(find(vowelInfo.t>=0.05 & vowelInfo.t<0.45));
f2real=vowelInfo.f2(find(vowelInfo.t>=0.05 & vowelInfo.t<0.45));

f1est=interp1(taxis2,b(:,5),ti)';
f2est=interp1(taxis2,b(:,6),ti)';

f1error=rms((f1est-f1real)./f1real);
f2error=rms((f2est-f2real)./f2real);

if nargout==2
    varargout{1}=[f1error];
    varargout{2}=[f2error];
end

%% Spectrogram with estimated formants overlaid

if (toPlot)
    figure;    
    [s,f,t]=spectrogram(sig(:,1),64,56,1024,fs);  
    imagesc(t,f,20*log10(abs(s)));  hold on;
    axis xy;    

    plot(vowelInfo.t,vowelInfo.f1,'g','LineWidth',1);
    plot(vowelInfo.t,vowelInfo.f2,'g','LineWidth',2);
    plot(taxis2,b(:,5:6),'w');
    % set(gca,'XLim',[taxis2(i1),taxis2(i2)]);
    set(gca,'YLim',[0,3000]);

    title(['F1Err = ',sprintf('%.2f',f1error*100),'%',...
        '; F2Err = ',sprintf('%.2f',f2error*100),'%']);
end
return