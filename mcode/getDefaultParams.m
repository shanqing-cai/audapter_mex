function p = getDefaultParams(sex,varargin)

switch sex
    case 'male'
        p.nLPC          = 13; 
        p.fn1           = 591;
        p.fn2           = 1314;
    case 'female'
        p.nLPC          = 11;	%SC-Mod(2008/02/08) Used to be 9
        p.fn1           = 675;
        p.fn2           = 1392;
    otherwise,
        error('specify sex (male / female');
end

p.aFact         = 1;
p.bFact         = 0.8;
p.gFact         = 1;
% p.nRampReps     = nRampReps;    % This is necessary because for the triphthong perturbation, both 
% phi and amp need to be updated nonlinearly for each rep of the ramp phase

p.downfact      = 4;  
p.dScale        = 1;
p.preempFact    = 0.98;% preemp factor 

p.sr            = 48000/p.downfact;

% framingxx
p.nWin          = 1;% 1 2 4  8 16 32 64 (max=p.framLen) Number of windows per frame  !!
p.frameLen      = 64/p.downfact;% must be a valid DMA Buffer size (64 128 256 ..)
p.nDelay        = 7;% the total process delay is: p.frameLen*p.nDelay/p.sr
p.frameShift    = p.frameLen/p.nWin;% 
p.bufLen        = (2*p.nDelay-1)*p.frameLen;
p.anaLen        = p.frameShift+2*(p.nDelay-1)*p.frameLen;
p.avgLen      = 8;    %ceil(p.sr/(f0*p.frameShift));

p.bCepsLift     = 1;

p.minVowelLen   = 60;

if (isequal(sex,'male'))
    p.cepsWinWidth  = 50;
elseif (isequal(sex,'female'))
    p.cepsWinWidth  = 30;    
end

% fmts
p.nFmts         = 2;%

% formant tracking 
p.nTracks       = 4;
p.bTrack        = 1;
p.bWeight       = 1; % weigthing (short time rms) of moving average formant estimate o
p.fmts_min      = [0;350;1200;2000;3000];
p.fmts_max      = [1500;3500;4500;5000;7000];
p.fn            = [500;1500;2500;3500;4500];
p.maxDelta      = 20;
p.trackIntroTime= 50;
p.trackBwStartW = 0.2;
p.trackMaStartW = 0;
p.trackFF       = 0.85;

%rms
%SC-Mod(2008/01/12)
p.rmsThresh     = 0.02*10^((getSPLTarg2('prod')-85)/20);
p.rmsRatioThresh= 1.3;% threshold for sibilant / vowel detection
p.rmsMeanPeak   = 6*p.rmsThresh;
p.rmsForgFact   = 0.95;% forgetting factor for rms computation

%transition detection
p.bDetect       = 1;
%p.fmtsForgFact  = 0.97;% formants forgetting factor [1 2 .. nFmts]
p.dFmtsForgFact = 0.93;% formants forgetting factor for derivate calculation
p.minDetected   = p.nWin;

% shifting
p.bShift        = 0;

p.bRatioShift   = 0;
p.bMelShift     = 1;

% p.shiftDirection=shifttype;

% if (isequal(p.shiftDirection,'inflate') | isequal(p.shiftDirection,'deflate') |...
%         isequal(p.shiftDirection,'f2_up') | isequal(p.shiftDirection,'f2_down') | ...
%         isequal(p.shiftDirection,'AccelDecel') |
%         isequal(p.shiftDirection,'DecelAccel'))
    p.gainAdapt     =0;
% else
%     p.bTriPertBaseF2=0;
%     p.gainAdapt     =1;
% end

if (isempty(findStringInCell(varargin,'ratio')))
    p.rMaxPert2FieldWidth=0.25;  % 0.2 or 0.25? 0.2 is more conservative. 0.25 is more optimistic. 
    p.pertRatio=p.rMaxPert2FieldWidth;
else
    p.rMaxPert2FieldWidth=varargin{findStringInCell(varargin,'ratio')+1};
    p.pertRatio=p.rMaxPert2FieldWidth;
end

p.fb=1; % Voice only;
p.lvNoise=0;

p.trialLen=2.5; %SC(2008/06/22)
p.rampLen=0.25; %SC(2008/06/22)
%% Triphthong-related variables: these are for the mel frequency space

p.F2Min=0;
p.F2Max=0;
%p.triF2MinMin=0;
%p.triF2MaxMax=0;
p.F1Min=0;
p.F1Max=0;
p.LBb=0;
p.LBk=0;
p.pertF2=zeros(1,257);
p.triPertL=zeros(1,257);
p.pertAmp=zeros(1,257);
p.pertPhi=zeros(1,257);
p.triTrajLen=0;
