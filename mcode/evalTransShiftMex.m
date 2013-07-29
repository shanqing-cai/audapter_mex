function evalTransShiftMex(varargin)
%% Config
tvDir='../testVowels';
sex='female';

%% Vowel info
vowelName=  {'IY','IH','EH','AE','AH','AA','AO','UW','UH','ER','AY','EY','AW','OW'};

vowelF1(1,:)=[288, 446, 650, 856, 795, 880, 570, 311, 480, 490, 1000,630, 1000,570];
vowelF1(2,:)=[288, 446, 650, 856, 795, 880, 570, 311, 480, 490, 500, 420, 420, 311];
vowelF2(1,:)=[2700,2500,2400,1900,1450,1250,800, 820, 1200,1600,1300,2000,1200,800];
vowelF2(2,:)=[2700,2500,2400,1900,1450,1250,800, 820, 1200,1600,2700,2700,820, 780];

%% Path configuration
addpath('./');
addpath('./commonmcode');
addpath('./commonmcode/graph');
addpath('../speechsyn');
addpath('../BIN/release/');

% Initialization
p=getDefaultParams(sex);
MexIO('init',p);

%% Run evaluation
d=dir(fullfile(tvDir,[sex(1),'_*.wav']));

nFiles=length(d);

err.F1=[];
err.F2=[];
err.F0Begin=[];
err.F0Mode=[];  % 0: const; -1: falling; 1: rising. 
err.bCepsLift=[];
err.F0s=[];
err.vowel={};
for bCepsLift=[0,1]
    for n=1:nFiles
        vowelInfo=load(fullfile(tvDir,strrep(d(n).name,'.wav','.mat')));
        
        idxVowel=find(vowelF1(1,:)==vowelInfo.f1(1) & vowelF1(2,:)==vowelInfo.f1(end) & ...
            vowelF2(1,:)==vowelInfo.f2(1) & vowelF2(2,:)==vowelInfo.f2(end));
        if ~isempty(idxVowel)
            err.vowel{length(err.vowel)+1}=vowelName(idxVowel);
        else
            err.vowel{length(err.vowel)+1}=[];
        end
        if vowelInfo.f0(end)==vowelInfo.f0(1)
            err.F0Mode=[err.F0Mode,0];
        elseif vowelInfo.f0(end)>vowelInfo.f0(1)
            err.F0Mode=[err.F0Mode,1];
        else
            err.F0Mode=[err.F0Mode,-1];
        end
        [f1err,f2err]=TransShiftDemo_CepsLift('fileName',fullfile(tvDir,d(n).name),'toPlot',0,'bCepsLift',bCepsLift,...
            'sex',sex);
        err.F1=[err.F1,f1err];
        err.F2=[err.F2,f2err];
        err.F0Begin=[err.F0Begin,vowelInfo.f0(1)];
        err.bCepsLift=[err.bCepsLift,bCepsLift];    
        
        if isempty(find(err.F0s==err.F0Begin(end)))
            err.F0s=[err.F0s,err.F0Begin(end)];
        end
        
        disp(['Vowel: bCepsLift=',num2str(bCepsLift),' - ',d(n).name]);
    end
end
err.F0s=sort(err.F0s);

% idx=findStringInCell(err.vowel,'EH');
% err.F1=err.F1(idx);
% err.F2=err.F2(idx);
% err.F0Begin=err.F0Begin(idx);
% err.F0Mode=err.F0Mode(idx);  % 0: const; -1: falling; 1: rising. 
% err.bCepsLift=err.bCepsLift(idx);
% err.vowel=err.vowel(idx);
%% Visualization
% load('err');

verr.cl0.f0const.f1=nan(1,length(err.F0s));
verr.cl0.f0const.f2=nan(1,length(err.F0s));
verr.cl0.f0change.f1=nan(1,length(err.F0s));
verr.cl0.f0change.f1=nan(1,length(err.F0s));
verr.cl1.f0const.f1=nan(1,length(err.F0s));
verr.cl1.f0const.f2=nan(1,length(err.F0s));
verr.cl1.f0change.f1=nan(1,length(err.F0s));
verr.cl1.f0change.f1=nan(1,length(err.F0s));

for n=1:length(err.F0s)
    verr.cl0.f0const.f1(n)=mean(err.F1(find(err.F0Begin==err.F0s(n) & err.bCepsLift==0 & err.F0Mode==0)));
    verr.cl0.f0const.f2(n)=mean(err.F2(find(err.F0Begin==err.F0s(n) & err.bCepsLift==0 & err.F0Mode==0)));
    verr.cl0.f0change.f1(n)=mean(err.F1(find(err.F0Begin==err.F0s(n) & err.bCepsLift==0 & err.F0Mode~=0)));
    verr.cl0.f0change.f2(n)=mean(err.F2(find(err.F0Begin==err.F0s(n) & err.bCepsLift==0 & err.F0Mode~=0)));
    verr.cl1.f0const.f1(n)=mean(err.F1(find(err.F0Begin==err.F0s(n) & err.bCepsLift==1 & err.F0Mode==0)));
    verr.cl1.f0const.f2(n)=mean(err.F2(find(err.F0Begin==err.F0s(n) & err.bCepsLift==1 & err.F0Mode==0)));
    verr.cl1.f0change.f1(n)=mean(err.F1(find(err.F0Begin==err.F0s(n) & err.bCepsLift==1 & err.F0Mode~=0)));
    verr.cl1.f0change.f2(n)=mean(err.F2(find(err.F0Begin==err.F0s(n) & err.bCepsLift==1 & err.F0Mode~=0)));     
end

figure('Position',[300,200,800,350]);
subplot(121);
plot(err.F0s,verr.cl0.f0const.f1,'r.-','MarkerSize',18,'LineWidth',1); hold on;
plot(err.F0s,verr.cl0.f0change.f1,'rs-','LineWidth',1);
plot(err.F0s,verr.cl1.f0const.f1,'b.-','MarkerSize',18,'LineWidth',1);
plot(err.F0s,verr.cl1.f0change.f1,'bs-','LineWidth',1);
legend({'No liftering, F0 constant','No liftering, F0 change','Liftering, F0 constant','Liftering, F0 change'},...
    'Location','Northwest','FontSize',11);
set(gca,'FontSize',12);
xlabel('F0 at onset (Hz)');
ylabel('Mean RMS fraction error');
% set(gca,'YLim',[0,0.14]);
title('F1');
subplot(122);
plot(err.F0s,verr.cl0.f0const.f2,'r.-','MarkerSize',18,'LineWidth',1); hold on;
plot(err.F0s,verr.cl0.f0change.f2,'rs-','LineWidth',1);
plot(err.F0s,verr.cl1.f0const.f2,'b.-','MarkerSize',18,'LineWidth',1);
plot(err.F0s,verr.cl1.f0change.f2,'bs-','LineWidth',1);
% set(gca,'YLim',[0,0.1]);
set(gca,'FontSize',12);
xlabel('F0 at onset (Hz)');
ylabel('Mean RMS fraction error');
title('F2');
return