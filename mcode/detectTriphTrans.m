function [k1,k2]=detectTriphTrans(f1,f2,i1,i2,varargin)
% i1, i2: the start and end of the main vocal region
% k1, k2: the start and end of the triphthong (/iau/) transition. 

% Debug
if (nargin==0)
    load fmtdata1;
    f1=b(:,5);  f2=b(:,6);    
    i1=iv1;     i2=iv2;
    varargin{1}='plot';
end

%% Parameters
if (~isempty(findStringInCell(varargin,'iao')))
    vowel='iao';
elseif (~isempty(findStringInCell(varargin,'iou')))
    vowel='iou';
elseif (~isempty(findStringInCell(varargin,'ao')))
    vowel='ao';
elseif (~isempty(findStringInCell(varargin,'ia')))
    vowel='ia';
elseif (~isempty(findStringInCell(varargin,'uai')))
    vowel='uai';
elseif (~isempty(findStringInCell(varargin,'a')))
    vowel='a';
elseif (~isempty(findStringInCell(varargin,'i')))
    vowel='i';
elseif (~isempty(findStringInCell(varargin,'u')))
	vowel='u';    
else
    vowel='iao';
end
    
fitlen=35;
if (isequal(vowel,'iao') | isequal(vowel,'ia') | isequal(vowel,'iou'))
    f11=[200,800];     %[200,800]
    f21=[800,3000];    %[1000,3000]
elseif (isequal(vowel,'ao') | isequal(vowel,'a'))
    f11=[400,1600];
    f21=[500,2500];
elseif (isequal(vowel,'uai') | isequal(vowel,'i') | isequal(vowel,'u'))
    f11=[200,800];
    f21=[300,1600];
    idxrms=findStringInCell(varargin,'rms');
    rms=varargin{idxrms+1};
end

fmtMinDiff=0.5; % Hz/samp
dFmtEdge1=0.5;
dFmtEdge2=1;
exitCntThresh=30;
minTransLen=60;

%% 1. Smoothing
f1v=f1(i1:i2);
f2v=f2(i1:i2);
if (isequal(vowel,'uai') | isequal(vowel,'i') | isequal(vowel,'u'))
    rms=rms(i1:i2);
    maxrms=max(rms);
    minrms=min(rms);
end

f1v=mva(f1v,11);
f2v=mva(f2v,11);

if (size(f1v,1)>size(f1v,2))    f1v=f1v';   end
if (size(f2v,1)>size(f2v,2))    f2v=f2v';   end

%% 2. Detection
kstart=[];
kend=[];

bTrans=0;
hfit=floor(fitlen/2);
df1s=[];
df2s=[];
exitCnt=0;

for n=1:length(f1v)
    f1s=f1v(max(1,n-hfit):min(n+hfit,length(f1v))); % For calculating the tangent
    f2s=f2v(max(1,n-hfit):min(n+hfit,length(f1v))); 
    b1=regress(f1s',[ones(length(f1s),1),(1:length(f1s))']);
    b2=regress(f2s',[ones(length(f2s),1),(1:length(f2s))']);

    df1=b1(2);   % Unit: Hz/samp
    df2=b2(2);

    df1s=[df1s,df1];
    df2s=[df2s,df2];
    
    if (isequal(vowel,'iao') | isequal(vowel,'ia') | isequal(vowel,'iou'))
        if (~bTrans)
            if (f1v(n)>f11(1) & f1v(n)<f11(2) & f2v(n)>f21(1) & f2v(n)<f21(2) & ...
                    df1>dFmtEdge1 & df2<dFmtEdge1 & df1-df2>fmtMinDiff)
                bTrans=1;
                kstart=[kstart,n];
                exitCnt=0;
            end
        else
            if (df2>dFmtEdge2)
                exitCnt=exitCnt+1;
            else
                exitCnt=0;
            end
            if (exitCnt==exitCntThresh)
                bTrans=0;
                kend=[kend,n-exitCntThresh];
                exitCnt=0;
            end
        end
    elseif (isequal(vowel,'ao') | isequal(vowel,'a'))
        if (~bTrans)
            if (f1v(n)>f11(1) & f1v(n)<f11(2) & f2v(n)>f21(1) & f2v(n)<f21(2))
                bTrans=1;
                kstart=[kstart,n];
                exitCnt=0;
            end
        else
            if (df2>dFmtEdge2)
                exitCnt=exitCnt+1;
            else
                exitCnt=0;
            end
            if (exitCnt==exitCntThresh)
                bTrans=0;
                kend=[kend,n-exitCntThresh];
                exitCnt=0;
            end
        end
    elseif (isequal(vowel,'i') | isequal(vowel,'u'))
        if (~bTrans)
            if (rms(n)>0.3*(maxrms-minrms)+minrms)
                bTrans=1;
                kstart=[kstart,n];
                exitCnt=0;
            end
        else
            if (rms(n)<0.3*(maxrms-minrms)+minrms)
                bTrans=0;
                kend=[kend,n-floor(exitCntThresh)];
            end
        end        
    elseif (isequal(vowel,'uai'))
        if (~bTrans)
            if (f1v(n)>f11(1) & f1v(n)<f11(2) & f2v(n)>f21(1) & f2v(n)<f21(2) & ...
                    df1>dFmtEdge1 & df2>dFmtEdge1)
                bTrans=1;
                kstart=[kstart,n];
                exitCnt=0;
            end
        else
            if (df1<dFmtEdge1 & df2<dFmtEdge1)
                exitCnt=exitCnt+1;
            elseif (rms(n)<0.3*(maxrms-minrms)+minrms)
                exitCnt=exitCnt+10;
            else
                exitCnt=0;
            end
            if (exitCnt==floor(exitCntThresh))
                bTrans=0;
                kend=[kend,n-floor(exitCntThresh)];
                exitCnt=0;
            end
        end
    end
end

% Debug
% figure;
% plot(df1s,'b');   hold on;
% plot(df2s,'g');
% plot(df1s-df2s,'r');

%% Pick the right interval
kstart=kstart(1:length(kend));

if (isempty(kend) | isempty(kstart))
    k1=NaN;
    k2=NaN;
    return
end

durs=kend-kstart;
[jnk,k]=max(durs);
k1=kstart(k)+i1-1;
k2=kend(k)+i1-1;

%% Visualization
if (~isempty(findStringInCell(varargin,'plot')))
    figure;
    plot(f1v);  hold on;
    plot(f2v);  
    ys=get(gca,'YLim');
    for n=1:length(kstart)
        plot([kstart(n),kstart(n)],[ys(1),ys(2)],'b');   
    end
    for n=1:length(kend)
        plot([kend(n),kend(n)],[ys(1),ys(2)],'r');
    end
end

return