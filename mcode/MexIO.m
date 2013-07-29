function varargout = MexIO(action,params,inFrame,varargin)
%

persistent p

toPrompt=0; % set to 1 when necessary during debugging

switch(action)
    case 'init',
        p=params;
  
        TransShiftMex(3,'srate',p.sr, toPrompt);
        TransShiftMex(3,'framelen',p.frameLen, toPrompt);
        
        TransShiftMex(3,'ndelay',p.nDelay, toPrompt);
        TransShiftMex(3,'nwin',p.nWin, toPrompt);
        TransShiftMex(3,'nlpc',p.nLPC, toPrompt);
        TransShiftMex(3,'nfmts',p.nFmts, toPrompt);
        TransShiftMex(3,'ntracks',p.nTracks, toPrompt);        
        TransShiftMex(3,'scale',p.dScale, toPrompt);
        TransShiftMex(3,'preemp',p.preempFact, toPrompt);
        TransShiftMex(3,'rmsthr',p.rmsThresh, toPrompt);
        TransShiftMex(3,'rmsratio',p.rmsRatioThresh, toPrompt);
        TransShiftMex(3,'rmsff',p.rmsForgFact, toPrompt);
        TransShiftMex(3,'dfmtsff',p.dFmtsForgFact, toPrompt);
        TransShiftMex(3,'bgainadapt',p.gainAdapt, toPrompt);
        TransShiftMex(3,'bshift',p.bShift, toPrompt);
        TransShiftMex(3,'btrack',p.bTrack, toPrompt);
        TransShiftMex(3,'bdetect',p.bDetect, toPrompt);      
        TransShiftMex(3,'avglen',p.avgLen, toPrompt);        
        TransShiftMex(3,'bweight',p.bWeight, toPrompt);    
        
        if (isfield(p,'minVowelLen'))
            TransShiftMex(3,'minvowellen',p.minVowelLen, toPrompt);
        end
        
        if (isfield(p,'bRatioShift'))
            TransShiftMex(3,'bratioshift',p.bRatioShift, toPrompt);
        end
        if (isfield(p,'bMelShift'))
            TransShiftMEx(3,'bmelshift',p.bMelShift, toPrompt);
        end
%% SC-Mod(2008/05/15) Cepstral lifting related
        if (isfield(p,'bCepsLift'))
            TransShiftMex(3,'bcepslift',p.bCepsLift, toPrompt);
        else
            TransShiftMex(3,'bcepslift',0, toPrompt);
        end
        if (isfield(p,'cepsWinWidth'))
            TransShiftMex(3,'cepswinwidth',p.cepsWinWidth, toPrompt);
        end        

%% SC-Mod(2008/04/04) Perturbatoin field related 
        if (isfield(p,'F2Min'))  % Mel
            TransShiftMex(3,'f2min',p.F2Min, toPrompt);
        end
        if (isfield(p,'F2Max'))  % Mel
            TransShiftMex(3,'f2max',p.F2Max, toPrompt);
        end
        if (isfield(p,'F1Min'))
            TransShiftMex(3,'f1min',p.F1Min, toPrompt);
        end
        if (isfield(p,'F1Max'))
            TransShiftMex(3,'f1max',p.F1Max, toPrompt);
        end
        if (isfield(p,'LBk'))
            TransShiftMex(3,'lbk',p.LBk, toPrompt);
        end
        if (isfield(p,'LBb'))
            TransShiftMex(3,'lbb',p.LBb, toPrompt);
        end
        if (isfield(p,'pertF2'))   % Mel, 257(=256+1) points
            TransShiftMex(3,'pertf2',p.pertF2, toPrompt);
        end
        if (isfield(p,'pertAmp'))   % Mel, 257 points
            TransShiftMex(3,'pertamp',p.pertAmp, toPrompt);
        end   
        if (isfield(p,'pertPhi'))   % Mel, 257 points
            TransShiftMex(3,'pertphi',p.pertPhi, toPrompt);
        end       
        
        if (isfield(p,'fb'))    % 2008/06/18
            TransShiftMex(3,'fb',p.fb, toPrompt);
        end       
        if (isfield(p,'trialLen'))  %SC(2008/06/22)
            TransShiftMex(3,'triallen',p.trialLen, toPrompt);
        else
            TransShiftMex(3,'triallen',2.5, toPrompt);
        end
        if (isfield(p,'rampLen'))  %SC(2008/06/22)
            TransShiftMex(3,'ramplen',p.rampLen, toPrompt);
        else
            TransShiftMex(3,'ramplen',0.05, toPrompt);
        end
        
        %SC(2008/07/16)
        if (isfield(p,'aFact'))
            TransShiftMex(3,'afact',p.aFact, toPrompt);
        else
            TransShiftMex(3,'afact',1, toPrompt);
        end
        if (isfield(p,'bFact'))
            TransShiftMex(3,'bfact',p.bFact, toPrompt);
        else
            TransShiftMex(3,'bfact',0.1, toPrompt);
        end
        if (isfield(p,'gFact'))
            TransShiftMex(3,'gfact',p.gFact, toPrompt);
        else
            TransShiftMex(3,'gfact',0.2, toPrompt);
        end
        
        if (isfield(p,'fn1'))
            TransShiftMex(3,'fn1',p.fn1, toPrompt);
        else
            TransShiftMex(3,'fn1',500, toPrompt);
        end
        if (isfield(p,'fn2'))
            TransShiftMex(3,'fn2',p.fn2, toPrompt);
        else
            TransShiftMex(3,'fn2',1500, toPrompt);
        end
        return;
%%            
    case 'process',

        TransShiftMex(5,inFrame);
        return;

    case 'getData',

        nout=nargout;

        [signalMat,dataMat]=TransShiftMex(4);       
        
        data=[];

        switch(nout)
            case 1,

%                 try
                data.signalIn       = signalMat(:,1);
                data.signalOut      = signalMat(:,2);

                data.intervals      = dataMat(:,1);
                data.rms            = dataMat(:,2:4);
                
                offS = 5;
                data.fmts           = dataMat(:,offS:offS+p.nTracks-1);
                data.rads           = dataMat(:,offS+p.nTracks:offS+2*p.nTracks-1);
                data.dfmts          = dataMat(:,offS+2*p.nTracks:offS+2*p.nTracks+1);
                data.sfmts          = dataMat(:,offS+2*p.nTracks+2:offS+2*p.nTracks+3);

                 offS = offS+2*p.nTracks+4;
                 data.ai             = dataMat(:,offS:offS+p.nLPC);
                data.params         = p;
                varargout(1)        = {data};

                return;

            case 2,
                varargout(1)        = {signalMat(:,1)};
                varargout(2)        = {signalMat(:,2)};
                return;

            case 3,

                varargout(1)        = {transdataMat(1:2,2)'};
                varargout(2)        = {transdataMat(1:2,3)'};
                varargout(3)        = {transdataMat(2,1)-transdataMat(1,1)};



                                
                return;

            otherwise,

        end


    otherwise,
        
        
    uiwait(errordlg(['No such action : ' action ],'!! Error !!'));


end
