function [burstLengthCGT,burstLengthFeingold,burstLengthMP,burstLengthOMP,burstLengthOMPMAGE,burstLengthOMPGEAR,burstLengthWT,burstLengthHilbert,diffPower,freqListCGT,freqListMP,freqListOMP,freqListOMPMAGE,freqListOMPGEAR,freqListWT,modListMP,modListOMP,modListOMPMAGE,modListOMPGEAR]=getBurstLengthRealData(subjectName,expDate,protocolName,electrodeNum,thresholdFraction,cVal,gammaFreqRangeHz,runMPAnalysisFlag,runOMPAnalysisFlag,runOMPMAGEAnalysisFlag,runOMPGEARAnalysisFlag,dictionarySizeMP,dictionarySizeOMP,dictionarySizeOMPMAGE,dictionarySizeOMPGEAR)

if ~exist('runMPAnalysisFlag','var');   runMPAnalysisFlag=1;            end
if ~exist('runMPAnalysisFlag','var');   runOMPAnalysisFlag=0;            end
if ~exist('runMPAnalysisFlag','var');   runOMPMAGEAnalysisFlag=0;            end

gridType = 'Microelectrode'; 
folderSourceString = '';
stimulusPeriodS=[0.5 2];
baselinePeriodS=[-1.5 0];

% CGT and Feingold and Hilbert
cgtGaborSD = 25/1000;
filterOrder=4;

% MP
maxIteration=50;
adaptiveDictionaryParam=0.9;
%dictionarySize=2500000;

% Get Data
folderName = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName,'segmentedData','LFP');
st = load(fullfile(folderName,['elec' num2str(electrodeNum) 'c' num2str(cVal) '.mat']));
analogData = st.analogData;
t  = load(fullfile(folderName,'lfpInfo.mat'));
timeVals = t.timeVals;

% Compute Threshold
diffPower = getChangeInPower(analogData,timeVals,stimulusPeriodS,baselinePeriodS,gammaFreqRangeHz);
thresholdFactor = diffPower*thresholdFraction;

% Get Burst lengths using CGT, Wavelet, Feingold and Hilbert
[burstLengthCGT,freqListCGT] = getBurstLengthCGT(analogData,timeVals,thresholdFactor,0,stimulusPeriodS,baselinePeriodS,gammaFreqRangeHz,cgtGaborSD);
[burstLengthWT,freqListWT] = getBurstLengthWavelet(analogData,timeVals,thresholdFactor,0,stimulusPeriodS,baselinePeriodS,gammaFreqRangeHz);
burstLengthFeingold = getBurstLengthFeingold(analogData,timeVals,thresholdFactor,0,stimulusPeriodS,baselinePeriodS,gammaFreqRangeHz,filterOrder);  
burstLengthHilbert = getBurstLengthHilbert(analogData,timeVals,thresholdFactor,0,stimulusPeriodS,baselinePeriodS,gammaFreqRangeHz,filterOrder);

if runMPAnalysisFlag
    % Save MP results
    folderNameMain = fullfile('data',subjectName,gridType,expDate,protocolName);
    tag = ['elec' num2str(electrodeNum) 'c' num2str(cVal)];
    
    % Output folder
    outputFolder = fullfile(folderNameMain,'mpAnalysis',tag);
    makeDirectory(outputFolder);
    
    % Save gaborInfo as a mat file
    gaborInfoFile = fullfile(outputFolder,['gaborInfo_G' num2str(gammaFreqRangeHz(1)) '-' num2str(gammaFreqRangeHz(2)) '_M' num2str(maxIteration) ...
        '_D' num2str(100*adaptiveDictionaryParam) '_R' num2str(dictionarySizeMP) '_S' num2str(1000*stimulusPeriodS(1)) '-' num2str(1000*stimulusPeriodS(2)) '.mat']);
    if exist(gaborInfoFile,'file')
        disp(['Opening saved file ' gaborInfoFile]);
        load(gaborInfoFile);
        [burstLengthMP,freqListMP,~,~,~,modListMP] = getBurstLength_all(analogData,timeVals,sqrt(thresholdFactor),0,stimulusPeriodS,baselinePeriodS,gammaFreqRangeHz,maxIteration,adaptiveDictionaryParam,dictionarySizeMP,gaborInfo,header,"MP");
    else
        [burstLengthMP,freqListMP,~,gaborInfo,header,modListMP] = getBurstLength_all(analogData,timeVals,sqrt(thresholdFactor),0,stimulusPeriodS,baselinePeriodS,gammaFreqRangeHz,maxIteration,adaptiveDictionaryParam,dictionarySizeMP,[],[],"MP");
        save(gaborInfoFile,'gaborInfo','header');
    end
else
    burstLengthMP=[];
    freqListMP=[];
    modListMP=[];
end
if runOMPAnalysisFlag
    % Save OMP results
    folderNameMain = fullfile('data',subjectName,gridType,expDate,protocolName);
    tag = ['elec' num2str(electrodeNum) 'c' num2str(cVal)];
    
    % Output folder
    outputFolder = fullfile(folderNameMain,'ompAnalysis',tag);
    makeDirectory(outputFolder);
    
    % Save gaborInfo as a mat file
    gaborInfoFile = fullfile(outputFolder,['gaborInfo_G' num2str(gammaFreqRangeHz(1)) '-' num2str(gammaFreqRangeHz(2)) '_M' num2str(maxIteration) ...
        '_D' num2str(100*adaptiveDictionaryParam) '_R' num2str(dictionarySizeOMP) '_S' num2str(1000*stimulusPeriodS(1)) '-' num2str(1000*stimulusPeriodS(2)) '.mat']);
    if exist(gaborInfoFile,'file')
        disp(['Opening saved file ' gaborInfoFile]);
        load(gaborInfoFile);
        [burstLengthOMP,freqListOMP,~,~,~,modListOMP] = getBurstLength_all(analogData,timeVals,sqrt(thresholdFactor),0,stimulusPeriodS,baselinePeriodS,gammaFreqRangeHz,maxIteration,adaptiveDictionaryParam,dictionarySizeOMP,gaborInfo,header,"OMP");
    else
        [burstLengthOMP,freqListOMP,~,gaborInfo,header,modListOMP] = getBurstLength_all(analogData,timeVals,sqrt(thresholdFactor),0,stimulusPeriodS,baselinePeriodS,gammaFreqRangeHz,maxIteration,adaptiveDictionaryParam,dictionarySizeOMP,[],[],"OMP");
        save(gaborInfoFile,'gaborInfo','header');
    end
else
    burstLengthOMP=[];
    freqListOMP=[];
    modListOMP=[];
end
if runOMPMAGEAnalysisFlag
    % Save OMP results
    folderNameMain = fullfile('data',subjectName,gridType,expDate,protocolName);
    tag = ['elec' num2str(electrodeNum) 'c' num2str(cVal)];
    
    % Output folder
    outputFolder = fullfile(folderNameMain,'ompmageAnalysis',tag);
    makeDirectory(outputFolder);
    
    % Save gaborInfo as a mat file
    gaborInfoFile = fullfile(outputFolder,['gaborInfo_G' num2str(gammaFreqRangeHz(1)) '-' num2str(gammaFreqRangeHz(2)) '_M' num2str(maxIteration) ...
        '_D' num2str(100*adaptiveDictionaryParam) '_R' num2str(dictionarySizeOMPMAGE) '_S' num2str(1000*stimulusPeriodS(1)) '-' num2str(1000*stimulusPeriodS(2)) '.mat']);
    if exist(gaborInfoFile,'file')
        disp(['Opening saved file ' gaborInfoFile]);
        load(gaborInfoFile);
        [burstLengthOMPMAGE,freqListOMPMAGE,~,~,~,modListOMPMAGE] = getBurstLength_all(analogData,timeVals,sqrt(thresholdFactor),0,stimulusPeriodS,baselinePeriodS,gammaFreqRangeHz,maxIteration,adaptiveDictionaryParam,dictionarySizeOMPMAGE,gaborInfo,header,"OMP-MAGE");
    else
        [burstLengthOMPMAGE,freqListOMPMAGE,~,gaborInfo,header,modListOMPMAGE] = getBurstLength_all(analogData,timeVals,sqrt(thresholdFactor),0,stimulusPeriodS,baselinePeriodS,gammaFreqRangeHz,maxIteration,adaptiveDictionaryParam,dictionarySizeOMPMAGE,[],[],"OMP-MAGE");
        save(gaborInfoFile,'gaborInfo','header');
    end
else
    burstLengthOMPMAGE=[];
    freqListOMPMAGE=[];
    modListOMPMAGE=[];
end
if runOMPGEARAnalysisFlag
    % Save OMP results
    folderNameMain = fullfile('data',subjectName,gridType,expDate,protocolName);
    tag = ['elec' num2str(electrodeNum) 'c' num2str(cVal)];
    
    % Output folder
    outputFolder = fullfile(folderNameMain,'ompgearAnalysis',tag);
    makeDirectory(outputFolder);
    
    % Save gaborInfo as a mat file
    gaborInfoFile = fullfile(outputFolder,['gaborInfo_G' num2str(gammaFreqRangeHz(1)) '-' num2str(gammaFreqRangeHz(2)) '_M' num2str(maxIteration) ...
        '_D' num2str(100*adaptiveDictionaryParam) '_R' num2str(dictionarySizeOMPMAGE) '_S' num2str(1000*stimulusPeriodS(1)) '-' num2str(1000*stimulusPeriodS(2)) '.mat']);
    if exist(gaborInfoFile,'file')
        disp(['Opening saved file ' gaborInfoFile]);
        load(gaborInfoFile);
        [burstLengthOMPGEAR,freqListOMPGEAR,~,~,~,modListOMPGEAR] = getBurstLength_all(analogData,timeVals,sqrt(thresholdFactor),0,stimulusPeriodS,baselinePeriodS,gammaFreqRangeHz,maxIteration,adaptiveDictionaryParam,dictionarySizeOMPGEAR,gaborInfo,header,"OMP-MAGE2");
    else
        [burstLengthOMPGEAR,freqListOMPGEAR,~,gaborInfo,header,modListOMPGEAR] = getBurstLength_all(analogData,timeVals,sqrt(thresholdFactor),0,stimulusPeriodS,baselinePeriodS,gammaFreqRangeHz,maxIteration,adaptiveDictionaryParam,dictionarySizeOMPGEAR,[],[],"OMP-MAGE2");
        save(gaborInfoFile,'gaborInfo','header');
    end
else
    burstLengthOMPGEAR=[];
    freqListOMPGEAR=[];
    modListOMPGEAR=[];
end
end