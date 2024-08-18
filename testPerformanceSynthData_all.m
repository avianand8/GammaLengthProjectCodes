% Algorithm Names MP,OMP,OMP-MAGE,HILBERT,CGT,FGLD,WAVELET
%algNames input added which is array of string containing algorithms to be used example 
%for MP,OMP ,OMP-MAGE and CGT us algnames =["MP","OMP","OMP-MAGE","CGT"] order of algorithm names does not matter.

function testPerformanceSynthData_all(thresholdFractionList,electrodeNum,numMeanBursts,dictionarySize,algNames)
% check for algorthims to plot
isMP =find(strcmp(algNames,"MP"))>0;
isOMP =find(strcmp(algNames,"OMP"))>0;
isOMPMAGE =find(strcmp(algNames,"OMP-MAGE"))>0;
isOMPGEAR =find(strcmp(algNames,"OMP-GEAR"))>0;
isHILBERT =find(strcmp(algNames,"HILBERT"))>0;
isCGT =find(strcmp(algNames,"CGT"))>0;
isFGLD =find(strcmp(algNames,"FGLD"))>0;
isWAVELET =find(strcmp(algNames,"WAVELET"))>0;
%%%%%%%%%%%%%%%

subjectName = 'alpa'; expDate = '120316'; protocolName = 'GRF_001';
gridType = 'Microelectrode'; folderSourceString = ''; cVal=100;

numThresholds = length(thresholdFractionList);

% BurstDataParameters
burstLenList = [ 0.05];
cvAmp=0.1;
stimulusPeriodS=[0.5 2];
baselinePeriodS=[-1.5 0];
gammaFreqRangeHz=[40 60];
numBurstLengths = length(burstLenList);
synthColorList = jet(numBurstLengths);
displayFlagBurst=0;

% CGT (Xing et al., 2014)
cgtGaborSDList = [12.5 25]/1000;
numCGTSDList = length(cgtGaborSDList);
displayFlagCGT=0;

% Wavelet
displayFlagWavelet=0;

% Feingold et al., 2015
filterOrderFeingold=4;
displayFlagFeingold=0;

% Hilbert
filterOrderHilbert=4;
displayFlagHilbert=0;

% MP
maxIteration=50; %#ok<*NASGU>
adaptiveDictionaryParam=0.9;
%dictionarySize=2500000;
displayFlagMP=0;
displayFlagOMP=0;
displayFlagOMPMAGE=0;
showMPResults=0;

% Initialize
if isCGT
    medianBurstLengthCGT = zeros(numBurstLengths,numCGTSDList,numThresholds);
    seBurstLengthCGT = zeros(numBurstLengths,numCGTSDList,numThresholds);
end
if isWAVELET
    medianBurstLengthWavelet = zeros(numBurstLengths,numThresholds);
    seBurstLengthWavelet = zeros(numBurstLengths,numThresholds);
end
if isFGLD
    medianBurstLengthFeingold = zeros(numBurstLengths,numThresholds);
    seBurstLengthFeingold = zeros(numBurstLengths,numThresholds);
end
if isHILBERT
    medianBurstLengthHilbert = zeros(numBurstLengths,numThresholds);
    seBurstLengthHilbert = zeros(numBurstLengths,numThresholds);
end
if isMP
    medianBurstLengthMP = zeros(numBurstLengths,numThresholds);
    seBurstLengthMP = zeros(numBurstLengths,numThresholds);
end
if isOMP
    medianBurstLengthOMP = zeros(numBurstLengths,numThresholds);
    seBurstLengthOMP = zeros(numBurstLengths,numThresholds);
end
if isOMPMAGE
    medianBurstLengthOMPMAGE = zeros(numBurstLengths,numThresholds);
    seBurstLengthOMPMAGE = zeros(numBurstLengths,numThresholds);
end
if isOMPGEAR
    medianBurstLengthOMPGEAR = zeros(numBurstLengths,numThresholds);
    seBurstLengthOMPGEAR = zeros(numBurstLengths,numThresholds);
end

for i=1:numBurstLengths
    burstLen=burstLenList(i);
    synthColorName = synthColorList(i,:);
    disp(['Burst Length: ' num2str(burstLen)]);
    
    [analogData,timeVals] = generateBurstData(subjectName,expDate,protocolName,gridType,folderSourceString,electrodeNum,cVal,burstLen,cvAmp,displayFlagBurst,synthColorName,stimulusPeriodS,gammaFreqRangeHz,numMeanBursts);
   
    diffPower = getChangeInPower(analogData,timeVals,stimulusPeriodS,baselinePeriodS,gammaFreqRangeHz);
    
    
    for ii=1:length(thresholdFractionList)
        thresholdFactor = diffPower*thresholdFractionList(ii);
        
        % Estimate burst length using CGT (Xing et al., 2012)
        if isCGT
            for j=1:numCGTSDList
                burstLengthCGT = getBurstLengthCGT(analogData,timeVals,thresholdFactor,displayFlagCGT,stimulusPeriodS,baselinePeriodS,gammaFreqRangeHz,cgtGaborSDList(j));
                [medianBurstLengthCGT(i,j,ii),seBurstLengthCGT(i,j,ii)] = getMedianAndSE(burstLengthCGT);
            end
        end
        
        if isWAVELET
             % Estimate burst length using Wavelet
            burstLengthWavelet = getBurstLengthWavelet(analogData,timeVals,thresholdFactor,displayFlagWavelet,stimulusPeriodS,baselinePeriodS,gammaFreqRangeHz);
            [medianBurstLengthWavelet(i,ii),seBurstLengthWavelet(i,ii)] = getMedianAndSE(burstLengthWavelet);
        end
        if isFGLD
            % Estimate burst length using Feingold et al., 2015
            burstLengthFeingold = getBurstLengthFeingold(analogData,timeVals,thresholdFactor,displayFlagFeingold,stimulusPeriodS,baselinePeriodS,gammaFreqRangeHz,filterOrderFeingold);
            [medianBurstLengthFeingold(i,ii),seBurstLengthFeingold(i,ii)] = getMedianAndSE(burstLengthFeingold);
        end
        if isHILBERT
            % Estimate burst length using Hilbert method
            burstLengthHilbert = getBurstLengthHilbert(analogData,timeVals,thresholdFactor,displayFlagHilbert,stimulusPeriodS,baselinePeriodS,gammaFreqRangeHz,filterOrderHilbert);
            [medianBurstLengthHilbert(i,ii),seBurstLengthHilbert(i,ii)] = getMedianAndSE(burstLengthHilbert);
        end
        
        % Estimate burst length using Stochastic MP
        if isMP
            if ii==1
                [burstLengthMP,~,~,gaborInfo,header] = getBurstLength_all(analogData,timeVals,sqrt(thresholdFactor),displayFlagMP,stimulusPeriodS,baselinePeriodS,gammaFreqRangeHz,maxIteration,adaptiveDictionaryParam,dictionarySize,[],[],"MP"); %#ok<*UNRCH>
            else
                burstLengthMP = getBurstLength_all(analogData,timeVals,sqrt(thresholdFactor),displayFlagMP,stimulusPeriodS,baselinePeriodS,gammaFreqRangeHz,maxIteration,adaptiveDictionaryParam,dictionarySize,gaborInfo,header,"MP");
            end
            [medianBurstLengthMP(i,ii),seBurstLengthMP(i,ii)] = getMedianAndSE(burstLengthMP);
        end
    

        if isOMP
            if ii==1
                [burstLengthOMP,~,~,gaborInfo1,header1] = getBurstLength_all(analogData,timeVals,sqrt(thresholdFactor),displayFlagOMP,stimulusPeriodS,baselinePeriodS,gammaFreqRangeHz,maxIteration,adaptiveDictionaryParam,dictionarySize,[],[],"OMP"); %#ok<*UNRCH>
            else
                burstLengthOMP = getBurstLength_all(analogData,timeVals,sqrt(thresholdFactor),displayFlagOMP,stimulusPeriodS,baselinePeriodS,gammaFreqRangeHz,maxIteration,adaptiveDictionaryParam,dictionarySize,gaborInfo1,header1,"OMP");
            end
            [medianBurstLengthOMP(i,ii),seBurstLengthOMP(i,ii)] = getMedianAndSE(burstLengthOMP);
        end
        if isOMPMAGE
            if ii==1
                [burstLengthOMPMAGE,~,~,gaborInfo2,header2] = getBurstLength_all(analogData,timeVals,sqrt(thresholdFactor),displayFlagOMPMAGE,stimulusPeriodS,baselinePeriodS,gammaFreqRangeHz,maxIteration,adaptiveDictionaryParam,dictionarySize,[],[],"OMP-MAGE"); %#ok<*UNRCH>
            else
                burstLengthOMPMAGE = getBurstLength_all(analogData,timeVals,sqrt(thresholdFactor),displayFlagOMPMAGE,stimulusPeriodS,baselinePeriodS,gammaFreqRangeHz,maxIteration,adaptiveDictionaryParam,dictionarySize,gaborInfo2,header2,"OMP-MAGE");
            end
            [medianBurstLengthOMPMAGE(i,ii),seBurstLengthOMPMAGE(i,ii)] = getMedianAndSE(burstLengthOMPMAGE);
        end
        if isOMPGEAR
            if ii==1
                [burstLengthOMPGEAR,~,~,gaborInfo22,header22] = getBurstLength_all(analogData,timeVals,sqrt(thresholdFactor),displayFlagOMPMAGE,stimulusPeriodS,baselinePeriodS,gammaFreqRangeHz,maxIteration,adaptiveDictionaryParam,dictionarySize,[],[],"OMP-GEAR"); %#ok<*UNRCH>
            else
                burstLengthOMPGEAR = getBurstLength_all(analogData,timeVals,sqrt(thresholdFactor),displayFlagOMPMAGE,stimulusPeriodS,baselinePeriodS,gammaFreqRangeHz,maxIteration,adaptiveDictionaryParam,dictionarySize,gaborInfo22,header22,"OMP-GEAR");
            end
            [medianBurstLengthOMPGEAR(i,ii),seBurstLengthOMPGEAR(i,ii)] = getMedianAndSE(burstLengthOMPGEAR);
        end
end

      
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%XingDataX = [0.025 0.05 0.075 0.1 0.125];
%XingDataY = [0.115 0.122 0.14 0.155 0.165];
caseFolder = ['PlotsGammaLength_' 'DSize_' num2str(dictionarySize)  num2str(numMeanBursts)];
if ~isfolder(caseFolder)
    mkdir(caseFolder);
end
%legend creation 
lgnd=[];
    if isCGT
        lgnd=[lgnd "CGT"];
    end
    if isWAVELET
        lgnd=[lgnd "WAVELET"];
    end
    if isHILBERT
        lgnd=[lgnd "HILBERT"];
    end
    if isMP
        lgnd=[lgnd "MP"];
    end
    if isOMP
        lgnd=[lgnd "OMP"];
    end
    if isOMPMAGE
        lgnd=[lgnd "OMP-MAGE"];
    end
    if isOMPGEAR
        lgnd=[lgnd "OMP-GEAR"];
    end
    lgnd=[lgnd "Injected Bursts"]
for ii=1:numThresholds
    figure(ii);
     hold on;
    
    colorNames = jet(numCGTSDList+7);
    if isCGT
        for i=1:numCGTSDList
            errorbar(burstLenList,squeeze(medianBurstLengthCGT(:,i,ii)),squeeze(seBurstLengthCGT(:,i,ii)),'color',colorNames(i,:),'LineWidth',2);
        end
    end
    if isWAVELET
        errorbar(burstLenList,medianBurstLengthWavelet(:,ii),seBurstLengthWavelet(:,ii),'color',colorNames(numCGTSDList+1,:),'LineWidth',2);
    end
    if isFGLD
        errorbar(burstLenList,medianBurstLengthFeingold(:,ii),seBurstLengthFeingold(:,ii),'color',colorNames(numCGTSDList+2,:),'LineWidth',2);
    end
    if isHILBERT
        errorbar(burstLenList,medianBurstLengthHilbert(:,ii),seBurstLengthHilbert(:,ii),'color',colorNames(numCGTSDList+3,:),'LineWidth',2);
    end
    if isMP
        errorbar(burstLenList,medianBurstLengthMP(:,ii),seBurstLengthMP(:,ii),'color',[0.8 0 0],'LineWidth',2);
    end
    if isOMP
        errorbar(burstLenList,medianBurstLengthOMP(:,ii),seBurstLengthOMP(:,ii),'color',[0 0 0.8],'LineWidth',2);
    end
     if isOMPMAGE
        errorbar(burstLenList,medianBurstLengthOMPMAGE(:,ii),seBurstLengthOMPMAGE(:,ii),'color',[0 0.8 0],'LineWidth',2);
    end
     if isOMPGEAR
        errorbar(burstLenList,medianBurstLengthOMPGEAR(:,ii),seBurstLengthOMPGEAR(:,ii),'color',[0 0.4 0],'LineWidth',2);
    end
    grid on;
    box on;
    xs=0:0.1:1;
    plot(xs,xs,'k','LineWidth',2);
    axis([0 1 0 1]);
    legend(lgnd);
    saveas(figure(ii),[caseFolder '/thF_' num2str(thresholdFractionList(ii)) '_.fig']);
end

    
end

function [medianAllBurstLength,seAllBurstLength,medianFirstBurstLength,seFirstBurstLength] = getMedianAndSE(lengthList)

allLengths=[]; firstLengths=[];
for j=1:length(lengthList)
    allLengths=cat(1,allLengths,lengthList{j}(:));
    if ~isempty(lengthList{j})
        firstLengths=cat(1,firstLengths,lengthList{j}(1));
    end
end
if ~isempty(allLengths)
    medianAllBurstLength = median(allLengths);
    seAllBurstLength = getSEMedian(allLengths);
    medianFirstBurstLength = median(firstLengths);
    seFirstBurstLength = getSEMedian(firstLengths);
else
    medianAllBurstLength = 0;
    seAllBurstLength = 0;
    medianFirstBurstLength = 0;
    seFirstBurstLength = 0;
end
end
