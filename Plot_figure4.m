% This part requires the full dataset
% subjectName = 'alpa'; expDate = '120316'; protocolName = 'GRF_001'; plotLocation = [1 2 3];
% subjectName = 'kesari'; expDate = '230716'; protocolName = 'GRF_002'; plotLocation = [4 5 6];
% x=load([subjectName 'MicroelectrodeRFData']);
% electrodeList = x.highRMSElectrodes;

% Data is provided for these three electrodes
subjectName1 = 'alpa'; expDate1 = '120316'; protocolName1 = 'GRF_001'; plotLocation = [1 2 3];
subjectName2 = 'kesari'; expDate2 = '230716'; protocolName2 = 'GRF_002';
electrodeList_alpa  = [2 3 4 5 6 7 12 13 15 16 21 22 23 24 25 26 29 30 31 32 33 34 35 36 40 41 42 43 44 45 46 50 52 53 54 55 56 61 62 63 64 65 66 69 71 72 73 74 75 76 79 80 81 82 83 84 85 86 89 90 91 92 93 94 95]; %86-low, 83-medium, 29-high;
electrodeList_kesari = [22 34 40 42 43 49 50 51 52 53 54 55 59 60 61 62 63 65 66 69 70 71 72 73 74 75 76 79 80 81 82 83 84 85 89 90 92 93 94];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cVal = 100; gammaFreqRangeHz = [40 60]; thresholdFraction = 0.5;

numElectrodes_kesari = length(electrodeList_kesari);
numElectrodes_alpa = length(electrodeList_alpa);

diffPower_kesari = zeros(1,numElectrodes_alpa);
diffPower_kesari = zeros(1,numElectrodes_kesari);

methodNames{1} = 'HILBERT';
methodNames{2} = 'MP';
methodNames{3} = 'OMP';
methodNames{4} = 'OMP-MAGE';
methodNames{5} = 'OMP-GEAR';
m1 = "HILBERT";
m2 = "MP";
m3 = "OMP";
m4 = "OMP-MAGE";
m5 = "OMP-GEAR";
dictionarySizeMP = 2500000;
dictionarySizeOMP = 1500000;
dictionarySizeOMPMAGE=1500000;
dictionarySizeOMPGEAR=1500000;

runMPAnalysisFlag=1;
runOMPAnalysisFlag=1;
runOMPMAGEAnalysisFlag=1;
runOMPGEARAnalysisFlag=1;
numMethods = length(methodNames);

color_MP = [0.8 0 0];
color_OMP = [0 0 0.8];
color_OMPMAGE = [0 0.8 0];
color_OMPGEAR = [0 0.4 0];
color_HILBERT = [0.33 1 0.67];
col_map = [color_HILBERT;color_MP;color_OMP;color_OMPMAGE;color_OMPGEAR];

symbol_MP = 'x';
symbol_OMP = '^';
symbol_OMPMAGE = '>';
symbol_OMPGEAR = '<';
symbol_HILBERT = 'o';

allBurstLengthsAllElectrodes_alpa = cell(1,numMethods);
allBurstFreqsAllElectrodes_alpa = cell(1,numMethods);
allBurstLengthsAllElectrodes_kesari = cell(1,numMethods);
allBurstFreqsAllElectrodes_kesari = cell(1,numMethods);
numBursts_alpa = zeros(numMethods,numElectrodes_alpa);
numBursts_kesari = zeros(numMethods,numElectrodes_kesari);
allModsMP_alpa=[];
allNormalizedModsMP_alpa=[];
allModsOMP_alpa=[];
allNormalizedModsOMP_alpa=[];
allModsOMPMAGE_alpa=[];
allNormalizedModsOMPMAGE_alpa=[];
allModsOMPGEAR_alpa=[];
allNormalizedModsOMPGEAR_alpa=[];
allModsMP_kesari=[];
allNormalizedModsMP_kesari=[];
allModsOMP_kesari=[];
allNormalizedModsOMP_kesari=[];
allModsOMPMAGE_kesari=[];
allNormalizedModsOMPMAGE_kesari=[];
allModsOMPGEAR_kesari=[];
allNormalizedModsOMPGEAR_kesari=[];

% Filter lengths that are too long
lengthLimit = 0.8;
%%%% Alpa analysis %%%
for i=1:numElectrodes_alpa
    electrodeNum=electrodeList_alpa(i);
    disp(i);
    
    [~,~,burstLengths_alpa{2},burstLengths_alpa{3},burstLengths_alpa{4},burstLengths_alpa{5},~,burstLengths_alpa{1},diffPower_alpa(i),~,freqLists_alpa{2},freqLists_alpa{3},freqLists_alpa{4},freqLists_alpa{5},~,modListMP_alpa,modListOMP_alpa,modListOMPMAGE_alpa,modListOMPGEAR_alpa]=getBurstLengthRealData(subjectName1,expDate1,protocolName1,electrodeNum,thresholdFraction,cVal,gammaFreqRangeHz,runMPAnalysisFlag,runOMPAnalysisFlag,runOMPMAGEAnalysisFlag,runOMPGEARAnalysisFlag,dictionarySizeMP,dictionarySizeOMP,dictionarySizeOMPMAGE,dictionarySizeOMPGEAR);
%     burstLengths_alpa = burstLengths_alpa(~cellfun('isempty',burstLengths_alpa));
%     freqLists_alpa = freqLists_alpa(~cellfun('isempty',freqLists_alpa));
    allBurstLengthsSingleElectrodes_alpa=cell(1,numMethods);
    allBurstFreqsSingleElectrode_alpa=cell(1,numMethods);
    numTrials=length(burstLengths_alpa{1});
    numMethods=size(burstLengths_alpa,2);
    if runMPAnalysisFlag
        allModsMPSingleElectrode_alpa=[];
        for iii=1:numTrials
            allModsMPSingleElectrode_alpa=cat(1,allModsMPSingleElectrode_alpa,modListMP_alpa{iii});
        end
        
        allModsMP_alpa=cat(1,allModsMP_alpa,allModsMPSingleElectrode_alpa);
        allNormalizedModsMP_alpa=cat(1,allNormalizedModsMP_alpa,allModsMPSingleElectrode_alpa/max(allModsMPSingleElectrode_alpa));
    end
    allModsOMPSingleElectrode_alpa=[];
    for iii=1:numTrials
        allModsOMPSingleElectrode_alpa=cat(1,allModsOMPSingleElectrode_alpa,modListOMP_alpa{iii});
    end
    
    allModsOMP_alpa=cat(1,allModsOMP_alpa,allModsOMPSingleElectrode_alpa);
    allNormalizedModsOMP_alpa=cat(1,allNormalizedModsOMP_alpa,allModsOMPSingleElectrode_alpa/max(allModsOMPSingleElectrode_alpa));
    allModsOMPMAGESingleElectrode_alpa=[];
    for iii=1:numTrials
        allModsOMPMAGESingleElectrode_alpa=cat(1,allModsOMPMAGESingleElectrode_alpa,modListOMPMAGE_alpa{iii});
    end
    
    allModsOMPMAGE_alpa=cat(1,allModsOMPMAGE_alpa,allModsOMPMAGESingleElectrode_alpa);
    allNormalizedModsOMPMAGE_alpa=cat(1,allNormalizedModsOMPMAGE_alpa,allModsOMPMAGESingleElectrode_alpa/max(allModsOMPMAGESingleElectrode_alpa));
    allModsOMPGEARSingleElectrode_alpa=[];
    for iii=1:numTrials
        allModsOMPGEARSingleElectrode_alpa=cat(1,allModsOMPGEARSingleElectrode_alpa,modListOMPGEAR_alpa{iii});
    end
    
    allModsOMPGEAR_alpa=cat(1,allModsOMPGEAR_alpa,allModsOMPGEARSingleElectrode_alpa);
    allNormalizedModsOMPGEAR_alpa=cat(1,allNormalizedModsOMPGEAR_alpa,allModsOMPGEARSingleElectrode_alpa/max(allModsOMPGEARSingleElectrode_alpa));
    
    
    for ii=1:numMethods
        for iii=1:numTrials
            allBurstLengthsSingleElectrodes_alpa{ii}=cat(1,allBurstLengthsSingleElectrodes_alpa{ii},burstLengths_alpa{ii}{iii}(:));
            if ii ~=1
                allBurstFreqsSingleElectrode_alpa{ii}=cat(1,allBurstFreqsSingleElectrode_alpa{ii},freqLists_alpa{ii}{iii}(:));
            end
        end
        allBurstLengthsAllElectrodes_alpa{ii}=cat(1,allBurstLengthsAllElectrodes_alpa{ii},allBurstLengthsSingleElectrodes_alpa{ii});
        badIndices = find(allBurstLengthsSingleElectrodes_alpa{ii}>lengthLimit);

        if ii ~= 1
            allBurstFreqsAllElectrodes_alpa{ii}=cat(1,allBurstFreqsAllElectrodes_alpa{ii},allBurstFreqsSingleElectrode_alpa{ii});
        end
        burstLengthsTMP = allBurstLengthsSingleElectrodes_alpa{ii};
        burstLengthsTMP(badIndices)=[];
        numBursts_alpa(ii,i)=length(burstLengthsTMP)/numTrials;
        medianBurstLength_alpa(ii,i)=median(burstLengthsTMP);
    end
% allBurstFreqsAllElectrodes = allBurstFreqsAllElectrodes(~cellfun('isempty',allBurstFreqsAllElectrodes));
end

%%% kesari analysis %%
for i=1:numElectrodes_kesari
    electrodeNum=electrodeList_kesari(i);
    disp(i);
    
    [~,~,burstLengths_kesari{2},burstLengths_kesari{3},burstLengths_kesari{4},burstLengths_kesari{5},~,burstLengths_kesari{1},diffPower_kesari(i),~,freqLists_kesari{2},freqLists_kesari{3},freqLists_kesari{4},freqLists_kesari{5},~,modListMP_kesari,modListOMP_kesari,modListOMPMAGE_kesari,modListOMPGEAR_kesari]=getBurstLengthRealData(subjectName2,expDate2,protocolName2,electrodeNum,thresholdFraction,cVal,gammaFreqRangeHz,runMPAnalysisFlag,runOMPAnalysisFlag,runOMPMAGEAnalysisFlag,runOMPGEARAnalysisFlag,dictionarySizeMP,dictionarySizeOMP,dictionarySizeOMPMAGE,dictionarySizeOMPGEAR);
%     burstLengths_alpa = burstLengths_alpa(~cellfun('isempty',burstLengths_alpa));
%     freqLists_alpa = freqLists_alpa(~cellfun('isempty',freqLists_alpa));
    allBurstLengthsSingleElectrodes_kesari=cell(1,numMethods);
    allBurstFreqsSingleElectrode_kesari=cell(1,numMethods);
    numTrials=length(burstLengths_kesari{1});
    numMethods=size(burstLengths_kesari,2);
    if runMPAnalysisFlag
        allModsMPSingleElectrode_kesari=[];
        for iii=1:numTrials
            allModsMPSingleElectrode_kesari=cat(1,allModsMPSingleElectrode_kesari,modListMP_kesari{iii});
        end
        
        allModsMP_kesari=cat(1,allModsMP_kesari,allModsMPSingleElectrode_kesari);
        allNormalizedModsMP_kesari=cat(1,allNormalizedModsMP_kesari,allModsMPSingleElectrode_kesari/max(allModsMPSingleElectrode_kesari));
    end
    allModsOMPSingleElectrode_kesari=[];
    for iii=1:numTrials
        allModsOMPSingleElectrode_kesari=cat(1,allModsOMPSingleElectrode_kesari,modListOMP_kesari{iii});
    end
    
    allModsOMP_kesari=cat(1,allModsOMP_kesari,allModsOMPSingleElectrode_kesari);
    allNormalizedModsOMP_kesari=cat(1,allNormalizedModsOMP_kesari,allModsOMPSingleElectrode_kesari/max(allModsOMPSingleElectrode_kesari));
    allModsOMPMAGESingleElectrode_kesari=[];
    for iii=1:numTrials
        allModsOMPMAGESingleElectrode_kesari=cat(1,allModsOMPMAGESingleElectrode_kesari,modListOMPMAGE_kesari{iii});
    end
    
    allModsOMPMAGE_kesari=cat(1,allModsOMPMAGE_kesari,allModsOMPMAGESingleElectrode_kesari);
    allNormalizedModsOMPMAGE_kesari=cat(1,allNormalizedModsOMPMAGE_kesari,allModsOMPMAGESingleElectrode_kesari/max(allModsOMPMAGESingleElectrode_kesari));
    allModsOMPGEARSingleElectrode_kesari=[];
    for iii=1:numTrials
        allModsOMPGEARSingleElectrode_kesari=cat(1,allModsOMPGEARSingleElectrode_kesari,modListOMPGEAR_kesari{iii});
    end
    
    allModsOMPGEAR_kesari=cat(1,allModsOMPGEAR_kesari,allModsOMPGEARSingleElectrode_kesari);
    allNormalizedModsOMPGEAR_kesari=cat(1,allNormalizedModsOMPGEAR_kesari,allModsOMPGEARSingleElectrode_kesari/max(allModsOMPGEARSingleElectrode_kesari));
    
    for ii=1:numMethods
        for iii=1:numTrials
            allBurstLengthsSingleElectrodes_kesari{ii}=cat(1,allBurstLengthsSingleElectrodes_kesari{ii},burstLengths_kesari{ii}{iii}(:));
            if ii ~=1
                allBurstFreqsSingleElectrode_kesari{ii}=cat(1,allBurstFreqsSingleElectrode_kesari{ii},freqLists_kesari{ii}{iii}(:));
            end
        end
        allBurstLengthsAllElectrodes_kesari{ii}=cat(1,allBurstLengthsAllElectrodes_kesari{ii},allBurstLengthsSingleElectrodes_kesari{ii});
        badIndices = find(allBurstLengthsSingleElectrodes_kesari{ii}>lengthLimit);

        if ii ~= 1
            allBurstFreqsAllElectrodes_kesari{ii}=cat(1,allBurstFreqsAllElectrodes_kesari{ii},allBurstFreqsSingleElectrode_kesari{ii});
        end
        burstLengthsTMP = allBurstLengthsSingleElectrodes_kesari{ii};
        burstLengthsTMP(badIndices)=[];
        numBursts_kesari(ii,i)=length(burstLengthsTMP)/numTrials;
        medianBurstLength_kesari(ii,i)=median(burstLengthsTMP);
    end
% allBurstFreqsAllElectrodes = allBurstFreqsAllElectrodes(~cellfun('isempty',allBurstFreqsAllElectrodes));
end

%%%% PLotting  Alpa  data %%%
t = tiledlayout(2,3);
ax1 = nexttile;
%%% plot on axes 1 and properties %%%%
hold on;
plot(diffPower_alpa,medianBurstLength_alpa(1,:),'color',color_HILBERT,'marker',symbol_HILBERT,'linestyle','none');
plot(diffPower_alpa,medianBurstLength_alpa(2,:),'color',color_MP,'marker',symbol_MP,'linestyle','none');
plot(diffPower_alpa,medianBurstLength_alpa(3,:),'color',color_OMP,'marker',symbol_OMP,'linestyle','none');
plot(diffPower_alpa,medianBurstLength_alpa(4,:),'color',color_OMPMAGE,'marker',symbol_OMPMAGE,'linestyle','none');
plot(diffPower_alpa,medianBurstLength_alpa(5,:),'color',color_OMPGEAR,'marker',symbol_OMPGEAR,'linestyle','none');
hold off;
set(ax1,'Fontsize',20)
set(get(ax1,'xlabel'),'string','Change in Power(dB)','Fontsize',20);
set(get(ax1,'ylabel'),'string','Burst Duration(s)','Fontsize',20);
xlim([4 40])
ylim([0 0.5])

%%% plot on ax2 and properties %%
histC = 0.025:0.025:2;
x_HILBERT_alpa = allBurstLengthsAllElectrodes_alpa{1}(allBurstLengthsAllElectrodes_alpa{1}<lengthLimit);
x_MP_alpa = allBurstLengthsAllElectrodes_alpa{2}(allBurstLengthsAllElectrodes_alpa{2}<lengthLimit);
x_OMP_alpa = allBurstLengthsAllElectrodes_alpa{3}(allBurstLengthsAllElectrodes_alpa{3}<lengthLimit);
x_OMPMAGE_alpa = allBurstLengthsAllElectrodes_alpa{4}(allBurstLengthsAllElectrodes_alpa{4}<lengthLimit);
x_OMPGEAR_alpa = allBurstLengthsAllElectrodes_alpa{5}(allBurstLengthsAllElectrodes_alpa{5}<lengthLimit);
ax2 = nexttile;
hold on;
histVals_HILBERT_alpa = hist(x_HILBERT_alpa,histC);
histVals_MP_alpa = hist(x_MP_alpa,histC);
histVals_OMP_alpa = hist(x_OMP_alpa,histC);
histVals_OMPMAGE_alpa = hist(x_OMPMAGE_alpa,histC);
histVals_OMPGEAR_alpa = hist(x_OMPGEAR_alpa,histC);
plot(histC,histVals_HILBERT_alpa/sum(histVals_HILBERT_alpa),'color',color_HILBERT,'LineWidth',2);
plot(histC,histVals_MP_alpa/sum(histVals_MP_alpa),'color',color_MP,'LineWidth',2);
plot(histC,histVals_OMP_alpa/sum(histVals_OMP_alpa),'color',color_OMP,'LineWidth',2);
plot(histC,histVals_OMPMAGE_alpa/sum(histVals_OMPMAGE_alpa),'color',color_OMPMAGE,'LineWidth',2);
plot(histC,histVals_OMPGEAR_alpa/sum(histVals_OMPGEAR_alpa),'color',color_OMPGEAR,'LineWidth',2);
xlim([0 0.82])
hold off
set(get(ax2,'xlabel'),'string','Burst Duration(s)','Fontsize',20);
set(get(ax2,'ylabel'),'string','Fraction of Bursts(s)','Fontsize',20);
legend('HILBERT','MP','OMP','OMP-MAGE','OMP-GEAR')
set(ax2,'Fontsize',20)


%%% plotting ax3 %%% 
ax3 = nexttile;
f_MP_alpa = round(allBurstFreqsAllElectrodes_alpa{2}(allBurstLengthsAllElectrodes_alpa{2}<lengthLimit),1);
f_MP_unique_alpa = unique(round(allBurstFreqsAllElectrodes_alpa{2}(allBurstLengthsAllElectrodes_alpa{2}<lengthLimit),1));
[medianLength_MP_alpa,seLengthMP_alpa] = getmlse(f_MP_alpa,f_MP_unique_alpa,x_MP_alpa);
hold on;
errorbar(f_MP_unique_alpa,medianLength_MP_alpa,seLengthMP_alpa,'Color',color_MP,'LineWidth',2);

f_OMP_alpa = round(allBurstFreqsAllElectrodes_alpa{3}(allBurstLengthsAllElectrodes_alpa{3}<lengthLimit),1);
f_OMP_unique_alpa = unique(round(allBurstFreqsAllElectrodes_alpa{3}(allBurstLengthsAllElectrodes_alpa{3}<lengthLimit),1));
[medianLength_OMP_alpa,seLengthOMP_alpa] = getmlse(f_OMP_alpa,f_OMP_unique_alpa,x_OMP_alpa);
errorbar(f_OMP_unique_alpa,medianLength_OMP_alpa,seLengthOMP_alpa,'Color',color_OMP,'LineWidth',2);
f_OMPMAGE_alpa = round(allBurstFreqsAllElectrodes_alpa{4}(allBurstLengthsAllElectrodes_alpa{4}<lengthLimit),1);
f_OMPMAGE_unique_alpa = unique(round(allBurstFreqsAllElectrodes_alpa{4}(allBurstLengthsAllElectrodes_alpa{4}<lengthLimit),1));
[medianLength_OMPMAGE_alpa,seLengthOMPMAGE_alpa] = getmlse(f_OMPMAGE_alpa,f_OMPMAGE_unique_alpa,x_OMPMAGE_alpa);
errorbar(f_OMPMAGE_unique_alpa,medianLength_OMPMAGE_alpa,seLengthOMPMAGE_alpa,'Color',color_OMPMAGE,'LineWidth',2);
f_OMPGEAR_alpa = round(allBurstFreqsAllElectrodes_alpa{5}(allBurstLengthsAllElectrodes_alpa{5}<lengthLimit),1);
f_OMPGEAR_unique_alpa = unique(round(allBurstFreqsAllElectrodes_alpa{5}(allBurstLengthsAllElectrodes_alpa{5}<lengthLimit),1));
[medianLength_OMPGEAR_alpa,seLengthOMPGEAR_alpa] = getmlse(f_OMPGEAR_alpa,f_OMPGEAR_unique_alpa,x_OMPGEAR_alpa);
errorbar(f_OMPGEAR_unique_alpa,medianLength_OMPGEAR_alpa,seLengthOMPGEAR_alpa,'Color',color_OMPGEAR,'LineWidth',2);
hold off
set(get(ax3,'xlabel'),'string','Burst Centre Freq(Hz)','Fontsize',20);
set(get(ax3,'ylabel'),'string','Burst Duration(s)','Fontsize',20);
set(ax3,'Fontsize',20)




%%%% Plotting kesari data %%%%

ax5 = nexttile;
%%% plot on axes 5 and properties %%%%
hold on;
plot(diffPower_kesari,medianBurstLength_kesari(1,:),'color',color_HILBERT,'marker',symbol_HILBERT,'linestyle','none');
plot(diffPower_kesari,medianBurstLength_kesari(2,:),'color',color_MP,'marker',symbol_MP,'linestyle','none');
plot(diffPower_kesari,medianBurstLength_kesari(3,:),'color',color_OMP,'marker',symbol_OMP,'linestyle','none');
plot(diffPower_kesari,medianBurstLength_kesari(4,:),'color',color_OMPMAGE,'marker',symbol_OMPMAGE,'linestyle','none');
plot(diffPower_kesari,medianBurstLength_kesari(5,:),'color',color_OMPGEAR,'marker',symbol_OMPGEAR,'linestyle','none');
hold off;
set(get(ax5,'xlabel'),'string','Change in Power(dB)','Fontsize',20);
set(get(ax5,'ylabel'),'string','Burst Duration(s)','Fontsize',20);
set(ax5,'Fontsize',20)
xlim([4 40])
ylim([0 0.5])

%%% plot on ax6 and properties %%
histC = 0.025:0.025:2;
x_HILBERT_kesari = allBurstLengthsAllElectrodes_kesari{1}(allBurstLengthsAllElectrodes_kesari{1}<lengthLimit);
x_MP_kesari = allBurstLengthsAllElectrodes_kesari{2}(allBurstLengthsAllElectrodes_kesari{2}<lengthLimit);
x_OMP_kesari = allBurstLengthsAllElectrodes_kesari{3}(allBurstLengthsAllElectrodes_kesari{3}<lengthLimit);
x_OMPMAGE_kesari = allBurstLengthsAllElectrodes_kesari{4}(allBurstLengthsAllElectrodes_kesari{4}<lengthLimit);
x_OMPGEAR_kesari = allBurstLengthsAllElectrodes_kesari{5}(allBurstLengthsAllElectrodes_kesari{5}<lengthLimit);
ax6 = nexttile;
hold on;
histVals_HILBERT_kesari = hist(x_HILBERT_kesari,histC);
histVals_MP_kesari = hist(x_MP_kesari,histC);
histVals_OMP_kesari = hist(x_OMP_kesari,histC);
histVals_OMPMAGE_kesari = hist(x_OMPMAGE_kesari,histC);
histVals_OMPGEAR_kesari = hist(x_OMPGEAR_kesari,histC);
plot(histC,histVals_HILBERT_kesari/sum(histVals_HILBERT_kesari),'color',color_HILBERT,'LineWidth',2);
plot(histC,histVals_MP_kesari/sum(histVals_MP_kesari),'color',color_MP,'LineWidth',2);
plot(histC,histVals_OMP_kesari/sum(histVals_OMP_kesari),'color',color_OMP,'LineWidth',2);
plot(histC,histVals_OMPMAGE_kesari/sum(histVals_OMPMAGE_kesari),'color',color_OMPMAGE,'LineWidth',2);
plot(histC,histVals_OMPGEAR_kesari/sum(histVals_OMPGEAR_kesari),'color',color_OMPGEAR,'LineWidth',2);
xlim([0 0.82])
hold off
set(get(ax6,'xlabel'),'string','Burst Duration(s)','Fontsize',20);
set(get(ax6,'ylabel'),'string','Fraction of Bursts(s)','Fontsize',20);
legend('HILBERT','MP','OMP','OMP-MAGE','OMP-GEAR')
set(ax6,'Fontsize',20)

%%% plotting ax7 %%% 
ax7 = nexttile;
f_MP_kesari = round(allBurstFreqsAllElectrodes_kesari{2}(allBurstLengthsAllElectrodes_kesari{2}<lengthLimit),1);
f_MP_unique_kesari = unique(round(allBurstFreqsAllElectrodes_kesari{2}(allBurstLengthsAllElectrodes_kesari{2}<lengthLimit),1));
[medianLength_MP_kesari,seLengthMP_kesari] = getmlse(f_MP_kesari,f_MP_unique_kesari,x_MP_kesari);
hold on;
errorbar(f_MP_unique_kesari,medianLength_MP_kesari,seLengthMP_kesari,'Color',color_MP,'LineWidth',2);

f_OMP_kesari = round(allBurstFreqsAllElectrodes_kesari{3}(allBurstLengthsAllElectrodes_kesari{3}<lengthLimit),1);
f_OMP_unique_kesari = unique(round(allBurstFreqsAllElectrodes_kesari{3}(allBurstLengthsAllElectrodes_kesari{3}<lengthLimit),1));
[medianLength_OMP_kesari,seLengthOMP_kesari] = getmlse(f_OMP_kesari,f_OMP_unique_kesari,x_OMP_kesari);
errorbar(f_OMP_unique_kesari,medianLength_OMP_kesari,seLengthOMP_kesari,'Color',color_OMP,'LineWidth',2);
f_OMPMAGE_kesari = round(allBurstFreqsAllElectrodes_kesari{4}(allBurstLengthsAllElectrodes_kesari{4}<lengthLimit),1);
f_OMPMAGE_unique_kesari = unique(round(allBurstFreqsAllElectrodes_kesari{4}(allBurstLengthsAllElectrodes_kesari{4}<lengthLimit),1));
[medianLength_OMPMAGE_kesari,seLengthOMPMAGE_kesari] = getmlse(f_OMPMAGE_kesari,f_OMPMAGE_unique_kesari,x_OMPMAGE_kesari);
errorbar(f_OMPMAGE_unique_kesari,medianLength_OMPMAGE_kesari,seLengthOMPMAGE_kesari,'Color',color_OMPMAGE,'LineWidth',2);
f_OMPGEAR_kesari = round(allBurstFreqsAllElectrodes_kesari{5}(allBurstLengthsAllElectrodes_kesari{5}<lengthLimit),1);
f_OMPGEAR_unique_kesari = unique(round(allBurstFreqsAllElectrodes_kesari{5}(allBurstLengthsAllElectrodes_kesari{5}<lengthLimit),1));
[medianLength_OMPGEAR_kesari,seLengthOMPGEAR_kesari] = getmlse(f_OMPGEAR_kesari,f_OMPGEAR_unique_kesari,x_OMPGEAR_kesari);
errorbar(f_OMPGEAR_unique_kesari,medianLength_OMPGEAR_kesari,seLengthOMPGEAR_kesari,'Color',color_OMPGEAR,'LineWidth',2);
hold off
set(get(ax7,'xlabel'),'string','Burst Centre Freq(Hz)','Fontsize',20);
set(get(ax7,'ylabel'),'string','Burst Duration(s)','Fontsize',20);
set(ax7,'Fontsize',20)

%%% plotting violin plot alpa %%

%ax4 = nexttile;
y_alpa = [medianBurstLength_alpa(1,:),medianBurstLength_alpa(2,:),medianBurstLength_alpa(3,:),medianBurstLength_alpa(4,:),medianBurstLength_alpa(5,:)];
g_alpa = [repmat(m1,1,size(medianBurstLength_alpa(1,:),2)) repmat(m2,1,size(medianBurstLength_alpa(2,:),2)) repmat(m3,1,size(medianBurstLength_alpa(3,:),2)) repmat(m4,1,size(medianBurstLength_alpa(4,:),2)) repmat(m5,1,size(medianBurstLength_alpa(5,:),2)) ];
x_alpa =  [ones(size(medianBurstLength_alpa(1,:))) 2*ones(size(medianBurstLength_alpa(2,:))) 3*ones(size(medianBurstLength_alpa(3,:))) 4*ones(size(medianBurstLength_alpa(4,:))) 5*ones(size(medianBurstLength_alpa(5,:)))];
g_alpa = cellstr(g_alpa);
pl_alpa = gramm('x',g_alpa','y',y_alpa','color',g_alpa');
pl_alpa.set_order_options('x', 0,'color',0);
pl_alpa.stat_violin('fill','transparent','dodge',0.75,'half',false);
pl_alpa.stat_boxplot('width',0.2, 'notch', true, 'dodge', 0.75)
pl_alpa.set_color_options('n_color',5,'map',col_map);
pl_alpa.set_text_options('base_size',20)
pl_alpa.axe_property('XTickLabelRotation',60)
%pl_alpa.set_names('x','Algorithms','y','Median Burst Length')
figure;f_alpa = pl_alpa.fig(pl_alpa.draw());


%%% plotting violin plot kesari %%
clear pl_alpa;
%ax8 = nexttile;
y_kesari = [medianBurstLength_kesari(1,:),medianBurstLength_kesari(2,:),medianBurstLength_kesari(3,:),medianBurstLength_kesari(4,:),medianBurstLength_kesari(5,:)];
g_kesari = [repmat(m1,1,size(medianBurstLength_kesari(1,:),2)) repmat(m2,1,size(medianBurstLength_kesari(2,:),2)) repmat(m3,1,size(medianBurstLength_kesari(3,:),2)) repmat(m4,1,size(medianBurstLength_kesari(4,:),2)) repmat(m5,1,size(medianBurstLength_kesari(5,:),2)) ];
x_kesari =  [ones(size(medianBurstLength_kesari(1,:))) 2*ones(size(medianBurstLength_kesari(2,:))) 3*ones(size(medianBurstLength_kesari(3,:))) 4*ones(size(medianBurstLength_kesari(4,:))) 5*ones(size(medianBurstLength_kesari(5,:)))];
g_kesari = cellstr(g_kesari);
pl_kesari = gramm('x',g_kesari','y',y_kesari','color',g_kesari');
pl_kesari.set_order_options('x', 0,'color',0);

pl_kesari.stat_violin('fill','transparent','dodge',0.75,'half',false);
pl_kesari.stat_boxplot('width',0.2, 'notch', true, 'dodge', 0.75)
pl_kesari.set_color_options('n_color',5,'map',col_map);
pl_kesari.set_text_options('base_size',20)
pl_kesari.axe_property('XTickLabelRotation',60)
pl_kesari.set_names('x','Algorithms','y','Median Burst Length')

figure;pl_kesari.draw();









