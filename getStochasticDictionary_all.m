% This program uses MP version 3.1 by Piotr Durka's group. This is the
% original stochastic dictionary code used in their 2001 paper.
% Modification: original code modified to include OMP and OMP-MAGE with
% option to select any of the three is algName 

function [gaborInfo,header] = getStochasticDictionary_all(data,timeVals,maxIteration,adaptiveDictionaryParam,dictionarySize,algName)
addpath ./Mage_try/
if ~exist('maxIteration','var');           maxIteration=50;             end
if ~exist('adaptiveDictionaryParam','var');adaptiveDictionaryParam=0.9; end
if ~exist('dictionarySize','var');      dictionarySize=[];              end
if ~exist('algName','var') ;            algName = "MP" ;                 end
Fs=round(1/(timeVals(2)-timeVals(1)));          % Sampling Frequency
numTrials = size(data,1);
sigLen = size(data,2);

%%%%%%%%%%%%%%%%%%%%%_________MP Parameters _________%%%%%%%%%%%%%%%%%%%
recAc=100; % Pecentage of Maximum reconstructed Energy (Value should be less than 100)

%%%%%%%%%%%%%%%%%%_______ Run MP on single trials ________%%%%%%%%%%%%%%
gaborInfo = zeros(numTrials,maxIteration,7);
header = zeros(numTrials,8);

for i=1:numTrials
    disp(['Trial ' num2str(i) ' of ' num2str(numTrials)]);
    if algName == "MP"
        % Save data as ASCII files
        sig=data(i,:)'; 
        save('sig.txt','sig','-ascii');
        clear sig;
        
        % Write Command File
        fp=fopen('commands.txt','w');
        if ~isempty(dictionarySize)
            fprintf(fp,['reinit -R ' num2str(dictionarySize) ' \n']);
        end
        fprintf(fp,['set -O ' num2str(sigLen) ' -M ' num2str(maxIteration) ' -E ' num2str(recAc) ' -F ' num2str(Fs) ' -D ' num2str(adaptiveDictionaryParam) ' \n']);
        fprintf(fp,'loadasc -O sig.txt\n');
        fprintf(fp,'mp\n'); % Run MP
        fprintf(fp,'save -S mpresults\n');
        fprintf(fp,'exit');        % Exit shell
        fclose(fp);
        
        % Run MP
        system(['./mp31' ' <commands.txt'],'-echo');
        
        % Read Data
        [gaborInfo(i,:,:), header(i,:)]=readbook('mpresults',0);
        
        % Delete tmp files
        %delete sig.txt;                         % Deletes the signal file
        %delete commands.txt;                    % Deletes the commamnds file created
        delete mpresults;                       % Deletes mpresults
    end
    if algName == "OMP"
        sig=data(i,:)'; 
        [a,h] = Omp_v1(sig,timeVals,i,dictionarySize,maxIteration);
        gaborInfo(i,:,:)=a;
        header(i,:) = h;
    end
    if algName == "OMP-MAGE"
        sig=data(i,:)'; 
        [a,h] = Omp_mage(sig,dictionarySize,Fs,maxIteration);
        gaborInfo(i,:,:)=a;
        header(i,:) = h;
    end
end
end