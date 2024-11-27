function MASC_Initialize(dataset,printFig,model,gridNum,pRecov,nSubj,nTrials,n,m)
%MASC_Initialize(dataset,printFig,model,gridNum,pRecov,nSubj,nTrials,n,m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MASC = Multi-Attribute Search and Choice Model %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initializes simulation of MASC (either with %%%%%%
%%% real options and parameters or generates them) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% Inputs:
%%% dataset -> Phone (= 1, default), Hotel (= 2), none = 0 (generate data)
%%% printFig -> whether figures should be printed or not (default = 0 = no)
%%% model -> search rule that should be used; the search rules are...
%%%     ... 1 = myopic: the "core" model of the project (default) 
%%%     ... 2 = rand1: random search with direct re-fixations
%%%     ... 3 = rand2: random search without direct re-fixations
%%%     ... 4 = attval: accumulated attribute values determine search
%%%     ... 5 = optval: accumulated option values determine search
%%%     ... 6 = attunc: uncertainty about attribute values determine search
%%%     ... 7 = optunc: uncertainty about option values determine search
%%%     ... 8 = attweight: weight of attributes determine search
%%%     ... 9 = myopic2: the "core" rule + consider choosing other options
%%%     ... 10 = myopicA: the "core" rule but absolute threshold
%%% gridNum -> number grid search steps when model was fit (default = 21)
%%% pRecov -> whether simulations are for parameter recovery (default = 0)
%%% nSubj -> number of subjects (if dataset = 0; default = 64)
%%% nTrials -> number of trials (if dataset = 0; default = 120)
%%% n -> number of options (if dataset = 0; default = 2)
%%% m -> number of attributes (if dataset = 0; default = 3)
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% sebastian.gluth@uni-hamburg.de %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath Search_Rules/ %this is where the search rule scripts are stored

%assign default values if not specified otherwise
if exist('dataset','var') == 0
    dataset = 1;
end
if exist('printFig','var') == 0
    printFig = 0;
end
if exist('model','var') == 0
    model = 1;
end
if exist('gridNum','var') == 0
    gridNum = 21;
end
if exist('pRecov','var') == 0
    pRecov = 0;
end
if dataset == 0
    if exist('nSubj','var') == 0
        nSubj = 64;
    end
    if exist('nTrials','var') == 0
        nTrials = 120;  %NOTE: must be divisible by 3
    end
    if exist('n','var') == 0
        n = 2;
    end
    if exist('m','var') == 0
        m = 3;
    end
end

%%% --- %%% --- %%%
% Load or create dataset and parameters
%%% --- %%% --- %%%

%load/create dataset
if dataset > 0
    namesCondition = {'Phone','Hotel'};
    preDat = load(['studyData_',namesCondition{dataset},'.mat']);
    dat = preDat.studyData;
    
    %settings
    attValues = dat.attValues;
    nSubj = size(attValues,4);
    nTrials = size(attValues,3);
    n = size(attValues,1); %number of options
    m = size(attValues,2); %number of attributes
end

%fixed parameters
parameters.lambdaPrior = 1; %prior precision (true precision of standardized attribute values is always 1)
parameters.thresh = zeros(nSubj,1)+0.01; %initial threshold (currently fixed, but could be free; therefore defined per subj)
parameters.thresh1 = zeros(1,nSubj); %for absolute threshold variant
parameters.thresh2 = zeros(nSubj,1)+0.01; %for absolute threshold variant

%load/create free parameters
if dataset > 0
    parameters.w = dat.attWeights; %attribute weights (taken from rating tasks, so actually fixed)
    modelFit = load(['./ManuscriptWorkspaces/ws_ModelFit_',namesCondition{dataset},'_model',int2str(model),'_gridNum',int2str(gridNum)]);
    if pRecov == 1 %at a bit noise for parameter recovery (only for illustration purposes)
        pRecovNoise = 0.1; %amount of noise in relation to group SD of parameter estimates
        jitterRecov = randn(size(modelFit.bestParams)).*repmat(std(modelFit.bestParams)*pRecovNoise,nSubj,1);
        modelFit.bestParams = modelFit.bestParams+jitterRecov;
    end
    parameters.sn = zeros(nSubj,m)+modelFit.bestParams(:,1); %sampling noise
    parameters.threshInc = modelFit.bestParams(:,2); %threshold increase
    parameters.threshDec = modelFit.bestParams(:,2); %threshold decrease (for absolute model)
    if ((model < 2) || (model > 3)) %no third parameter for RANDOM I and II models
        parameters.searchSense = modelFit.bestParams(:,3); %search sensitivity
    end
else %if dataset does not exist, create it here together with free parameters
    [dat,parameters] = MASC_Create(model,nSubj,nTrials,n,m,parameters); %create new dataset
    attValues = dat.attValues;
end

%%% --- %%% --- %%%
% Simulate the model
%%% --- %%% --- %%%

%prepare loop
settings = struct('n',n,'m',m,'maxSteps',100); %max number of steps (currently: fixations)
choice = zeros(nTrials,nSubj);
RT = zeros(nTrials,nSubj);
allFix = zeros(settings.maxSteps,nTrials,nSubj);

%loop over subjects and trials
if model < 10 %run different model script for relative vs. absolute threshold
    parfor s = 1:nSubj %for s = 1:nSubj
        for t = 1:nTrials
            [choice(t,s),RT(t,s),allFix(:,t,s)] = MASC_Model(settings,parameters,attValues(:,:,t,s),s,model);
        end
    end
else %absoulte threshold models
    parfor s = 1:nSubj % for s = 1:nSubj
        for t = 1:nTrials
            [choice(t,s),RT(t,s),allFix(:,t,s)] = MASC_ModelA(settings,parameters,attValues(:,:,t,s),s);
        end
    end
end
        
%%% --- %%% --- %%%
% Analyze the model
%%% --- %%% --- %%%
if pRecov == 1 %if dataset is for parameter recovery, save summary statistics here
    studyData.attWeights = dat.attWeights;
    studyData.attValues = dat.attValues;
    studyData.sumStats.choice = choice;
    studyData.sumStats.RT = RT;
    propFix = nan+zeros(n*m,nTrials,nSubj);
    for s = 1:nSubj
        for t = 1:nTrials
            propFix(:,t,s) = histcounts(allFix(:,t,s),1:(n*m+1))./sum(histcounts(allFix(:,t,s),1:(n*m+1)));
        end
    end
    studyData.sumStats.propFix = propFix;
    save([pwd,'/studyData_',namesCondition{dataset},'_pRecov.mat'],'studyData'); %not sure why but pwd required to see saved files
else
    % MASC_Analyze(settings,parameters,attValues,choice,RT,allFix,dataset); %for general analyses of any type of dataset
    significantEffects = MASC_AnalyzeStudy(printFig,settings,parameters,dat,choice,RT,allFix,dataset); %for analyses specific to n = 2, m = 3
end

end