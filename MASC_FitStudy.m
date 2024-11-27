function MASC_FitStudy(dataset,model,gridNum,pRecov,ID)
%MASC_FitStudy(dataset,model,gridNum,pRecov)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MASC = Multi-Attribute Search and Choice Model %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Fit MASC to actual study data via grid search %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% Inputs:
%%% dataset -> Phone (= 1, default), Hotel (= 2)
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
%%% gridNum -> number of steps of grid search per parameter (default = 21)
%%% pRecov -> use simulations for parameter recovery (default = 0)
%%% ID -> ID to add to the results file. Useful if want to run the same fit
%%% multiple times.
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% sebastian.gluth@uni-hamburg.de %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath Search_Rules/ %this is where the search rule scripts are stored

%assign default values if not specified otherwise
if exist('dataset','var') == 0
    dataset = 1;
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

%%% --- %%% --- %%%
% Specifications
%%% --- %%% --- %%%

%load dataset
namesCondition = {'Phone','Hotel'};
if pRecov == 1
    preDat = load(['studyData_',namesCondition{dataset},'_pRecov.mat']);
else
    preDat = load(['studyData_',namesCondition{dataset},'_wFix.mat']);
end
dat = preDat.studyData;

%settings
attValues = dat.attValues;
nSubj = size(attValues,4);
nTrials = size(attValues,3);
n = size(attValues,1); %number of options
m = size(attValues,2); %number of attributes
settings = struct('n',n,'m',m,'maxSteps',100); %max number of steps (currently: fixations)

%fixed parameters
parameters.lambdaPrior = 1; %prior precision (true precision of standardized attribute values is always 1)
parameters.w = dat.attWeights; %attribute weights (taken from rating tasks)
parameters.thresh = zeros(nSubj,1)+0.01; %initial threshold (currently fixed, but could be free; therefore defined per subj)
parameters.thresh1 = zeros(1,nSubj); %for absolute threshold variant
parameters.thresh2 = zeros(nSubj,1)+0.01; %for absolute threshold variant

%%% --- %%% --- %%%
% Grid search loop
%%% --- %%% --- %%%

sumStats = dat.sumStats; %summary statistics of dataset

%parameter ranges
paramRange.sn = linspace(.001,3.001,gridNum);
paramRange.threshInc = linspace(0,.09,gridNum); %relative decision thresholds increase
paramRange.threshDec = linspace(0,.27,gridNum); %absolute decision thresholds decrease
if (model == 2) || (model == 3) %basic models assume random attention allocation, so no free parameter
    paramRange.searchSense = nan+zeros(1,gridNum);
else
    paramRange.searchSense = linspace(0,10,gridNum);
end

%pre-allocation and running the loop
fitChoice = nan+zeros(gridNum,gridNum,gridNum,nSubj);
fitFixNum = nan+zeros(gridNum,gridNum,gridNum,nSubj);
fitFixProp = nan+zeros(gridNum,gridNum,gridNum,nSubj);


for p1 = 1:gridNum
    parameters.sn = zeros(nSubj,settings.m)+paramRange.sn(p1); %sampling noise
    for p2 = 1:gridNum
        parameters.threshInc = zeros(nSubj,1)+paramRange.threshInc(p2); %relative decision thresholds increase
        parameters.threshDec = zeros(nSubj,1)+paramRange.threshDec(p2); %absoulte decision thresholds decrease
        for p3 = 1:gridNum
            parameters.searchSense = zeros(nSubj,1)+paramRange.searchSense(p3);

            %simulate the model
            choice = zeros(nTrials,nSubj);
            RT = zeros(nTrials,nSubj);
            allFix = zeros(settings.maxSteps,nTrials,nSubj);

            %loop over subjects and trials
            if model < 10 %relative threshold models
                if ((model == 2) || (model == 3)) && (p3 > 1) %don't simulate random search models more than once
                    choice(:) = nan;
                    RT(:) = nan;
                    allFix(:) = nan;
                else

                    [choice, RT, allFix] = switchModelCall(settings,parameters,attValues,nTrials,nSubj,model);


                end
            else %absoulte threshold models
                for s = 1:nSubj % for s = 1:nSubj
                    for t = 1:nTrials
                        [choice(t,s),RT(t,s),allFix(:,t,s)] = MASC_ModelA(settings,parameters,attValues(:,:,t,s),s);
                    end
                end
            end

            %assess model fit via summary statistics
            for s = 1:nSubj
                fitChoice(p1,p2,p3,s) = mean(choice(:,s)==sumStats.choice(:,s));
                fitFixNum(p1,p2,p3,s) = 1-nanmean(abs(RT(:,s)-sumStats.RT(:,s))./max(sumStats.RT(:,s)));
                propFix = nan+zeros(n*m,nTrials);
                for t = 1:nTrials
                    propFix(:,t) = histcounts(allFix(:,t,s),1:(n*m+1))./sum(histcounts(allFix(:,t,s),1:(n*m+1)));
                    %  fitFix(p1,p2,p3,s) = dl_distance_mex(sumStats.allFix(:,t,s), allFix(:,t,s));
                end
                fitFixProp(p1,p2,p3,s) = 1-nanmean(nanmean(abs(propFix-sumStats.propFix(:,:,s))));

            end
        end
    end
    save([pwd,'/Workspaces/ws_ModelFit_currentRun']) %saves current status within outer loop (for breakdowns)
    disp(['Model fitting: done with outer loop ',int2str(p1)])
end

%%% --- %%% --- %%%
% Determine best fit
%%% --- %%% --- %%%

%find overall best fit
fitAll = (fitChoice+fitFixNum+fitFixProp)./3;
% fitAll = log(max(eps,fitChoice))+log(max(eps,fitFixNum))+log(max(eps,fitFixProp));
fitAllAvg = nanmean(fitAll,4);
[f1,f2,f3] = ind2sub(size(fitAllAvg),find(fitAllAvg==max(fitAllAvg(1:numel(fitAllAvg)))));
bestParamsAvg = [paramRange.sn(f1),paramRange.threshInc(f2),paramRange.searchSense(f3)];

%find best fit per participant
bestParams = zeros(nSubj,length(bestParamsAvg));
for s = 1:nSubj
    fitS = fitAll(:,:,:,s);
    [f1,f2,f3] = ind2sub(size(fitS),find(fitS==max(fitS(1:numel(fitS)))));
    if model < 10
        bestParams(s,:) = [paramRange.sn(f1),paramRange.threshInc(f2),paramRange.searchSense(f3)];
    else
        bestParams(s,:) = [paramRange.sn(f1),paramRange.threshDec(f2),paramRange.searchSense(f3)];
    end
end

%%% --- %%% --- %%%
% Save results (not sure why but pwd required to see saved files)
%%% --- %%% --- %%%
if pRecov == 1
    save([pwd,'/ParameterRecovery/ws_ModelFit_pRecov_',int2str(ID),'_',namesCondition{dataset},'_model',int2str(model),'_gridNum',int2str(gridNum)])
else
    save([pwd,'/Workspaces/ws_ModelFit_',int2str(ID),'_',namesCondition{dataset},'_model',int2str(model),'_gridNum',int2str(gridNum)])
end


end