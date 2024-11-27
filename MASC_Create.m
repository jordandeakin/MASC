function [dat,parameters] = MASC_Create(model,nSubj,nTrials,n,m,parameters)
%MASC_Create(model,nSubj,nTrials,n,m,parameters)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MASC = Multi-Attribute Search and Choice Model %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Creates a set of options and parameters %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% Inputs:
%%% model -> search rule
%%% nSubj -> number of subjects (if dataset = 0; default = 64)
%%% nTrials -> number of trials (if dataset = 0; default = 128)
%%% n -> number of options (if dataset = 0; default = 2)
%%% m -> number of attributes (if dataset = 0; default = 3)
%%% parameters -> structure with fixed parameters; free will be added
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% sebastian.gluth@uni-hamburg.de %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%start with parameters (because you need w to define difficulty)
for s = 1:nSubj
    pre_w = betarnd(3/4,3/4,1,m); %note: in case of m = 3, about 7% of attributes will have a weight < 1%, which may be plausible
    parameters.w(:,s) = pre_w./sum(pre_w); %attribute weights
    parameters.sn(s,:) = repmat(0+rand*3,1,m); %sampling noise
    parameters.threshInc(s,1) = .04+randn*.015; %threshold increase
    parameters.threshDec(s,1) = .04+randn*.015; %threshold decrease (for absolute model)
    if ((model < 2) || (model > 3)) %no third parameter for RANDOM I and II models
        parameters.searchSense(s,1) = (rand*10); %search sensitivity
    end
end

%now, create choice sets
attValues = zeros(n,m,nTrials,nSubj); %attribute values
valDiff12 = zeros(n,nTrials,nSubj);
difficulty = zeros(nTrials,nSubj); %1 = easy, 2 = medium, 3 = hard
for s = 1:nSubj
    for t = 1:nTrials
        isNoBrainer = 1;
           attSuggest = randn(n,m)*sqrt(1/parameters.lambdaPrior);
         
        if m >1
        while (isNoBrainer==1)
            attSuggest = randn(n,m)*sqrt(1/parameters.lambdaPrior);
            isNoBrainer = any(sum(attSuggest==max(attSuggest),2)==m); %is there an option that is best w.r.t. all attributes?
        end
        end
        attValues(:,:,t,s) = attSuggest;
        optValues = attValues(:,:,t,s)*parameters.w(:,s);
        valDiff12(t,s) = [1,-1,zeros(1,n-2)]*sort(optValues,'descend'); %value difference (best vs. 2nd-best)
        % valDiff1M(t,s) = max(optValues(:,t,s))-(sum(optValues(:,t,s))-max(optValues(:,t,s)))./(n-1); %value difference (best vs. mean)
    end
    difficulty(:,s) = 1+sum(quantile(valDiff12(:,s),[1/3,2/3])>valDiff12(:,s),2);
end
dat.attValues = attValues;
dat.difficulty = difficulty;

end