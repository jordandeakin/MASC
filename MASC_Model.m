function [choice,RT,allFix] = MASC_Model(settings,parameters,attValues,s,model)
%[choice,RT,allFix] = MASC_Model(settings,parameters,attValues,s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MASC = Multi-Attribute Search and Choice Model %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Basic MASC model function %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Inputs (defaults):
%%% settings -> number of options and attributes, etc.
%%% parameters -> attribute weights, SD, etc.
%%% attValues -> true attribute (mean) values
%%% s -> current subject ID (required for individual parameters)
%%% model -> search-rule model that should be used
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% sebastian.gluth@uni-hamburg.de %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% --- %%% --- %%%
% Specifications
%%% --- %%% --- %%%

%settings
n = settings.n; %number of options
m = settings.m; %number of attributes
maxSteps = settings.maxSteps; %max number of steps (currently: fixations)

%parameters
w = parameters.w(:,s); %attribute weights
w2 = w.^2; %squared attribute weights (could be pre-allocated)
lambdaPrior = parameters.lambdaPrior; %prior precision
sn = parameters.sn(s,:); %sampling noise (SD)
sv = (sn.^2); %sampling variance (could be pre-allocated)
sp = 1./sv; %sampling precision (could be pre-allocated)
thresh = parameters.thresh(s); %threshold for choosing
threshInc = parameters.threshInc(s); %increase in threshold per cycle
searchSense = parameters.searchSense(s); %sensitivity of search

%prior belief distributions
attPrecision = zeros(n,m)+lambdaPrior; %precision of attribute BD
attMean = zeros(n,m); %mean of attribute BD

%%% --- %%% --- %%%
% Simulate decision
%%% --- %%% --- %%%

t = 0; %cycle counter
decisionMade = 0;
OAPmatrix = reshape(1:n*m,n,m); %"option-attribute pair"
allFix = nan+zeros(maxSteps,1);
currentFix = nan; %currently fixated OAP ("option-attribute pair")
while (decisionMade == 0) && (t < maxSteps)

    %part 1: generate fixation
    if model == 1
        [transitionMatrix] = MASC_SearchRule_myopic(n,m,w,w2,sp,thresh,searchSense,attPrecision,attMean);
    elseif model == 2
        transitionMatrix = MASC_SearchRule_rand1(n,m);
    elseif model == 3
        transitionMatrix = MASC_SearchRule_rand2(n,m,currentFix);
    elseif model == 4
        transitionMatrix = MASC_SearchRule_attval(searchSense,attMean);
    elseif model == 5
        transitionMatrix = MASC_SearchRule_optval(m,w,searchSense,attMean);
    elseif model == 6
        transitionMatrix = MASC_SearchRule_attunc(searchSense,attPrecision);
    elseif model == 7
        transitionMatrix = MASC_SearchRule_optunc(m,w2,searchSense,attPrecision);
    elseif model == 8
        transitionMatrix = MASC_SearchRule_attweight(n,w,searchSense);
    elseif model == 9
        transitionMatrix = MASC_SearchRule_myopic2(n,m,w,w2,sp,thresh,searchSense,attPrecision,attMean);
    end
    currentFix = sum(rand>cumsum(transitionMatrix(1:n*m)))+1;
    [~,jFix] = find(OAPmatrix==currentFix);

    %part 2: draw a sample and update belief distributions
    currentSample = attValues(currentFix)+randn*sn(jFix);
    newPrecision = attPrecision(currentFix)+sp(jFix); %necessary to use both the old and new precision in the mean-update (next line)
    attMean(currentFix) = (currentSample*sp(jFix) + attMean(currentFix)*attPrecision(currentFix))./newPrecision;
    attPrecision(currentFix) = newPrecision;
    optMean = attMean*w; %(note: could be faster to only update the currently fixated option, esp. in the case of many options)
    optVar = (1./attPrecision)*w2; %(note: could be faster to only update the currently fixated option, esp. in the case of many options)

    %part 3: check whether threshold has been crossed
    currentBest = optMean==max(optMean); %identify currently best option
    if sum(currentBest)==1 %only makes sense to check threshold if there is one option better than all the others
        overlaps = normcdf(0,optMean(currentBest)-optMean(~currentBest),sqrt(optVar(currentBest)+optVar(~currentBest)));
        decisionMade = sum(thresh>overlaps)==(n-1);
    end

    %part 4: update threshold and cycle counter
    thresh = thresh + threshInc;
    t = t + 1;

    allFix(t) = currentFix;
end

%store information
choice = find(currentBest,1);
RT = t;

end