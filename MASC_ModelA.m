function [choice,RT,allFix] = MASC_ModelA(settings,parameters,attValues,s)
%[choice,RT,allFix] = MASC_ModelA(settings,parameters,attValues,s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MASC = Multi-Attribute Search and Choice Model %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Absolute threshold version of MASC  %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% Inputs (defaults):
%%% settings -> number of options and attributes, etc.
%%% parameters -> attribute weights, SD, etc.
%%% attValues -> true attribute (mean) values
%%% s -> current subject ID (required for individual parameters)
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
% % sn = parameters.sn(s); %sampling noise (SD)
sn = parameters.sn(s,:); %sampling noise (SD)
sv = (sn.^2); %sampling variance (could be pre-allocated)
sp = 1./sv; %sampling precision (could be pre-allocated)
thresh1 = parameters.thresh1(s);
thresh2 = parameters.thresh(s); %threshold for choosing
threshDec = parameters.threshDec(s); %decrease in threshold per cycle
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
    transitionMatrix = MASC_SearchRule_myopicA(n,m,w,w2,sp,thresh1,thresh2,searchSense,attPrecision,attMean);
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
    %     overlap = normcdf(thresh1,optMean(iFix),sqrt(optVar(iFix)));
    overlap = min(normcdf(thresh1,optMean,sqrt(optVar)));
    decisionMade = thresh2>overlap;

    %part 4: update threshold and cycle counter
    thresh1 = thresh1 - threshDec;
    t = t + 1;
    allFix(t) = currentFix;
end

%store information
choice = find(min(normcdf(thresh1,optMean,sqrt(optVar)))==normcdf(thresh1,optMean,sqrt(optVar)),1);
if length(choice)>1
    choice = randi(length(choice),1);
end
RT = t;

end