function exampleSimulation()
close all
%% This function shows the example input expected by the function MASC_Model.
% parameters - a struct array containing the parameters for the model
% (including the attribute weights). 
% settings - a struct array containing the following fields (n = number of
% options, m = number of attrbutes, maxSteps = maximum iterations for the
% model).
% attValues = the attribute values for the presented options. 
% s = subject
% model = search rule to be used. See MASC_Model for different search
% rules. 

% MASC_Model(settings,parameters,attValues,s,model)

% settings struct
n = 2; % nOptions
m = 3; % nAttributes
settings.n = n;
settings.m = m;
settings.maxSteps = 100;

%fixed parameters
% Each entry corresponds to one subject.
model = 1; 
nSubj = 1;
nTrials = 100;
parameters.lambdaPrior = 1; %prior precision (true precision of standardized attribute values is always 1)
parameters.thresh = zeros(nSubj,1)+0.01; %initial threshold (currently fixed, but could be free; therefore defined per subj)
parameters.thresh1 = zeros(1,nSubj); %for absolute threshold variant
parameters.thresh2 = zeros(nSubj,1)+0.01; %for absolute threshold variant


% MASC_Create creates random data and parameters in the format expected
% by MASC. Adapt to your own data accordingly (below).
% This function adds the following fields
% w = attribute weights (here randomly sampled)
% sn = sampling noise parameter - each row is a subject, each column is the
% sampling noise for each attribute.
% threshInc = threshold increase parameter.
% thresDec = threshold decrease parameter used for absolute threshold
% variant.
% searchSense - Search Sensitivity Parameter. 
[dat,parameters] = MASC_Create(model,nSubj,nTrials,n,m,parameters); %create new dataset

% If you want to set different parameters
% In this case, we are simulating one participant so we will specify one
% value.
parameters.sn = ones(1,m); % Sampling noise = 1 for all attributes.
parameters.searchSense = 3; 
parameters.threshInc = 0.01;


%% Adding your own data. 
% Attvalues is a n x m x nTrials x nSubj matrix.
attValues = dat.attValues;

% The weights are a m x nSubj matrix. 
% This should be added to the parameters struct as parameters.w.
weights = parameters.w;

% The model outputs the choice, number of fixations and the fixation paths
% for each trial.
% The numbers in allFix (fixation paths) correspond to the numbers in
% OAPMatrix, where each row is an option and each column is an attribute. 
OAPmatrix = reshape(1:n*m,n,m); %"option-attribute pair"



for s = 1:nSubj
for t = 1:nTrials
[choice(t,s), nFix(t,s), allFix(:,t,s)] = MASC_Model(settings,parameters,attValues,s,model);

% Calculating the proportion of fixations to each attribute. 
% 1 3 5 = Option 1, 2 4 6 = Option 2. 
propFix(:,t,s) = histcounts(allFix(:,t,s),1:(n*m+1))./sum(histcounts(allFix(:,t,s),1:(n*m+1)));
end

%% Plotting distribution of number of fixations. 
figure(1)
nexttile(s)
for iChoice = 1:n
histogram(nFix(choice(:,s) == iChoice,s),'BinWidth',1);
hold on
end
title(sprintf('Subject %d',s))
sgtitle("Predicted Number of Fixations")
if s == 1
l = legend(strcat(repmat({'Choice = Option'}, 1, n),string(1:n)),'NumColumns',n);
l.Layout.Tile = 'south';
end


%% Plotting proportion of fixations to each OAP
mPropFix = mean(propFix(:,:,s),2);
figure(3)
nexttile(s)
for iOpt = 1:n
    propToPlot(iOpt,:) = mPropFix(OAPmatrix(iOpt,:))';
end
h = heatmap(propToPlot);
h.XDisplayLabels = strcat(repmat({'Att'}, 1, m),string(1:m));
h.YDisplayLabels = strcat(repmat({'Opt'}, 1, n),string(1:n));
h.ColorbarVisible = 'off';
title(sprintf('Subject %d',s))
end
sgtitle('Proportion of Fixations to Each OAP')






%% Plotting pChoice
figure(2)
for iChoice = 1:n
y(:,iChoice) = mean(choice == iChoice);
end
bar(1:nSubj, y,'stacked');
xlabel('Subject')
ylabel('p(Choice)')
legend(strcat(repmat({'Choice = Option'}, 1, n),string(1:n)))
sgtitle('Predicted Choice Probability')









