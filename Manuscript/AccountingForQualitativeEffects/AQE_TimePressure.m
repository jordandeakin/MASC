function AQE_TimePressure()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function was used in the 'Accounting for Qualitative Effects in
% Previous Work' section.  It simulates how the following change with time
% pressure (manipulated through the threshold parameter:
%%% Payne Index
%%% Attribute Variance
%%% Option Variance
%%% Probability of fixating highest weight attribute
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Presettings for MASC %%
n = 2; % Number of Options
m = 3; % Number of Attributes
settings.m = m;
settings.n = n;
settings.maxSteps = 100; % Maximum timesteps for the model.
nSubj = 100; nTrials = 200; % 100 simulated participants, 200 trials each.
matrixOAP = reshape(1:n*m,n,m);
model = 1; % Which search rule to use. 1 = default myopic search rule. See MASC_Model for details of other search rules.

%% Load in existing dataset for reproducibility or create one if it doesn't exist.
if exist('AQE_TimePressure_Dat.mat','file') == 2
    f = load('AQE_TimePressure_Dat.mat');
    dat = f.dat; parameters = f.parameters;
else
    [dat,parameters] = MASC_Simulate(nSubj,nTrials,n,m);
    save('AQE_TimePressure_Dat','dat','parameters');
end

%% Get the attribute values and weights.
attValues = dat.attValues;
w = dat.attWeights;


%% Number of threshold parameters and levels of alpha (search sensitivity) to test.
nThresh = 10;
nSearch = 5;
threshToTest = linspace(.001,.2,nThresh);
senseToTest = linspace(0,3,nSearch);

%% Preallocation.
propFixAttr = nan(m,nTrials,nSubj);
propFixOpt = nan(n,nTrials,nSubj);
[pMostImp, pLeastImp, tempAttr, tempOpt] = deal(nan(nTrials,nSubj));
[pFixMax, pFixMin, attrVar, optVar, PI, CC, RTs, RR] = deal(nan(nSubj,nThresh));
tiledlayout(2,4)


for iSearch = 1:length(senseToTest)

    % Set search sensitivity.
    parameters.searchSense = zeros(nSubj,1) + senseToTest(iSearch);

    % Set threshold parameters.
    for iThresh = 1:length(threshToTest)
        parameters.thresh = zeros(nSubj,1) + threshToTest(iThresh);

        % This function adapts the function call for easier parallel use
        % with slurm.
        [choice, RT, allFix] = switchModelCall(settings,parameters,attValues, nTrials, nSubj, model);

        % For each subject, find the least and most important attribute
        % based on the weights.
        for iSub = 1:nSubj
            [~,mostImportant] = max(w(:,iSub));
            [~, leastImportant] = min(w(:,iSub));

            for iTrial = 1:nTrials
                nFix = sum(~isnan(allFix(:,iTrial,iSub)));

                % Calculate the proportion of fixations on each attribute.
                for iAttr = 1:m
                    propFixAttr(iAttr,iTrial,iSub) =  sum(ismember(allFix(:,iTrial,iSub),matrixOAP(:,iAttr))) / nFix;
                end

                % And on each option.
                for iOpt = 1:n
                    propFixOpt(iOpt,iTrial,iSub) = sum(ismember(allFix(:,iTrial,iSub),matrixOAP(iOpt,:))) /nFix;
                end

                % Calculate proportion of fixations on most vs least
                % important attrbute.
                pMostImp(iTrial,iSub) = propFixAttr(mostImportant,iTrial,iSub);
                pLeastImp(iTrial,iSub) = propFixAttr(leastImportant,iTrial,iSub);
            end

            % Find the variance in the proportion of fixations on each
            % option/attrbute.
            tempAttr(:,iSub) = var(propFixAttr(:,:,iSub),1);
            tempOpt(:,iSub) = var(propFixOpt(:,:,iSub),1);
        end


        % Averaging
        pFixMax(:,iThresh) = mean(pMostImp);
        pFixMin(:,iThresh) = mean(pLeastImp);


        attrVar(:,iThresh) = mean(tempAttr);
        optVar(:,iThresh) = mean(tempOpt);


        % Calculate the Payne Index
        [~, PayneIndex] = H7(allFix,nSubj,nTrials,matrixOAP,w);
        PI(:,iThresh) = PayneIndex;
        CC(:,iThresh) = getChoiceConsistency(dat,choice);
        RTs(:,iThresh) = mean(RT);
        RR(:,iThresh) = CC(:,iThresh)./RTs(:,iThresh);

    end







    %% Payne Index
    nexttile(1)
    plot(threshToTest, mean(PI));
    hold on


    %% Attribute Variance
    nexttile(2)
    plot(threshToTest,mean(attrVar));
    hold on


    %% Options Variance
    nexttile(3)
    plot(threshToTest,mean(optVar));
    hold on

    %% Time Spent on Most Important Attribute
    nexttile(4)
    plot(threshToTest, mean(pFixMax));
    hold on

    %% Choice Consistency
    nexttile(5)
    plot(threshToTest,mean(CC))
    hold on

    %% Number of Fixations
    nexttile(6)
    plot(threshToTest,mean(RTs))
    hold on

    %% Reward Rate
    nexttile(7)
    plot(threshToTest,mean(RR))
    hold on

    drawnow
end

savefig(gcf,'WeightsAttention')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define Color Map
startColor = [194 218 184]/255;
endColor = [1 50 32]/255;
colors = [linspace(startColor(1), endColor(1), nSearch)', ...
    linspace(startColor(2), endColor(2), nSearch)', ...
    linspace(startColor(3), endColor(3), nSearch)'];

labels = {'Payne Index','Attribute Variance','Option Variance','p(Fix = Most Important)','Choice Consistency','Number of Fixations','Reward Rate'};
letters = {'A','B','C','D','E','F','G'};


for iPlot = 1:length(labels)
    nexttile(iPlot)
    title(letters{iPlot});
    ylabel(labels{iPlot});
    xlabel('Threshold \theta');
    xlim([-.05 .25]);
    set(gca,'TitleHorizontalAlignment','left')

    h = flip(findobj(gca,'Type','Line'));
    for k = 1:length(h)
        set(h(k), 'Marker','o', 'MarkerFaceColor', colors(k,:), 'MarkerEdgeColor','none','LineWidth',2,'Color',colors(k,:));
    end
end

savefig(gcf,'WeightsAttention')


l = legend(string(senseToTest));
l.Title.String = "Search Sensitivity ";
l.Layout.Tile = 'south';
l.NumColumns = 5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [distWeight, PayneIndex] = H7(allFix,nSubj,nTrials,matrixOAP,w)
transAlt = nan+zeros(size(allFix));
transAtt = nan+zeros(size(allFix));
for ss = 1:nSubj
    for t = 1:nTrials
        if sum(~isnan(allFix(:,t,ss)))>1 %you need to have at least 2 fixations
            for ff = 1:sum(isnan(allFix(:,t,ss)))-1
                transAlt(ff,t,ss) = max(sum((matrixOAP==allFix(ff,t,ss))+(matrixOAP==allFix(ff+1,t,ss)),2))==2;
                transAtt(ff,t,ss) = (max(sum((matrixOAP==allFix(ff,t,ss))+(matrixOAP==allFix(ff+1,t,ss))))==2)&(transAlt(ff,t,ss)~=1);
            end
        end
    end
end

aw = abs(w);
distWeight = (max(aw)-max((aw~=max(aw)).*aw))'; %difference of highest and second-highest weight per subject
PayneIndex = reshape((nansum(nansum(transAlt))-nansum(nansum(transAtt))),nSubj,1)./...
    reshape((nansum(nansum(transAlt))+nansum(nansum(transAtt))),nSubj,1); %standard Payne Index
%b = gl