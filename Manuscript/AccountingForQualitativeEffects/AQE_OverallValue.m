function [betasNfix,tStats] = AQE_OverallValue()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function was used in the 'Accounting for Qualitative Effects in
% Previous Work' section. It computes the relationship between overall
% value and predicted RT whilst controlling for value difference using a
% GLM. It also varies the threshold parameter to test the hypothesis that
% the relationship between OV and RT is stronger when decisions are made less cautiously.
% (in MASC's case, higher threshold means less cautious.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Setting up figure and settings.
tiledlayout(1,3)
nSubj = 100;
nTrials = 200;
settings.m = 1;
settings.n = 5;
m = settings.m; n = settings.n;
settings.maxSteps = 100;


%% Simulate a dataset.
[dat,parameters] = MASC_Simulate(nSubj,nTrials,n,m);
f.dat = dat; f.parameters = parameters;


%% Range of threshold values to test.
threshToTest = linspace(.001,.2,3);


for iTest = 1:length(threshToTest)
    betasNfix = nan(nSubj,3);
    predicted_RT = nan(nTrials,nSubj);
    parameters = f.parameters;
    attValues = f.dat.attValues;
    w = f.parameters.w;
    parameters.thresh = zeros(nSubj,1) + threshToTest(iTest);
    %  parameters.threshInc =  zeros(nSubj,1) + .02;


    % Get overall value for all trials
    diffValue = nan(nTrials,nSubj);
    overallValue = nan(nTrials,nSubj);
    for s = 1:nSubj
        for t = 1:nTrials
            optValues = attValues(:,:,t,s)*w(:,s);
            diffValue(t,s) = max(optValues)-max(optValues(max(optValues)~=optValues)); %best vs. 2nd-best
            overallValue(t,s) = sum(optValues);
        end
    end


    %% Simulate Model
    [~, RT, ~] = switchModelCall(settings,parameters,attValues, nTrials, nSubj, 1);


    %% Get Betas
    for s = 1:nSubj
        zVD = (diffValue(:,s)-mean(diffValue(:,s)))./std(diffValue(:,s));
        zOV = (overallValue(:,s)-mean(overallValue(:,s)))./std(overallValue(:,s));
        betasNfix(s,:) = glmfit([zVD,zOV],RT(:,s));
        predicted_RT(:,s)= glmval(betasNfix(s,:)', [zVD, zOV], 'identity');
    end



    %% T-Test to check if betas are significantly different from zero and negative.
    [~,p,~,tStats] = ttest(betasNfix(:,3),zeros(nSubj,1),'tail','left');
    tStats.p = p;



    %% Plot one participant.
    nexttile(iTest)
    % Regression line when controlling for value difference.
    predicted_RT_valContr = glmval(betasNfix(end,:)', [repmat(mean(zVD),200,1), zOV], 'identity');

    ov = overallValue(:,end);
    plot(ov,predicted_RT(:,end)','.','Color',[.5 .5 .5],'HandleVisibility','off');
    hold on
    plot(ov,predicted_RT_valContr,'-','Color',[0 .75 0],'LineWidth',2);
    xlabel('Overall Value')
    ylabel('Predicted Number of Fixations');
    xlim([-6 6])
    [r,p] = corr(overallValue(:,end),predicted_RT(:,end));
    if p < .001
        legend(sprintf('r = %.3f, p < .001',r))
    else
        legend(sprintf('r = %.3f, p = %.3f',r,p));
    end
    title(sprintf('Threshold = %.2f',threshToTest(iTest)))
end
