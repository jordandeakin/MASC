function t = MASC_RewardRate_Grid(gridNum)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MASC = Multi-Attribute Search and Choice Model %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Uses parameter grid to assess relationship between parameters and
%%% number of fixations, choice consistency and reward rate.  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Fixed Parameters
nTrials = 200;
model = 1;
nSubj = gridNum^3;
n = 2;
m = 3;
parameters.lambdaPrior = 1;
settings.n = n; settings.m = m; settings.maxSteps = 100;
paramNames = {'Sampling Noise (\sigma)','Threshold Change (\Delta)','Search Sensitivity (\alpha)'};
fieldNames = {'sn','threshInc','searchSense'};


% Preallocation
[choice, RT] = deal(zeros(nTrials,nSubj));
allFix = zeros(settings.maxSteps,nTrials,nSubj);
consistency = nan(nTrials, nSubj);
propFix = nan+zeros(n*m,nTrials,nSubj);


% Use CreateRR to get some parameter values, weights and attribute values.
[dat,parameters] = MASC_Create(model,nSubj,nTrials,n,m,parameters);
parameters.thresh = zeros(nSubj,1)+0.01; %initial threshold (currently fixed, but could be free; therefore defined per subj)
parameters.thresh1 = zeros(1,nSubj); %for absolute threshold variant
parameters.thresh2 = zeros(nSubj,1)+0.01; %for absolute threshold variant
dat.attWeights = parameters.w;
attValues = dat.attValues;


% Creating equidistant parameter values.
sn = linspace(.001,3.001,gridNum);
threshInc = linspace(.01,.09,gridNum); % Small values of threshInc lead to strange results so starting at .01.
searchSense = linspace(0,10,gridNum);

% Creating 3d grid of parameter values.
[snGrid, thGrid, ssGrid] = ndgrid(sn, threshInc, searchSense);
samples = [snGrid(:), thGrid(:), ssGrid(:)];
% Add some noise...
samples = samples + 0.001 * randn(size(samples));
% Negative parameter values have caused problems... flip any that have
% become negative due to the noise added.
samples(samples<0) = -samples(samples<0);

% Assign new sampled parameters.
parameters.sn = repmat(samples(:,1),1,3);
parameters.threshInc = samples(:,2);
parameters.searchSense = samples(:,3);


% Simulate the model with the given attribute values/parameters
parfor s = 1:nSubj
    for t = 1:nTrials
        [choice(t,s),RT(t,s),allFix(:,t,s)] = MASC_Model(settings,parameters,attValues(:,:,t,s),s,1);
    end
end

dat.sumStats.choice = choice;
dat.sumStats.RT = RT;


% Calculate fixation proportion and consistency
for s = 1:nSubj
    optValues = zeros(n,nTrials); %option values
    for t = 1:nTrials
        propFix(:,t,s) = histcounts(allFix(:,t,s),1:(n*m+1))./sum(histcounts(allFix(:,t,s),1:(n*m+1)));
        optValues(:,t) = dat.attValues(:,:,t,s)*dat.attWeights(:,s);
    end
    consistency(:,s) = dat.sumStats.choice(:,s)==((1:n)*(repmat(max(optValues),n,1)==optValues))';
end

dat.sumStats.propFix = propFix;


RT = mean(dat.sumStats.RT,'omitnan')';
CC = mean(consistency)';
RR = CC./RT; %reward rate

% Remove outlier participants
incl_s = true(nSubj,1);
CC = CC(incl_s);
RT = RT(incl_s);
RR = RR(incl_s);


%% GLM Fit
predictors = [parameters.sn(incl_s,1) parameters.threshInc(incl_s) parameters.searchSense(incl_s)];
[bCC,~,stats_CC] = glmfit(predictors,CC);
[bRT,~,stats_RT] = glmfit(predictors,RT);
[bRR,~,stats_RR] = glmfit(predictors,RR);

res.bCC.pred = bCC;
res.bCC.stats = stats_CC;

res.bRT.pred = bRT;
res.bRT.stats = stats_RT;

res.bRR.pred = bRR;
res.bRR.stats = stats_RR;


for iParam = 1:3

    currPred = predictors(:,iParam);
    minVal = min(currPred); maxVal = max(currPred);
    xPadding = (maxVal - minVal) * .1;

    % Correlations
    [r,p] = corr(currPred,CC);
    res.(fieldNames{iParam}).CC = [r, p];
    [r,p] = corr(currPred,RT);
    res.(fieldNames{iParam}).RT = [r, p];
    [r,p] = corr(currPred,RR);
    res.(fieldNames{iParam}).RR = [r, p];


    %draw subplots for figure (which has been initiated outside the loop)
    subplot(3,3,iParam*3-2);hold on;
    if ismember(iParam*3-2,1:3)
        title('Choice Consistency');
    end
    xlabel(paramNames{iParam});ylabel('Consistency');%set(gca,'YTick',.5:.1:1)
    plot(currPred,CC,'.','MarkerSize',3,'Color',[.5 .5 .5],'LineWidth',3);
    hold on
    plot([minVal maxVal],bCC(1)+bCC(iParam+1)*[minVal, maxVal],...
        '-','LineWidth',2,'Color','g')
    xlim([minVal-xPadding maxVal+xPadding])
    ylim([.5 1])
    yticks([.5 .75 1]);


    subplot(3,3,iParam*3-2+1);hold on;
    if ismember(iParam*3-2+1,1:3)
        title('Number of Fixations');
    end
    xlabel(paramNames{iParam});ylabel('Number of Fixations');%set(gca,'YTick',2:4:20)
    plot(currPred,RT,'.','MarkerSize',3,'Color',[.5 .5 .5],'LineWidth',3);
    hold on
    plot([minVal maxVal],bRT(1)+bRT(iParam+1)*[minVal, maxVal],...
        '-','LineWidth',2,'Color','g')
    xlim([minVal-xPadding maxVal+xPadding])
    ylim([0 40])
    yticks([0 20 40]);


    subplot(3,3,iParam*3-2+2);hold on;
    if ismember(iParam*3-2+2,1:3)
        title('Reward Rate');
    end
    xlabel(paramNames{iParam});ylabel('Reward Rate');%set(gca,'YTick',.02:.04:.2)
    plot(currPred,RR,'.','MarkerSize',3,'Color',[.5 .5 .5],'LineWidth',3);
    hold on
    plot([minVal maxVal],bRR(1)+bRR(iParam+1)*[minVal, maxVal],...
        '-','LineWidth',2,'Color','g')
    xlim([minVal-xPadding maxVal+xPadding])
    ylim([0 .4]);
    yticks([0 .2 .4])

end




set(findall(groot,'Type','Axes'),'FontName','Calibri','FontSize',14)

t = table();
t.names = ['Intercept',fieldNames]';
t.betaCC = stats_CC.beta;

t.df = repmat(stats_CC.dfe,4,1);
t.tCC = stats_CC.t;
t.pCC = stats_CC.p;


t.betaRT = stats_RT.beta;
t.tRT = stats_RT.t;
t.pRT = stats_RT.p;

t.betaRR = stats_RR.beta;
t.tRR = stats_RR.t;
t.pRR = stats_RR.p;



end

