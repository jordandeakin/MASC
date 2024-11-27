function Extra_HicksLaw()
%%
% This function was used in the 'Accounting for Qualitative Effects in
% Previous Work' section. It calculates number of fixations (and
% approximated RT) across varying N. It plots this along with the function
% predicted by Hick's Law. 
%%

%% Colormap
clf
startColor = [194 218 184]/255;
HL = @(a,b,n) a + b * log2(n-1);


plotApprox = 0; % If plotApprox == 1, function will plot approximate RT by sampling from an exGauss


mu = 139.5530;
sigma = 63.1433;
tau = 91.0133;

% Simulate with 20 options and then add options...
nSubj = 100;
nTrials = 200;
maxOpt = 20;
[dat,parameters] = MASC_Simulate(nSubj,nTrials,maxOpt,1);
parameters.sn(parameters.sn < 0) = -parameters.sn(parameters.sn<0);
parameters.sn = zeros(nSubj,1) + 1;
parameters.threshInc =  zeros(nSubj,1) + .01; % Longer RTs.
parameters.searchSense = zeros(nSubj,1) + 5;


options = 2:1:maxOpt;
approxRT = deal(zeros(length(options),1));
settings.m = 1;
settings.maxSteps = 100;
[mRT, mApproxRT] = deal(nan(length(options),1));


for iOpt = 1:length(options)
    settings.n = options(iOpt);
    attValues = dat.attValues(1:options(iOpt),:,:,:);

    %% Simulate Model
        % This function switches the model call based on if you are using slurm or
% not (model is the same)
        [~, RT, ~] = switchModelCall(settings,parameters,attValues, nTrials, nSubj, 1);



    %% Sample from Ex-Gauss to get approximate RT.
    for iSubj = 1:nSubj
        for iTrial = 1:nTrials
            currRT = RT(iTrial,iSubj);
            gaussComp = normrnd(mu,sigma,currRT,1);
            expComp = exprnd(tau,currRT,1);
            approxRT(iTrial,iSubj) = sum(gaussComp + expComp);
        end
    end





    mApproxRT(iOpt) = mean(mean(approxRT));
    mRT(iOpt) = mean(mean(RT));


    if plotApprox == 0
        pd = polyfit(log2(options(1:iOpt)-1), mRT(1:iOpt),1);
        y = HL(pd(2),pd(1),1:maxOpt);
        plot(options(1:iOpt),mRT(1:iOpt),'o','MarkerFaceColor',[.5 .5 .5],'MarkerEdgeColor','none');hold on
        plot(1:maxOpt,y,'-','Color',startColor,'LineWidth',2);
        xlabel('Number of Options (N)');
        ylabel('Number of Fixations')
        hold off
        drawnow
    else
        pd = polyfit(log2(options(1:iOpt)-1), mApproxRT(1:iOpt),1);
        y = HL(pd(2),pd(1),1:maxOpt);
        plot(options(1:iOpt),mApproxRT(1:iOpt),'o','MarkerFaceColor',[.5 .5 .5],'MarkerEdgeColor','none');hold on
        plot(1:maxOpt,y,'-','Color',startColor,'LineWidth',2);
        xlabel('Number of Options (N)');
        ylabel('Reaction Time (ms)');
        hold off
        drawnow



    end
    legend('Number of Fixations','Hick''s Law: a + b * log_{2}(N-1)')
end
   