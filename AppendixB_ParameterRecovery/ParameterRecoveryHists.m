function ParameterRecoveryHists()
figure

paramNames= {'Sampling Noise \sigma',	'Threshold Increase \Delta',	'Search Sensitivity \alpha'};
datasets = {'Phone','Hotel'};
colors = [0 0 1; 0 .75 0];
tiledlayout(2,3)
plotIdx = [1:3; 4:6];

nSimulations = 100; % We did 100 parameter recovery simulations. 
nParam = length(paramNames);


[r,p] = deal(zeros(nSimulations,length(nParam),2));


for iDat = 1:2
    
    % Load in the true parameters. 
    trueFile = load(sprintf('ws_ModelFit_1_%s_model1_gridNum31',datasets{iDat}));
    trueParams = trueFile.bestParams;

    % Get recovered parameters
    recovFile = load(sprintf('RecoveredParameters_%s',datasets{iDat}));
    recoveredParams = recovFile.recovered;
    

    for iRepeat = 1:nSimulations
        currentSim = squeeze(recoveredParams(iRepeat,:,:));
        for iParam = 1:nParam
   [r(iRepeat,iParam,iDat), p(iRepeat,iParam,iDat)] = corr(trueParams(:,iParam), currentSim(:,iParam));
        end
    end

    for iParam = 1:3
        nexttile(plotIdx(iDat,iParam))
        h = histogram(r(:,iParam,iDat),'FaceColor',colors(iDat,:),'NumBins',10);
   
        xline(mean(r(:,iParam,iDat),'all'),'LineWidth',3,'Color',colors(iDat,:)*.5,'HandleVisibility','off')
mu = round(mean(r(:,iParam,iDat),'all'),3);
title({paramNames{iParam}, strcat('r = ', num2str(mu))});

        % title({paramNames{iParam},strcat('$\overline{r}$ = ',mu)}, 'Interpreter', 'latex');
    end
    drawnow
end


for iParam = 1:6
    nexttile(iParam)
    xlabel('Correlation Coefficient');
    xlim([.4 1]);
   % legend(legText{iParam,:});
end
