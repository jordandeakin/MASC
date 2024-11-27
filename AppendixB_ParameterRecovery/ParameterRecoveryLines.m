function ParameterRecoveryLines()


paramNames= {'Sampling Noise \sigma',	'Threshold Increase \Delta',	'Search Sensitivity \alpha'};
datasets = {'Phone','Hotel'};
colors = [0 0 1; 0 .75 0];
tiledlayout(2,3)
plotIdx = [1:3; 4:6];

nSimulations = 100; % We did 100 parameter recovery simulations.
nParam = length(paramNames);



for iDat = 1:2

    % Load in the true parameters.
    trueFile = load(sprintf('ws_ModelFit_1_%s_model1_gridNum31',datasets{iDat}));
    trueParams = trueFile.bestParams;

    % Get recovered parameters
    recovFile = load(sprintf('RecoveredParameters_%s',datasets{iDat}));
    recoveredParams = recovFile.recovered;


    % Just for display purposes - adding jitter so that the points don't
    % overlap
    pRecovNoise = 0.1; %amount of noise in relation to group SD of parameter estimates
    jitterRecov = randn(size(trueFile.bestParams)).*repmat(std(trueFile.bestParams)*pRecovNoise,63,1);
    trueJit = trueFile.bestParams+jitterRecov;
    trueJit(trueJit<0) = trueFile.bestParams(trueJit<0);




    for iParam = 1:nParam
        x = trueFile.bestParams(:,iParam);
        xx = linspace(0,10,1000);
        nexttile(plotIdx(iDat, iParam))

        for i = 1:nSimulations
            y = recoveredParams(i,:,iParam);
            pd = polyfit(x,y,1);
            plotY(i,:) = polyval(pd, xx); % area plot 
        end



        mi = min(plotY, [], 1); mi = mi';
        ma = max(plotY, [], 1); ma = ma';

        % Plot the area and average of recovered parameters
        fill([sort(xx) fliplr(sort(xx))], [sort(mi); flipud(sort(ma))], [0.5 0.5 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        hold on
        plot(trueJit(:,iParam), mean(recoveredParams(:,:,iParam)),'.','Color',colors(iDat,:));
        drawnow
    end
end

% Limits for plotting
lims = [0 3.5; 0 .07; 0 10];
for i = 1:3
    nexttile(i)
    xlim(lims(i,:));
    ylim(xlim);
    xlabel('True Parameter Value')
    ylabel('Recovered Parameter Value');
    title(paramNames{i})

    nexttile(i+3)
    xlim(lims(i,:));
    ylim(xlim)
    xlabel('True Parameter Value');
    ylabel('Recovered Parameter Value');
    title(paramNames{i})
end




