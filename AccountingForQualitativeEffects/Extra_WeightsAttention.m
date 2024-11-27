function Extra_WeightsAttention()
%%
% This function was used in the 'Accounting for Qualitative Effects in
% Previous Work' section. It computes the relationship between attribute
% weights and attention given varying weight differences. It also plots the
% development of fixations over time when weight difference = 0.5 - as
% (roughly) reported in March & Gluth, 2024 & X. Yang et al., 2024).


% Range of noise levels to test.
noiseToTest = 0.5:.5:3;
weights = 0:.05:1;


% Setting up figure
startColor = [194 218 184]/255;
endColor = [1 50 32]/255;
colors = [linspace(startColor(1), endColor(1), length(noiseToTest))', ...
    linspace(startColor(2), endColor(2), length(noiseToTest))', ...
    linspace(startColor(3), endColor(3), length(noiseToTest))'];
tiledlayout(1,2)
nexttile(2)
plot([-1 1],[-1 1],'-','Color',[.7 .7 .7],'HandleVisibility','off');
hold on

% Load saved attribute weights & set fixed parameters
f = load('Extra_WeightsAttention_Dat');
[n,m,nTrials,nSubj] = size(f.attValues);
settings.m = m; settings.n = m; settings.maxSteps = 100;
attValues = f.attValues;
matrixOAP = reshape(1:n*m,n,m);
parameters.lambdaPrior = 1; %prior precision (true precision of standardized attribute values is always 1)
parameters.thresh = zeros(nSubj,1)+0.01; %initial threshold (currently fixed, but could be free; therefore defined per subj)
parameters.searchSense = zeros(nSubj,1) + 10; % High search sensitivity.




for iT = 1:length(noiseToTest)


    parameters.threshInc = zeros(nSubj,1) + .05;
    parameters.sn = zeros(nSubj,2)+noiseToTest(iT);
    [attentionDiff, weightDiff] = deal(nan(length(weights),1));

    for iW = 1:length(weights)

        %% Range of different weight differences.
        parameters.w(1,:) = zeros(nSubj,1) + weights(iW);
        parameters.w(2,:) =(1-parameters.w(1,:));
        weightDiff(iW) = parameters.w(1,1) - parameters.w(2,1);
        pAtt = nan(nTrials, n*m, nSubj);


        %% Simulate Model

        % This function switches the model call based on if you are using slurm or
        % not (model is the same)
        [~, RT, allFix] = switchModelCall(settings,parameters,attValues, nTrials, nSubj, 1);

        %% Proportion of fixations spent on each attribute.
        for s = 1:nSubj
            for t = 1:nTrials
                currTrial = allFix(:,t,s);
                currTrial = currTrial(~isnan(currTrial));
                pAtt(t,:,s) = histcounts(currTrial,1:(n*m+1))./sum(histcounts(currTrial,1:(n*m+1)));

            end
        end

        pAttM = squeeze(mean(pAtt, 1))';

        pAtt1 = pAttM(:,1) + pAttM(:,2);
        pAtt2 = pAttM(:,3) + pAttM(:,4);


        attentionDiff(iW) = mean(pAtt1 - pAtt2);


        %% If weight difference == 0.5 (as in March & Gluth, 2024 & X.Yang et al., 2024, then plot the development of fixations over time.
        if weightDiff(iW) == .5
            nexttile(1)
            maxSteps = 100;
            fixAttWeight = zeros(1,maxSteps,nTrials,nSubj);
            for s = 1:nSubj
                for t = 1:nTrials
                    for j = 1:m %loop over attributes
                        fixAttWeight(j,:,t,s) = ismember(allFix(:,t,s),matrixOAP(:,j))...
                            ./~isnan(allFix(:,t,s));
                    end
                end
            end
            fixAttWeight = reshape(nanmean(fixAttWeight,3),m,maxSteps,nSubj);

            % When weight difference == 0.5, the first attribute is the
            % most important.
            FixAttr1 = reshape(fixAttWeight(1,:,:),maxSteps,nSubj)';

            plot(nanmean(FixAttr1(:,1:max(RT(:)))),'-o','Color',colors(iT,:),'MarkerFaceColor',colors(iT,:),'MarkerEdgeColor','none')
            hold on
            xlabel('Fixation Number')%,'FontWeight','bold');
            ylabel('p(Fix)_{Most Important}')%,'FontWeight','bold')
        end

        drawnow
    end


    %% Plot the function relating weight difference to attention difference.
    nexttile(2)
    plot(weightDiff,attentionDiff,'-o','Color',colors(iT,:),'MarkerFaceColor',colors(iT,:),'MarkerEdgeColor','none');
    hold on
    ylim([-1 1])
    xlabel('Weight_{Att_{1}} - Weight_{Att_{2}}')%,'FontWeight','bold')
    ylabel('p(Fix)_{Att_{1}} - p(Fix)_{Att_{2}}')%,'FontWeight','bold')
    l = legend(string(noiseToTest),'NumColumns',6);
    l.Layout.Tile = 'south';
    l.Title.String = 'Sampling Noise';
end
