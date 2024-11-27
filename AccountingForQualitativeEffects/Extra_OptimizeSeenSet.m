function Extra_OptimizeSeenSet()
%%
% This function was used in the 'Accounting for Qualitative Effects in
% Previous Work' section. It computes the probability of having seen all
% options and the probability of choosing the best option (best being
% defined as sum of weighted attribute values) as a function of increasing
% N. MASC predicts that the relationship between set size and the
% probability of choosing the best option is dependent on noise, therefore
% three different noise levels are simulated.
%%

%% Define the noise levels to be tested.
noise = [.01 1 3];
maxOpt = 20;
options = 2:2:maxOpt;

%% Setting up figure (colors, tiles)
clf
nn = length(noise);
startColor = [194 218 184]/255;
endColor = [1 50 32]/255;
colors = [linspace(startColor(1), endColor(1), nn)', ...
    linspace(startColor(2), endColor(2), nn)', ...
    linspace(startColor(3), endColor(3), nn)'];
tiledlayout(1,2);


% Simulates a set of parameters and attribute values/weights.
% We first simulate 20 options and then choose N. 
nTrials = 200;
nSubj = 100;
m = 1; % Single Attribute
[dat,parameters] = MASC_Simulate(nSubj,nTrials,maxOpt,m);
settings.m = m;
settings.maxSteps = 100;


% For each noise level
for iNoise = 1:length(noise)
    parameters.sn = zeros(nSubj,m) + noise(iNoise);


   
    [nSeen, choseBest] = deal(nan(nSubj,length(options)));
    

    % In each iteration of the loop we append more options...

    for iOpt = 1:length(options)

        isSeen = nan(options(iOpt),nTrials);
        settings.n = options(iOpt);
        n = settings.n;
        attValues = dat.attValues(1:options(iOpt),:,:,:);

        matrixOAP = reshape(1:n*m,n,m); %matrix of option-attribute pairs





% This function switches the model call based on if you are using slurm or
% not (model is the same)
        [choice, ~, allFix] = switchModelCall(settings,parameters,attValues, nTrials, nSubj,1);


        for iSubj = 1:nSubj
            bestChosen = nan(nTrials,1);
            % Which Items were seen?
            for iSeen = 1:options(iOpt)
                fixs = allFix(:,:,iSubj);
                isSeen(iSeen,:) = sum(ismember(fixs,matrixOAP(iSeen,:))) > 0;
            end


            % Was the best item chosen out of those that were seen?
            for t = 1:nTrials
                optValues = attValues(1:options(iOpt),:,t,iSubj)*parameters.w(:,iSubj);
                bestSeen = find(optValues == max(optValues(logical(isSeen(:,t)))));
                bestChosen(t) = choice(t,iSubj)==bestSeen;
            end



            nSeen(iSubj,iOpt) = mean(mean(isSeen));
            choseBest(iSubj,iOpt) = mean(bestChosen);
        end

    end
    nexttile(1)
    plot(options(1:iOpt),mean(nSeen),'o-','Color',startColor,'MarkerFaceColor',colors(iNoise,:),'MarkerEdgeColor',colors(iNoise,:))
    hold on
    % xlabel('Number of Options (N)')
    ylabel('Proportion of Items Seen')
    ylim([0 1])
    xlim([0 22])
    xlabel('Number of Options (N)');

    nexttile(2)
    plot(options(1:iOpt),mean(choseBest),'o-','Color',startColor,'MarkerFaceColor',colors(iNoise,:),'MarkerEdgeColor',colors(iNoise,:));
    hold on
    %  xlabel('Number of Options (N)')
    ylabel('p(Choice = Best Seen Option)');
    ylim([0 1])
    xlim([0 22]);
    xlabel('Number of Options (N)');
    drawnow

end

t = findall(groot,'Type','TiledLayout');
%t.XLabel.String = 'Number of Options (N)';
t.XLabel.FontWeight = 'bold';


l = legend(string(noise));
l.Layout.Tile = 'south';
l.NumColumns = length(noise);
l.Title.String = 'Sampling Noise \sigma';