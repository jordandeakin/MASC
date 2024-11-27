function [dat,parameters] = MASC_Simulate(nSubj,nTrials,n,m)

addpath Search_Rules/ %this is where the search rule scripts are stored


%fixed parameters
parameters.lambdaPrior = 1; %prior precision (true precision of standardized attribute values is always 1)
parameters.thresh = zeros(nSubj,1)+0.01; %initial threshold (currently fixed, but could be free; therefore defined per subj)
parameters.thresh1 = zeros(1,nSubj); %for absolute threshold variant
parameters.thresh2 = zeros(nSubj,1)+0.01; %for absolute threshold variant

    [dat,parameters] = MASC_Create(1,nSubj,nTrials,n,m,parameters); %create new dataset
    
end