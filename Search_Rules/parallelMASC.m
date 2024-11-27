function [choice,RT,allFix] = parallelMASC(settings,parameters,attValues,model,nTrials,nSubj)

%inFunction = runSingleMASC;


[choice, RT] = deal(nan(nTrials,nSubj));
allFix = nan(settings.maxSteps, nTrials,nSubj);
futures = parallel.FevalFuture.empty(nSubj, 0);

for s = 1:nSubj

    futures(s)  = parfeval(@runSingleMASC,3,settings,parameters,attValues(:,:,:,s),s,model,nTrials);
end



for iter = 1:nSubj
    [completedIdx, choiceS, RTS, allFixS] = fetchNext(futures);
    choice(:, completedIdx) = choiceS;
    RT(:, completedIdx) = RTS;
    allFix(:, :, completedIdx) = allFixS;

end




end