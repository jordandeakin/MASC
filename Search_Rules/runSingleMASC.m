function [choice, RT, allFix] = runSingleMASC(settings,parameters,attValues,s,model,nTrials)


choice = nan(nTrials,1); RT = choice; allFix = nan(settings.maxSteps,nTrials);



                        for t = 1:nTrials
                            [choice(t),RT(t),allFix(:,t)] = MASC_FastModel(settings,parameters,attValues(:,:,t),s,model);
                        end
    end