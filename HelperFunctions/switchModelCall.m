function [choice, RT, allFix] = switchModelCall(settings,parameters,attValues,nTrials,nSubj,model)

if ~isempty(getenv('SLURM_JOB_ID'))

    % parallelMASC function to be used with slurm.
    [choice,RT,allFix] = parallelMASC(settings,parameters,attValues,1,nTrials,nSubj);

else
     for s = 1:nSubj
        for t = 1:nTrials
            [choice(t,s),RT(t,s),allFix(:,t,s)] = MASC_Model(settings,parameters,attValues(:,:,t,s),s,model);
        end
    end

end
end