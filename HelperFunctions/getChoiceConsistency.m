function CC = getChoiceConsistency(studyData, choice)
[n,m,nTrials,nSub] = size(studyData.attValues);

difficulty = studyData.difficulty;
for s = 1:nSub
    for t = 1:nTrials
        optValues(:,t) = studyData.attValues(:,:,t,s)*studyData.attWeights(:,s);
    end
    consistent(:,s) = choice(:,s)==((1:n)*(repmat(max(optValues),n,1) == optValues))';
end
CC = mean(consistent);
% CC = [(sum(consistent.*(difficulty==1))./sum(difficulty==1))',...
%     (sum(consistent.*(difficulty==2))./sum(difficulty==2))',...
%     (sum(consistent.*(difficulty==3))./sum(difficulty==3))'];
% CC = mean(CC)*100;