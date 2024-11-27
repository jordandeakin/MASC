function transitionMatrix = MASC_SearchRule_myopic2(n,m,w,w2,sp,thresh,searchSense,attPrecision,attMean)
%transitionMatrix = MASC_SearchRule_myopic2(n,m,w,w2,sp,thresh,searchSense,attPrecision,attMean)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MASC = Multi-Attribute Search and Choice Model %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Search Rule: more complex myopic search %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% For definition of inputs, see other BISMAD scripts
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% sebastian.gluth@uni-hamburg.de %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%updates of non-fixated elements (can be done outside of loop)
% % newPrecisions = attPrecision+sp;
newPrecisions = attPrecision+repmat(sp,n,1);
optMeans = attMean*w;
optVarsOld = (1./attPrecision)*w2;

%loop over every OAP
myopicScore = zeros(n,m);
for i = 1:n
    not_i = i~=(1:n)';
    for j = 1:m
        not_j = j~=(1:m);
        optVarNew = (1./newPrecisions(i,j))*w2(j)+(1./attPrecision(i,not_j))*w2(not_j);
        optMeanThresh = [optMeans(not_i)-norminv(thresh,0,sqrt(optVarNew+optVarsOld(not_i))),... %value required to choose i
                         optMeans(not_i)+norminv(thresh,0,sqrt(optVarNew+optVarsOld(not_i)))];   %value required to choose other options
        sampleThresh = (newPrecisions(i,j)./w(j)*(optMeanThresh-attMean(i,not_j)*w(not_j))-attPrecision(i,j).*attMean(i,j))./sp(j);
        attSD = sqrt(1./attPrecision(i,j));
        myopicScore(i,j) = normcdf(attMean(i,j),max(sampleThresh(:,1)),attSD)+1-normcdf(attMean(i,j),max(sampleThresh(:,2)),attSD);
    end
end
myopicScore = myopicScore./sum(sum(myopicScore));
transitionMatrix = exp(searchSense*myopicScore)./sum(sum(exp(searchSense*myopicScore)));

end