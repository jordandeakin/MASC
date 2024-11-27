function transitionMatrix = MASC_SearchRule_myopicA(n,m,w,w2,sp,thresh1,thresh2,searchSense,attPrecision,attMean)
%transitionMatrix = MASC_SearchRule_myopicA(n,m,w,w2,sp,thresh,searchSense,attPrecision,attMean)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MASC = Multi-Attribute Search and Choice Model %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Search Rule: myopic rule for absolute threshold %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% For definition of inputs, see other MASC scripts
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% sebastian.gluth@uni-hamburg.de %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%updates of non-fixated elements (can be done outside of loop)
newPrecisions = attPrecision+repmat(sp,n,1);

%loop over every OAP
myopicScore = zeros(n,m);
for i = 1:n
    for j = 1:m
        not_j = j~=(1:m);
        optVarNew = (1./newPrecisions(i,j))*w2(j)+(1./attPrecision(i,not_j))*w2(not_j);
        optMeanThresh = thresh1-norminv(thresh2,0,sqrt(optVarNew));
        sampleThresh = (newPrecisions(i,j)./w(j)*(optMeanThresh-attMean(i,not_j)*w(not_j))-attPrecision(i,j).*attMean(i,j))./sp(j);
        attSD = sqrt(1./attPrecision(i,j));
        myopicScore(i,j) = normcdf(attMean(i,j),sampleThresh,attSD);
    end
end
myopicScore = myopicScore./sum(sum(myopicScore));
transitionMatrix = exp(searchSense*myopicScore)./sum(sum(exp(searchSense*myopicScore)));

end