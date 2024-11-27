function transitionMatrix = MASC_SearchRule_myopic(n,m,w,w2,sp,thresh,searchSense,attPrecision,attMean)
%transitionMatrix = MASC_SearchRule_myopic(n,m,w,w2,sp,thresh,searchSense,attPrecision,attMean)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MASC = Multi-Attribute Search and Choice Model %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Search Rule: standard search rule of MASC %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% For definition of inputs, see other MASC scripts
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% sebastian.gluth@uni-hamburg.de %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%updates of non-fixated elements (can be done outside of loop)
newPrecisions = attPrecision+repmat(sp,n,1);
optMeans = attMean*w;
optVarsOld = (1./attPrecision)*w2;

%loop over every OAP
myopicScore = zeros(n,m);

if m > 1
for i = 1:n
    not_i = i~=(1:n)';
    for j = 1:m
        not_j = j~=(1:m);
        optVarNew = (1./newPrecisions(i,j))*w2(j)+(1./attPrecision(i,not_j))*w2(not_j);
        optMeanThresh = optMeans(not_i)-norminv(thresh,0,sqrt(optVarNew+optVarsOld(not_i)));
        sampleThresh = (newPrecisions(i,j)./w(j)*(optMeanThresh-attMean(i,not_j)*w(not_j))-attPrecision(i,j).*attMean(i,j))./sp(j);
        attSD = sqrt(1./attPrecision(i,j));
        myopicScore(i,j) = normcdf(attMean(i,j),max(sampleThresh),attSD);
    end
end

else
for i = 1:n
    not_i = i~=(1:n)';
    for j = 1:m
        optVarNew = (1./newPrecisions(i,j))*w2(j);
        optMeanThresh = optMeans(not_i)-norminv(thresh,0,sqrt(optVarNew+optVarsOld(not_i)));
        sampleThresh = (newPrecisions(i,j)./w(j)*(optMeanThresh)-attPrecision(i,j).*attMean(i,j))./sp(j);
        attSD = sqrt(1./attPrecision(i,j));
        myopicScore(i,j) = normcdf(attMean(i,j),max(sampleThresh),attSD);
    end
end


end

myopicScore(myopicScore == 0) = realmin;
myopicScoreS = myopicScore./sum(sum(myopicScore));
transitionMatrix = exp(searchSense*myopicScoreS)./sum(sum(exp(searchSense*myopicScoreS)));

end