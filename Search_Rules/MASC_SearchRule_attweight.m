function transitionMatrix = MASC_SearchRule_attweight(n,w,searchSense)
%transitionMatrix = MASC_SearchRule_attweight(n,w)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MASC = Multi-Attribute Search and Choice Model %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Search Rule: attribute weights %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% For definition of inputs, see other BISMAD scripts
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% sebastian.gluth@uni-hamburg.de %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%determine fixation probabilities based on attribute values
wOAP = repmat(w',n,1);
transitionMatrix = exp(searchSense*wOAP)./sum(sum(exp(searchSense*wOAP)));

end