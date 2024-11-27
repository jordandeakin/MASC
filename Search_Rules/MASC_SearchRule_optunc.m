function transitionMatrix = MASC_SearchRule_optunc(m,w2,searchSense,attPrecision)
%transitionMatrix = MASC_SearchRule_optunc(m,w2,searchSense,attPrecision)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MASC = Multi-Attribute Search and Choice Model %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Search Rule: option uncertainty %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% For definition of inputs, see other BISMAD scripts
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% sebastian.gluth@uni-hamburg.de %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%determine fixation probabilities based on attribute values
optVarOAP = repmat((1./attPrecision)*w2,1,m);
transitionMatrix = exp(searchSense*sqrt(optVarOAP))./sum(sum(exp(searchSense*sqrt(optVarOAP))));

end