function transitionMatrix = MASC_SearchRule_attunc(searchSense,attPrecision)
%transitionMatrix = MASC_SearchRule_attunc(searchSense,attPrecision)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MASC = Multi-Attribute Search and Choice Model %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Search Rule: attribute uncertainty %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% For definition of inputs, see other BISMAD scripts
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% sebastian.gluth@uni-hamburg.de %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%determine fixation probabilities based on attribute values
attSD = 1./sqrt(attPrecision);
transitionMatrix = exp(searchSense*attSD)./sum(sum(exp(searchSense*attSD)));

end