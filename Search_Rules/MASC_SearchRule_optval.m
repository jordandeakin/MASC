function transitionMatrix = MASC_SearchRule_optval(m,w,searchSense,attMean)
%transitionMatrix = MASC_SearchRule_optval(m,w,searchSense,attMean)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MASC = Multi-Attribute Search and Choice Model %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Search Rule: option values %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% For definition of inputs, see other BISMAD scripts
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% sebastian.gluth@uni-hamburg.de %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%determine fixation probabilities based on attribute values
optValOAP = repmat(attMean*w,1,m); %option value projected on all OAPs
transitionMatrix = exp(searchSense*optValOAP)./sum(sum(exp(searchSense*optValOAP)));

end