function transitionMatrix = MASC_SearchRule_rand1(n,m)
%transitionMatrix = MASC_SearchRule_rand1(n,m,currentFix)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MASC = Multi-Attribute Search and Choice Model %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Search Rule: random %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% Inputs (defaults):
%%% n -> number of options
%%% m -> number of attributes
%%% currentFix -> currently fixated OAP ("option-attribute pair")
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% sebastian.gluth@uni-hamburg.de %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%generate transition probabilities (all equal except currently fixated)
transitionMatrix = zeros(n,m)+1./(n*m);

end