function transitionMatrix = MASC_SearchRule_rand2(n,m,currentFix)
%transitionMatrix = MASC_SearchRule_rand2(n,m,currentFix)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MASC = Multi-Attribute Search and Choice Model %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Search Rule: random without re-fixations %%%%%%%%%
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
numOAP = n*m; %number of OAPs
if isnan(currentFix)
    transitionMatrix = zeros(n,m)+1./(numOAP);
else
    transitionMatrix = zeros(n,m)+1./(numOAP-1);
    transitionMatrix(currentFix) = 0;
end

end