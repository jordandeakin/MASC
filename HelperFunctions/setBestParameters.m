function f = setBestParameters(f,m)

if f.model < 10

[r,c] = size(f.bestParams);
if c == 3
f.parameters.sn = zeros(f.nSubj,m) + f.bestParams(:,1);
f.parameters.threshInc = f.bestParams(:,2);
f.parameters.searchSense = f.bestParams(:,3);
else
f.parameters.threshInc = f.bestParams(:,1);
f.parameters.searchSense = f.bestParams(:,2);

end


else

f.parameters.sn = zeros(f.nSubj,m) + f.bestParams(:,1);
f.parameters.threshDec = f.bestParams(:,2);
f.parameters.searchSense = f.bestParams(:,3);



end
end