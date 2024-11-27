function significantEffects = MASC_AnalyzeStudy(printFig,settings,parameters,dat,choice,RT,allFix,dataset)
%MASC_AnalyzeStudy(settings,parameters,attValues,choice,RT,allFix,dataset)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MASC = Multi-Attribute Search and Choice Model %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Analysis and illustration of model predictions %%%
%%% tailored to compare the model with study data %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% Inputs (defaults):
%%% printFig -> whether figures should be printed or not
%%% settings -> number of options and attributes, etc.
%%% parameters -> attribute weights, SD, etc.
%%% dat -> dataset information (attribute values, difficulty)
%%% choice -> chosen option
%%% RT -> number of cycles to reach choice ("Response Time")
%%% allFix -> fixations of option-attribute pairs (OAP)
%%% dataset -> HotelHotel (~= 2, default) or PhonePhone (= 2) dataset
%%%
%%% Outputs:
%%% significantEffects -> indicates which effects (H1-H7) are significant
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% sebastian.gluth@uni-hamburg.de %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

namesCondition = {'Phone','Hotel','Simulation'};
dataset(dataset==0) = 3;
dsName = namesCondition{dataset};

%settings
attValues = dat.attValues;
difficulty = dat.difficulty;
nSubj = size(attValues,4);
nTrials = size(attValues,3);
n = settings.n; %number of options
m = settings.m; %number of attributes
maxSteps = settings.maxSteps; %max number of steps

%parameters
w = parameters.w; %attribute weights


%%% --- %%% --- %%% --- %%% --- %%% --- %%% --- %%% --- %%% --- %%% ---
% 1 & 2) Consistency and RT for 3 difficulty levels
%%% --- %%% --- %%% --- %%% --- %%% --- %%% --- %%% --- %%% --- %%% ---

%define choice Consistency
consistentChoice = zeros(nTrials,nSubj);
for s = 1:nSubj
    optValues = zeros(n,nTrials); %option values
    for t = 1:nTrials
        optValues(:,t) = attValues(:,:,t,s)*w(:,s);
    end
    consistentChoice(:,s) = choice(:,s)==((1:n)*(repmat(max(optValues),n,1)==optValues))'; %whether the option with the highest value was chosen
end

%difficulty levels
CC_l = [(sum(consistentChoice.*(difficulty==1))./sum(difficulty==1))',...
        (sum(consistentChoice.*(difficulty==2))./sum(difficulty==2))',...
        (sum(consistentChoice.*(difficulty==3))./sum(difficulty==3))'];
CC_ci = norminv(.975)*std(CC_l)./sqrt(nSubj);
RT_l = [(sum(RT.*(difficulty==1))./sum(difficulty==1))',...
        (sum(RT.*(difficulty==2))./sum(difficulty==2))',...
        (sum(RT.*(difficulty==3))./sum(difficulty==3))'];
RT_ci = norminv(.975)*std(RT_l)./sqrt(nSubj);

figure;
subplot(1,2,1);hold on;title('MASC: Choice consistency');xlim([.5,3.5]);ylim([.3,1])
plot(mean(CC_l),'k.-','MarkerSize',20,'LineWidth',3)
plot(mean(CC_l)+CC_ci,'k-')
plot(mean(CC_l)-CC_ci,'k-')
plot([.6,3.4],[.5,.5],'k--')
set(gca,'XTick',1:3,'XTickLabel',{'Easy','Medium','Hard'})
xlabel('Difficulty');ylabel('Consistency')
subplot(1,2,2);hold on;title('MASC: Number of fixations');xlim([.5,3.5]);ylim([3,21])
plot(mean(RT_l),'b.-','MarkerSize',20,'LineWidth',3)
plot(mean(RT_l)+RT_ci,'b-')
plot(mean(RT_l)-RT_ci,'b-')
set(gca,'XTick',1:3,'XTickLabel',{'Easy','Medium','Hard'},'YTick',0:5:20)
xlabel('Difficulty');ylabel('RT in ms')
set(gcf,'PaperPosition',[1/4,1/4,24,12])
if printFig == 1
    print(gcf,'-dpng','-r300',[pwd,'/Figures/',dsName,'_Choice_and_RT'])
end

%stats
[~,p12] = ttest(CC_l(:,1),CC_l(:,2),'tail','right');
[~,p13] = ttest(CC_l(:,1),CC_l(:,3),'tail','right');
[~,p23] = ttest(CC_l(:,2),CC_l(:,3),'tail','right');
significantEffects.test1 = (p12<.05)&(p13<.05)&(p23<.05);
[~,p12] = ttest(RT_l(:,1),RT_l(:,2),'tail','left');
[~,p13] = ttest(RT_l(:,1),RT_l(:,3),'tail','left');
[~,p23] = ttest(RT_l(:,2),RT_l(:,3),'tail','left');
significantEffects.test2 = (p12<.05)&(p13<.05)&(p23<.05);


%%% --- %%% --- %%% --- %%% --- %%% --- %%% --- %%% --- %%% --- %%% ---
% 3) Association of attention and choice
%%% --- %%% --- %%% --- %%% --- %%% --- %%% --- %%% --- %%% --- %%% ---

matrixOAP = reshape(1:n*m,n,m); %matrix of option-attribute pairs
maxFix = sum(mean(mean(~isnan(allFix),2),3)>.2); %maximum fixation number to be shown (current rule: 20% of trials should include a fixation)
ncat = 8; %number of categories 
catPropFix = zeros(nSubj,ncat);
propFix1 = nan+zeros(nSubj,nTrials);
for s = 1:nSubj
    for t = 1:nTrials
        nfix = sum(~isnan(allFix(:,t,s))); %number of fixations in current trial
        if nfix > 1 %otherwise, the "NoLast" analysis can't be done
            nfix1 = sum(ismember(allFix(:,t,s),matrixOAP(1,:))); %number of fixations on first option
            propFix1(s,t) = nfix1/nfix;
        end
    end
    sortedProp = sortrows([propFix1(s,:)',(1:nTrials)'],1);
    for c = 1:ncat
        catTrials = (nTrials/ncat*(c-1)+1):(nTrials/ncat*c);
        catPropFix(s,c) = mean(choice(sortedProp(catTrials,2),s)==1);
    end
end
[~,~,ciPropFix] = ttest(catPropFix); %confidence intervals

figure;hold on;title('MASC: Gaze time and choice');xlim([.5,.5+ncat]);ylim([0,1])
plot(mean(catPropFix),'.-','LineWidth',3,'Color',[.2,.6,.6],'MarkerSize',20);
plot(ciPropFix(1,:),'-','Color',[.2,.6,.6]);
plot(ciPropFix(2,:),'-','Color',[.2,.6,.6]);
plot([.6,.4+ncat],[1/2,1/2],'k--')
xlabel('Proportion of gaze time on option 1');ylabel('P(choice option 1)')
set(gcf,'PaperPosition',[1/4,1/4,18,12])
if printFig == 1
    print(gcf,'-dpng','-r300',[pwd,'/Figures/',dsName,'_Gaze_Choice'])
end

%stats
betas = zeros(nSubj,3);
for s = 1:nSubj
    valDiff = zeros(nTrials,1); %control for value difference in this test
    for t = 1:nTrials
        optValues = attValues(:,:,t,s)*w(:,s);
        valDiff(t) = optValues(1)-mean(optValues(2:end));
    end
    warning('off')
    lastwarn('')
    betas(s,:) = glmfit([propFix1(s,:)',valDiff],choice(:,s)==1,'binomial');
    warningMessage = lastwarn;
    if ~isempty(warningMessage)
        betas(s,:) = nan;
    end
end
[~,p] = ttest(betas(:,2),zeros(nSubj,1),'tail','right');
significantEffects.test3 = p<.05;


%%% --- %%% --- %%% --- %%% --- %%% --- %%% --- %%% --- %%% --- %%% ---
% 4 & 5) Fixation development and value difference; ASE
%%% --- %%% --- %%% --- %%% --- %%% --- %%% --- %%% --- %%% --- %%% ---

betas = zeros(nSubj,4);
valFixCorr = zeros(nSubj,maxSteps)+nan;
for s = 1:nSubj
    glmFixNum = [];
    glmValDiff = [];
    glmFix1 = [];
    for t = 1:nTrials
        optValues = attValues(:,:,t,s)*w(:,s);
        valDiff = optValues(1)-mean(optValues(2:end));
        nfix = sum(~isnan(allFix(:,t,s)));
        glmFix1 = [glmFix1;ismember(allFix(1:nfix,t,s),matrixOAP(1,:))];
        glmValDiff = [glmValDiff;zeros(nfix,1)+valDiff];
        glmFixNum = [glmFixNum;(1:nfix)'];
    end
    %for figure
    maxFixGLM = sum(sum(repmat(1:max(glmFixNum),length(glmFixNum),1)==repmat(glmFixNum,1,max(glmFixNum)))>=(nTrials*.2));
    for f = 1:maxFixGLM
        valFixCorr(s,f) = corr(glmValDiff(f==glmFixNum),glmFix1(f==glmFixNum));
    end

    %stats
    glmFixNum = (glmFixNum-mean(glmFixNum))./std(glmFixNum);
    glmValDiff = (glmValDiff-mean(glmValDiff))./std(glmValDiff);
    glmInteract = glmFixNum.*glmValDiff;
    betas(s,:) = glmfit([glmFixNum,glmValDiff,glmInteract],glmFix1,'binomial');
end
vFC_ci = norminv(.975)*nanstd(valFixCorr)./sqrt(nSubj);

figure;hold on;title('Fixation and value difference');xlim([.5,.5+20]);ylim([-.5,.5])
plot(nanmean(valFixCorr),'-','LineWidth',3,'Color',[0,0,.5]);
plot(nanmean(valFixCorr)+vFC_ci,'-','Color',[0,0,.5])
plot(nanmean(valFixCorr)-vFC_ci,'-','Color',[0,0,.5])
plot([.6,.4+20],[0,0],'k--')
xlabel('Fixation number');ylabel('Correlation coefficient')
set(gcf,'PaperPosition',[1/4,1/4,18,12])
if printFig == 1
    print(gcf,'-dpng','-r300',[pwd,'/Figures/',dsName,'_Fixations_ValueDifference'])
end

%stats
[~,p] = ttest(betas(:,3:4),zeros(nSubj,2),'tail','right');
significantEffects.test4 = sum(p<.05)==length(p);

%ASE = Attraction Search Effect
firstAttValuePos = nan+zeros(nSubj,nTrials);
stayFix = nan+zeros(nSubj,nTrials);
for s = 1:nSubj
    for t = 1:nTrials
        if ~isnan(allFix(2,t,s)) %you need at least the first two fixations
            currentAttValues = attValues(:,:,t,s);
            firstAttValuePos(s,t) = currentAttValues(allFix(1,t,s))>0; %is first attribute value positive
            stayFix(s,t) = min(sum(allFix(1,t,s)==matrixOAP,2)==sum(allFix(2,t,s)==matrixOAP,2))==1; %does 2nd fixation remain on same option?
        end
    end
end

ASE = [sum((stayFix==1).*(firstAttValuePos==1),2)./sum((firstAttValuePos==1),2),...
       sum((stayFix==1).*(firstAttValuePos==0),2)./sum((firstAttValuePos==0),2),...
       sum((stayFix==0).*(firstAttValuePos==1),2)./sum((firstAttValuePos==1),2),...
       sum((stayFix==0).*(firstAttValuePos==0),2)./sum((firstAttValuePos==0),2)];

ASE_ci = norminv(.975)*std(ASE)./sqrt(nSubj);

figure;hold on;title('MASC: Attraction Search Effect');xlim([.75,4.25]);ylim([0,1])
p1 = plot(1:2,mean(ASE(:,1:2)),'.-','Color',[.2,.2,.8],'MarkerSize',20,'LineWidth',3);
p2 = plot(3:4,mean(ASE(:,3:4)),'.-','Color',[.2,.8,.2],'MarkerSize',20,'LineWidth',3);
plot(1:2,mean(ASE(:,1:2)+ASE_ci(1:2)),'-','Color',[.2,.2,.8])
plot(1:2,mean(ASE(:,1:2)-ASE_ci(1:2)),'-','Color',[.2,.2,.8])
plot(3:4,mean(ASE(:,3:4)+ASE_ci(3:4)),'-','Color',[.2,.8,.2])
plot(3:4,mean(ASE(:,3:4)-ASE_ci(3:4)),'-','Color',[.2,.8,.2])
set(gca,'XTick',1:4,'XTickLabel',{'Positive','Negative','Positive','Negative'})
xlabel('Sign of 1st attribute');ylabel('P(transition)')
legend([p1,p2],{'Stay','Switch'}); legend boxoff
set(gcf,'PaperPosition',[1/4,1/4,12,12])
if printFig == 1
    print(gcf,'-dpng','-r300',[pwd,'/Figures/',dsName,'_ASE'])
end

%stats
[~,p] = ttest(ASE(:,1),ASE(:,2),'tail','right');
significantEffects.test5 = p<.05;


%%% --- %%% --- %%% --- %%% --- %%% --- %%% --- %%% --- %%% --- %%% ---
% 6) Fixation development and attribute weights
%%% --- %%% --- %%% --- %%% --- %%% --- %%% --- %%% --- %%% --- %%% ---

fixAttWeight = zeros(m,maxSteps,nTrials,nSubj);
for s = 1:nSubj
    [~,attRank] = sort(w(:,s),'descend'); %the ith entry of attRank is the ith-important attribute
    for t = 1:nTrials        
        for j = 1:m %loop over attributes (from most to least important)
            fixAttWeight(j,:,t,s) = ismember(allFix(:,t,s),matrixOAP(:,attRank(j)))...
                ./~isnan(allFix(:,t,s));
        end
    end
end
fixAttWeight = reshape(nanmean(fixAttWeight,3),m,maxSteps,nSubj);

FixAttrMax_s = reshape(fixAttWeight(1,:,:),maxSteps,nSubj)';
FixAttrMid_s = reshape(fixAttWeight(2,:,:),maxSteps,nSubj)';
FixAttrMin_s = reshape(fixAttWeight(3,:,:),maxSteps,nSubj)';
FixAttrMax_ci = norminv(.975)*nanstd(FixAttrMax_s)./sqrt(nSubj);
FixAttrMid_ci = norminv(.975)*nanstd(FixAttrMid_s)./sqrt(nSubj);
FixAttrMin_ci = norminv(.975)*nanstd(FixAttrMin_s)./sqrt(nSubj);

figure;hold on;title('MASC: Fixation by attribute weights');xlim([.5,.5+20]);ylim([0,1])
p1 = plot(nanmean(FixAttrMax_s),'.-','Color',[0,.25,0],'LineWidth',3,'MarkerSize',20);
p2 = plot(nanmean(FixAttrMid_s),'.-','Color',[0,.75,0],'LineWidth',3,'MarkerSize',20);
p3 = plot(nanmean(FixAttrMin_s),'.-','Color',[.25,1,.25],'LineWidth',3,'MarkerSize',20);
plot(nanmean(FixAttrMax_s)+FixAttrMax_ci,'-','Color',[0,.25,0])
plot(nanmean(FixAttrMax_s)-FixAttrMax_ci,'.-','Color',[0,.25,0])
plot(nanmean(FixAttrMid_s)+FixAttrMid_ci,'-','Color',[0,.75,0])
plot(nanmean(FixAttrMid_s)-FixAttrMid_ci,'.-','Color',[0,.75,0])
plot(nanmean(FixAttrMin_s)+FixAttrMin_ci,'-','Color',[.25,1,.25])
plot(nanmean(FixAttrMin_s)-FixAttrMin_ci,'.-','Color',[.25,1,.25])
plot([.6,.4+20],[1/3,1/3],'k--')
xlabel('Fixation number');ylabel('P(fix)')
legend([p1,p2,p3],{'Highest weight','Medium weight','Lowest weight'},'Location','NorthEast');legend boxoff
set(gcf,'PaperPosition',[1/4,1/4,18,12])
if printFig == 1
    print(gcf,'-dpng','-r300',[pwd,'/Figures/',dsName,'_Fixations_Weights'])
end

%stats
betas = zeros(nSubj,2);
for s = 1:nSubj
    glmFixNum = nan+zeros(nTrials,maxFix);
    glmHighestFix = nan+zeros(nTrials,maxFix);
    for t = 1:nTrials
        [~,attRank] = sort(w(:,s),'descend');
        %nfix1 = sum(~isnan(allFix(:,t,s)))-1;
        nfix = sum(~isnan(allFix(:,t,s)));
        glmFixNum(t,1:nfix) = 1:nfix;
        glmHighestFix(t,1:nfix) = ismember(allFix(1:nfix,t,s),matrixOAP(:,attRank(1)));
    end
    zFixNum = (glmFixNum(glmFixNum>0)-mean(glmFixNum(glmFixNum>0)))./std(glmFixNum(glmFixNum>0));
    betas(s,:) = glmfit(zFixNum,glmHighestFix(glmFixNum>0),'binomial');
    betas(s,:) = glmfit(glmFixNum(glmFixNum>0)-1,glmHighestFix(glmFixNum>0),'binomial');
end
chanceLevelIntercept = log((1/m)/(1-(1/m))); %chance level of the logistic intercept for m attributes
[~,p1] = ttest(betas(:,1),zeros(nSubj,1)+chanceLevelIntercept,'tail','right');
[~,p2] = ttest(betas(:,2),zeros(nSubj,1),'tail','left');
significantEffects.test6 = (p1<.05)&(p2<.05);


%%% --- %%% --- %%%
% 7) Correlation of weight extremeness and Payne Index
%%% --- %%% --- %%%

transAlt = nan+zeros(size(allFix));
transAtt = nan+zeros(size(allFix));
for s = 1:nSubj
    for t = 1:nTrials
        if sum(~isnan(allFix(:,t,s)))>1 %you need to have at least 2 fixations
            for f = 1:sum(isnan(allFix(:,t,s)))-1
                transAlt(f,t,s) = max(sum((matrixOAP==allFix(f,t,s))+(matrixOAP==allFix(f+1,t,s)),2))==2;
                transAtt(f,t,s) = (max(sum((matrixOAP==allFix(f,t,s))+(matrixOAP==allFix(f+1,t,s))))==2)&(transAlt(f,t,s)~=1);
            end
        end
    end
end

aw = abs(w);
distWeight = (max(aw)-max((aw~=max(aw)).*aw))'; %difference of highest and second-highest weight per subject
PayneIndex = reshape((nansum(nansum(transAlt))-nansum(nansum(transAtt))),nSubj,1)./...
             reshape((nansum(nansum(transAlt))+nansum(nansum(transAtt))),nSubj,1); %standard Payne Index
b = glmfit(distWeight,PayneIndex);

figure;hold on;title('MASC: Weights and Payne Index');xlim([0,.7]);ylim([0,.7])
plot(distWeight,PayneIndex,'k.','MarkerSize',20)
plot(distWeight,distWeight*b(2)+b(1),'-','Color',[.2,.6,.6],'LineWidth',3)
xlabel('1st weight - 2nd weight');ylabel('Payne Index')
set(gcf,'PaperPosition',[1/4,1/4,12,12])
if printFig == 1
    print(gcf,'-dpng','-r300',[pwd,'/Figures/',dsName,'_PayneIndex_Weights'])
end

%stats
[~,p] = corr(distWeight(~isnan(PayneIndex)),PayneIndex(~isnan(PayneIndex)),'tail','left');
significantEffects.test7 = p<.05;

close all

end