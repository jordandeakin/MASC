function res = runHypothesisTests_Data()
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preamble...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load in icons for phone and hotel.
% Imshow does some strange things with the icons so need to do some special
% stuff.
% phoneIcon = imread('PhoneIcon.png'); hotelIcon = imread('HotelIcon.png');
% phoneIcon = phoneIcon(:,:,3); maxRGB = max(phoneIcon(:));
% phoneIcon(phoneIcon == 0) = 255;
% phoneIcon = phoneIcon - maxRGB;
% phoneIcon(phoneIcon == 255-maxRGB) = 255;
% hotelIcon(hotelIcon == 0) = 255;
% % Position of Phone and Hotel Icons.
% phonePosition = [0.03125,0.799069767441861,0.048958333333333,0.125581395348837];
% hotelPosition = phonePosition;
% hotelPosition(2) = 0.3246511627907;

datasetName = {'Phone','Hotel'};
gridNum = 31;

make_lighter = @(color,factor) (color+([255 255 255]-color)*factor)/255;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop Over Datasets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iDat = 1:2
    % Loading in the model fit..
   % f = load(sprintf('./CurrentPaper/PaperFits/ws_ModelFit_1_%s_model1_gridNum%d.mat',datasetName{iDat},gridNum));
    f = load(sprintf('./Workspaces/ws_ModelFit_1_%s_model1_gridNum31.mat',datasetName{iDat}));
    d = load(sprintf('studyData_%s_wFix.mat',datasetName{iDat}));

    % Task Parameters
    nTrials = f.nTrials;
    [nSubj,~] = size(f.bestParams);
    n = f.n; mm = f.m;
    matrixOAP = reshape(1:n*mm,n,mm); %matrix of option-attribute pairs
    attValues = d.studyData.attValues;
    w = d.studyData.attWeights;
    difficulty = d.studyData.difficulty;
    dataChoice = d.studyData.sumStats.choice;
    dataRT = d.studyData.sumStats.RT;
    dataFix = d.studyData.sumStats.allFix;
    f = setBestParameters(f,mm);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% H1 & H2: Choice Consistency & Number of Fixations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [CC_l, CC_ci, RT_l, RT_ci] = H1(attValues,w,nTrials,nSubj,difficulty,dataChoice,dataRT,n);

    % Data
    % Choice Consistency
    figure(1);
    subplot(2,4,(iDat*4)-3)
    plot(CC_l','Color',make_lighter([0 0 0],.8))
    hold on
    errorbar(1:3,mean(CC_l),CC_ci,'ko-','LineWidth',2,'CapSize',0,'MarkerFaceColor','k','MarkerEdgeColor','none')
    yline(1/n,'--')
    xlim([0 4]); ylim([.2 1]);
    xticks(1:3);xticklabels({'Easy','Medium','Hard'}); xlabel('Difficulty','FontWeight','bold');
    title('Data: Choice Consistency'); ylabel('Consistency')

    % Number of Fixations
    subplot(2,4,(iDat*4)-2)
    plot(RT_l','Color',make_lighter([0 0 255],.8));
    hold on
    errorbar(1:3,mean(RT_l),RT_ci,'b','LineWidth',2,'CapSize',0,'MarkerFaceColor','b','MarkerEdgeColor','none')
    xlim([0 4]);ylim([0 20])
    xticks(1:3);xticklabels({'Easy','Medium','Hard'}); xlabel('Difficulty','FontWeight','bold');
    title('Data: Number of Fixations'); ylabel('Number of Fixations')

    [res.(datasetName{iDat}).data.H1, res.(datasetName{iDat}).data.H2] = computeH1_H2(CC_l,RT_l);

    % Model
    % Choice Consistency
    [choice,RT, allFix] = simulateMASC(f);
    [CC_l, CC_ci, RT_l, RT_ci] = H1(attValues,w,nTrials,nSubj,difficulty,choice,RT,n);
    subplot(2,4,(iDat*4)-1)
    plot(1:3,mean(CC_l),'ko-','LineWidth',2,'MarkerFaceColor','k','MarkerEdgeColor','none')
    hold on
    plot(1:3,mean(CC_l)+CC_ci,'-k'); plot(1:3,mean(CC_l)-CC_ci,'-k');
    xlim([0 4]); ylim([.2 1]);
    xticks(1:3);xticklabels({'Easy','Medium','Hard'}); xlabel('Difficulty','FontWeight','bold');
    yline(1/n,'--')
    title('MASC: Choice Consistency'); ylabel('Consistency')

    % Number of Fixations
    subplot(2,4,(iDat*4))
    plot(1:3,mean(RT_l),'o-b','LineWidth',2,'MarkerFaceColor','b','MarkerEdgeColor','none')
    hold on
    plot(1:3,mean(RT_l)+RT_ci,'-b'); plot(1:3,mean(RT_l)-RT_ci,'-b');
    xlim([0 4]);ylim([0 20])
    xticks(1:3);xticklabels({'Easy','Medium','Hard'}); xlabel('Difficulty','FontWeight','bold');
    title('MASC: Number of Fixations'); ylabel('Number of Fixations')
    [res.(datasetName{iDat}).model.H1, res.(datasetName{iDat}).model.H2] = computeH1_H2(CC_l,RT_l);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% H3: Dwell Time & Choice.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Data
    [catPropFix, ciPropFix,nCat,propFix1] = H3(dataFix, dataChoice, matrixOAP, nSubj, nTrials);
    figure(2)
    % Data
    subplot(2,2,(iDat*2)-1)
    plot(catPropFix','Color',make_lighter([0 128 128],.6));
    hold on
    errorbar(1:nCat,mean(catPropFix),ciPropFix,'o-','Color',[0 128 128]/255,'LineWidth',2,'MarkerFaceColor',[0 128 128]/255,'CapSize',0)
    xlim([0 9]); ylim([0 1])
    title('Data: Gaze Time & Choice'); xlabel('Binned Number of Fixations on Option 1'); ylabel('p(Choice = Option 1)');
    res.(datasetName{iDat}).data.H3 = computeH3(nSubj,nTrials, attValues, w, dataChoice,propFix1);


    % Model
    [catPropFix, ciPropFix,nCat, propFix1] = H3(allFix, choice, matrixOAP, nSubj, nTrials);
    subplot(2,2,iDat*2)
    plot(mean(catPropFix),'o-','Color',[0 128 128]/255,'LineWidth',2,'MarkerFaceColor',[0 128 128]/255)
    hold on
    plot(mean(catPropFix)+ciPropFix(1,:),'-','Color',[0 128 128]/255)
    plot(mean(catPropFix)-ciPropFix(1,:),'-','Color',[0 128 128]/255)
    xlim([0 9]); ylim([0 1]);
    yline(1/n,'--k')
    title('MASC: Gaze Time & Choice'); xlabel('Binned Number of Fixations on Option 1'); ylabel('p(Choice = Option 1)');
    res.(datasetName{iDat}).model.H3 = computeH3(nSubj,nTrials, attValues, w, choice,propFix1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% H4: Value Difference & Fixations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Data
    figure(3)
    subplot(2,4,(iDat*4)-3)
    [valFixCorr, vFC_ci,betasH4] = H4(nSubj,attValues,w,dataFix,matrixOAP,nTrials);
    plot(valFixCorr','Color',make_lighter([0 0 255],.8));
    hold on
    errorbar(1:20, mean(valFixCorr(:,1:20),'omitnan'),vFC_ci(1:20),'o-b','LineWidth',2,'CapSize',0,'MarkerFaceColor','b','MarkerEdgeColor','none');
    xlim([0 20]); ylim([-.5 .5]);
    xlabel('Fixation Number'); ylabel('Correlation Coefficient')
    title('Data: Fixations & Value Difference')
    res.(datasetName{iDat}).data.H4 = computeH4(betasH4);

    % Model
    subplot(2,4,(iDat*4)-1)
    [valFixCorr, vFC_ci,betasH4] = H4(nSubj,attValues,w,allFix,matrixOAP,nTrials);
    plot(mean(valFixCorr(:,1:20),'omitnan'),'b-o','MarkerFaceColor','b','MarkerEdgeColor','none','LineWidth',2);
    hold on
    plot(mean(valFixCorr(:,1:20),'omitnan')+vFC_ci(1:20),'b-');
    plot(mean(valFixCorr(:,1:20),'omitnan')-vFC_ci(1:20),'b-');
    xlim([0 20]); ylim([-.5 .5]);
    xlabel('Fixation Number'); ylabel('Correlation Coefficient')
    title('MASC: Fixations & Value Difference')
    res.(datasetName{iDat}).model.H4 = computeH4(betasH4);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% H5: Attraction Search Effect
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Data
    subplot(2,4,(iDat*4)-2)
    [ASE, ASE_ci] = H5(nSubj,nTrials,attValues,dataFix,matrixOAP);
    plot(1:2,ASE(:,1:2),'-','Color',make_lighter([.2,.2,.8]*255,.8),'HandleVisibility','off');
    hold on
    plot(3:4,ASE(:,3:4),'-','Color',make_lighter([.2,.8,.2]*255,.8),'HandleVisibility','off');
    errorbar(1:2,mean(ASE(:,1:2)),ASE_ci(1:2),'-o','Color',[.2,.2,.8],'LineWidth',2,'CapSize',0,'MarkerFaceColor',[.2,.2,.8],'MarkerEdgeColor','none');
    errorbar(3:4,mean(ASE(:,3:4)),ASE_ci(3:4),'-o','Color',[.2,.8,.2],'LineWidth',2,'CapSize',0,'MarkerFaceColor',[.2,.8,.2],'MarkerEdgeColor','none');
    set(gca,'XTick',1:4,'XTickLabel',{'Positive','Negative','Positive','Negative'})
    xlabel('Sign of 1st Attribute');ylabel('p(Transition)')
    legend({'Stay','Switch'}); legend boxoff
    title('Data: Attraction Search Effect');xlim([.75,4.25]);ylim([0,1])
    res.(datasetName{iDat}).data.H5 = computeH5(ASE);

    % Model
    subplot(2,4,(iDat*4))
    [ASE, ASE_ci] = H5(nSubj,nTrials,attValues,allFix,matrixOAP);
    plot(1:2,mean(ASE(:,1:2)),'.-','Color',[.2,.2,.8],'MarkerSize',20,'LineWidth',3);
    hold on
    plot(3:4,mean(ASE(:,3:4)),'.-','Color',[.2,.8,.2],'MarkerSize',20,'LineWidth',3);
    plot(1:2,mean(ASE(:,1:2)+ASE_ci(1:2)),'-','Color',[.2,.2,.8],'HandleVisibility','off')
    plot(1:2,mean(ASE(:,1:2)-ASE_ci(1:2)),'-','Color',[.2,.2,.8],'HandleVisibility','off')
    plot(3:4,mean(ASE(:,3:4)+ASE_ci(3:4)),'-','Color',[.2,.8,.2],'HandleVisibility','off')
    plot(3:4,mean(ASE(:,3:4)-ASE_ci(3:4)),'-','Color',[.2,.8,.2],'HandleVisibility','off')
    set(gca,'XTick',1:4,'XTickLabel',{'Positive','Negative','Positive','Negative'})
    xlabel('Sign of 1st Attribute');ylabel('p(Transition)')
    legend({'Stay','Switch'}); legend boxoff
    title('MASC: Attraction Search Effect');xlim([.75,4.25]);ylim([0,1])
    res.(datasetName{iDat}).model.H5 = computeH5(ASE);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% H6: Fixations Over Time by Weights
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(4)
    maxCol = [0, .25, 0];
    midCol = [0, .75, 0];
    minCol = [.25, 1, .25];

    % Data
    out = H6(w,nTrials,nSubj,dataFix,matrixOAP);
    subplot(2,2,(iDat*2)-1)
    plot(out.FixAttrMax_s',':','Color',make_lighter(maxCol*255,.5),'LineWidth',.25,'HandleVisibility','off');
    hold on
    plot(out.FixAttrMid_s',':','Color',make_lighter(midCol*255,.5),'LineWidth',.25,'HandleVisibility','off');
    plot(out.FixAttrMin_s',':','Color',make_lighter(minCol*255,.5),'LineWidth',.25,'HandleVisibility','off');
    errorbar(1:20, mean(out.FixAttrMax_s(:,1:20),'omitnan'),out.FixAttrMax_ci(1:20),'o-','Color',maxCol,'LineWidth',2,'CapSize',0,'MarkerFaceColor',maxCol,'MarkerEdgeColor',maxCol);
    errorbar(1:20, mean(out.FixAttrMid_s(:,1:20),'omitnan'),out.FixAttrMid_ci(1:20),'o-','Color',midCol,'LineWidth',2,'CapSize',0,'MarkerFaceColor',midCol,'MarkerEdgeColor',midCol);
    errorbar(1:20, mean(out.FixAttrMin_s(:,1:20),'omitnan'),out.FixAttrMin_ci(1:20),'o-','Color',minCol,'LineWidth',2,'CapSize',0,'MarkerFaceColor',minCol,'MarkerEdgeColor',minCol);
    xlim([0 20]);ylim([0 1])
    yline(1/mm,'--k')
    title('Data: Fixations & Attribute Weights'); xlabel('Fixation Number'); ylabel('p(Fixation)');
    legend({'Highest Weight','Medium Weight','Lowest Weight'},'Location','NorthEast');legend boxoff
    res.(datasetName{iDat}).data.H6 = computeH6(nSubj, nTrials,dataFix,w,matrixOAP);

    % Model
    out = H6(w,nTrials,nSubj,allFix,matrixOAP);
    subplot(2,2,(iDat*2))
    plot(mean(out.FixAttrMax_s(:,1:20),'omitnan'),'-o','Color',maxCol,'LineWidth',2,'MarkerFaceColor',maxCol,'MarkerEdgeColor','none');
    hold on
    plot(mean(out.FixAttrMax_s(:,1:20),'omitnan') + out.FixAttrMax_ci(1:20),'-','Color',maxCol,'LineWidth',2,'HandleVisibility','off');
    plot(mean(out.FixAttrMax_s(:,1:20),'omitnan') - out.FixAttrMax_ci(1:20),'-','Color',maxCol,'LineWidth',2,'HandleVisibility','off');

    plot(mean(out.FixAttrMid_s(:,1:20),'omitnan'),'-o','Color',midCol,'LineWidth',2,'MarkerFaceColor',midCol,'MarkerEdgeColor','none');
    plot(mean(out.FixAttrMid_s(:,1:20),'omitnan') + out.FixAttrMid_ci(1:20),'-','Color',midCol,'LineWidth',2,'HandleVisibility','off');
    plot(mean(out.FixAttrMid_s(:,1:20),'omitnan') - out.FixAttrMid_ci(1:20),'-','Color',midCol,'LineWidth',2,'HandleVisibility','off');

    plot(mean(out.FixAttrMin_s(:,1:20),'omitnan'),'-o','Color',minCol,'LineWidth',2,'MarkerFaceColor',minCol,'MarkerEdgeColor','none');
    plot(mean(out.FixAttrMin_s(:,1:20),'omitnan') + out.FixAttrMin_ci(1:20),'-','Color',minCol,'LineWidth',2,'HandleVisibility','off');
    plot(mean(out.FixAttrMin_s(:,1:20),'omitnan') - out.FixAttrMin_ci(1:20),'-','Color',minCol,'LineWidth',2,'HandleVisibility','off');
    xlim([0 20]);ylim([0 1])
    title('MASC: Fixations & Attribute Weights'); xlabel('Fixation Number'); ylabel('p(Fixation)');
    legend({'Highest Weight','Medium Weight','Lowest Weight'},'Location','NorthEast');legend boxoff
    res.(datasetName{iDat}).model.H6 = computeH6(nSubj, nTrials,allFix,w,matrixOAP);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% H7: Payne Index & Weight Difference
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(5)
    % Data
    [distWeight, PayneIndex, b] = H7(dataFix,nSubj,nTrials,matrixOAP,w);
    subplot(2,2,(iDat*2)-1)
    plot(distWeight,PayneIndex,'k.','MarkerSize',20);hold on
    plot(distWeight,distWeight*b(2)+b(1),'-','Color',[.2,.6,.6],'LineWidth',3)
    xlabel('1st Weight - 2nd Weight');ylabel('Payne Index')
    title('Data: Weights & Payne Index');xlim([0,1]);ylim([0,1])
    res.(datasetName{iDat}).data.H7 = computeH7(distWeight,PayneIndex);

    % Model
    [distWeight, PayneIndex, b] = H7(allFix,nSubj,nTrials,matrixOAP,w);
    subplot(2,2,(iDat*2))
    plot(distWeight,PayneIndex,'k.','MarkerSize',20); hold on
    plot(distWeight,distWeight*b(2)+b(1),'-','Color',[.2,.6,.6],'LineWidth',3)
    xlabel('1st Weight - 2nd Weight');ylabel('Payne Index')
    title('MASC: Weights & Payne Index');xlim([0,1]);ylim([0,1])
    res.(datasetName{iDat}).model.H7 = computeH7(distWeight,PayneIndex);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Adding Icons, Labels & Figure Printing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% names = {'RT_CC','GazeTime_Choice','ASE','FixationWeights','PayneIndex'};
% for i = 1:5
%     figure(i)
%     addIcons(phoneIcon,hotelIcon, phonePosition,hotelPosition)
%     set(gcf, 'WindowState', 'maximized');
%     print(gcf,'-dpng','-r300',[pwd,'/CurrentPaper/Figures/',names{i}])
%     savefig(gcf,[pwd,'/Figures/',names{i}]);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%












%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions Used Throughout
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% H1 & H2: Choice Consistency and Number Of Fixations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [CC_l, CC_ci, RT_l, RT_ci] = H1(attValues,w,nTrials,nSubj,difficulty,choice,RT,n)
        consistentChoice = zeros(nTrials,nSubj);
        for iSub = 1:nSubj
            optValues = zeros(n,nTrials); %option values
            for iTrial = 1:nTrials
                optValues(:,iTrial) = attValues(:,:,iTrial,iSub)*w(:,iSub);
            end
            consistentChoice(:,iSub) = choice(:,iSub)==((1:n)*(repmat(max(optValues),n,1)==optValues))'; %whether the option with the highest value was chosen
        end

        %difficulty levels
        CC_l = [(sum(consistentChoice.*(difficulty==1))./sum(difficulty==1))',...
            (sum(consistentChoice.*(difficulty==2))./sum(difficulty==2))',...
            (sum(consistentChoice.*(difficulty==3))./sum(difficulty==3))'];
        CC_ci = norminv(.975)*std(CC_l)./sqrt(nSubj);
        RT_l = [(sum(RT.*(difficulty==1),'omitnan')./sum(difficulty==1))',...
            (sum(RT.*(difficulty==2),'omitnan')./sum(difficulty==2))',...
            (sum(RT.*(difficulty==3),'omitnan')./sum(difficulty==3))'];
        RT_ci = norminv(.975)*std(RT_l)./sqrt(nSubj);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% H3: Dwell Time & Choice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [catPropFix, ciPropFix,ncat,propFix1] = H3(allFix, choice, matrixOAP,nSubj,nTrials)
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
ciPropFix = norminv(.975)*std(catPropFix)./sqrt(length(catPropFix));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% H4: Value Differences & Fixations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [valFixCorr, vFC_ci,betas] = H4(nSubj,attValues,w,allFix,matrixOAP,nTrials)
betas = zeros(nSubj,4);
valFixCorr = zeros(nSubj,100)+nan;
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
    for ff = 1:maxFixGLM
        valFixCorr(s,ff) = corr(glmValDiff(ff==glmFixNum),glmFix1(ff==glmFixNum));
    end

    %stats
    glmFixNum = (glmFixNum-mean(glmFixNum))./std(glmFixNum);
    glmValDiff = (glmValDiff-mean(glmValDiff))./std(glmValDiff);
    glmInteract = glmFixNum.*glmValDiff;
    betas(s,:) = glmfit([glmFixNum,glmValDiff,glmInteract],glmFix1,'binomial');
end
vFC_ci = norminv(.975)*nanstd(valFixCorr)./sqrt(nSubj);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% H5: Attraction Search Effect
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ASE, ASE_ci] = H5(nSubj,nTrials,attValues,allFixH5,matrixOAP)

firstAttValuePos = nan+zeros(nSubj,nTrials);
stayFix = nan+zeros(nSubj,nTrials);
for s = 1:nSubj
    for t = 1:nTrials
        if ~isnan(allFixH5(2,t,s)) %you need at least the first two fixations
            currentAttValues = attValues(:,:,t,s);
            firstAttValuePos(s,t) = currentAttValues(allFixH5(1,t,s))>0; %is first attribute value positive
            stayFix(s,t) = min(sum(allFixH5(1,t,s)==matrixOAP,2)==sum(allFixH5(2,t,s)==matrixOAP,2))==1; %does 2nd fixation remain on same option?
        end
    end
end

ASE = [sum((stayFix==1).*(firstAttValuePos==1),2)./sum((firstAttValuePos==1),2),...
    sum((stayFix==1).*(firstAttValuePos==0),2)./sum((firstAttValuePos==0),2),...
    sum((stayFix==0).*(firstAttValuePos==1),2)./sum((firstAttValuePos==1),2),...
    sum((stayFix==0).*(firstAttValuePos==0),2)./sum((firstAttValuePos==0),2)];

ASE_ci = norminv(.975)*std(ASE)./sqrt(nSubj);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% H6: Fixations over Time by Weights 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = H6(w,nTrials,nSubj,allFix,matrixOAP)
[~,mmH6] = size(matrixOAP);
maxSteps = 100;
fixAttWeight = zeros(mmH6,maxSteps,nTrials,nSubj);
for s = 1:nSubj
    [~,attRank] = sort(w(:,s),'descend'); %the ith entry of attRank is the ith-important attribute
    for t = 1:nTrials
        for j = 1:mmH6 %loop over attributes (from most to least important)
            fixAttWeight(j,:,t,s) = ismember(allFix(:,t,s),matrixOAP(:,attRank(j)))...
                ./~isnan(allFix(:,t,s));
        end
    end
end
fixAttWeight = reshape(mean(fixAttWeight,3,'omitnan'),mmH6,maxSteps,nSubj);

out.FixAttrMax_s = reshape(fixAttWeight(1,:,:),maxSteps,nSubj)';
out.FixAttrMid_s = reshape(fixAttWeight(2,:,:),maxSteps,nSubj)';
out.FixAttrMin_s = reshape(fixAttWeight(3,:,:),maxSteps,nSubj)';
out.FixAttrMax_ci = norminv(.975)*std(out.FixAttrMax_s,'omitnan')./sqrt(nSubj);
out.FixAttrMid_ci = norminv(.975)*std(out.FixAttrMid_s,'omitnan')./sqrt(nSubj);
out.FixAttrMin_ci = norminv(.975)*std(out.FixAttrMin_s,'omitnan')./sqrt(nSubj);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% H7: Payne Index & Weight Difference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [distWeight, PayneIndex, b] = H7(allFix,nSubj,nTrials,matrixOAP,w)
transAlt = nan+zeros(size(allFix));
transAtt = nan+zeros(size(allFix));
for s = 1:nSubj
    for t = 1:nTrials
        if sum(~isnan(allFix(:,t,s)))>1 %you need to have at least 2 fixations
            for ff = 1:sum(isnan(allFix(:,t,s)))-1
                transAlt(ff,t,s) = max(sum((matrixOAP==allFix(ff,t,s))+(matrixOAP==allFix(ff+1,t,s)),2))==2;
                transAtt(ff,t,s) = (max(sum((matrixOAP==allFix(ff,t,s))+(matrixOAP==allFix(ff+1,t,s))))==2)&(transAlt(ff,t,s)~=1);
            end
        end
    end
end

aw = abs(w);
distWeight = (max(aw)-max((aw~=max(aw)).*aw))'; %difference of highest and second-highest weight per subject
PayneIndex = reshape((nansum(nansum(transAlt))-nansum(nansum(transAtt))),nSubj,1)./...
    reshape((nansum(nansum(transAlt))+nansum(nansum(transAtt))),nSubj,1); %standard Payne Index
b = glmfit(distWeight,PayneIndex);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add Icons & Labels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function addIcons(phoneIcon,hotelIcon, phonePosition,hotelPosition)
axes('Position',phonePosition);
imshow(phoneIcon);
axis off;
axes('Position',hotelPosition);
imshow(hotelIcon);
axis off;

annotation(gcf,'textbox',...
    [0.129645833333333 0.950697674418605 0.0213958333333333 0.0437209302325584],...
    'String',{'A'},...
    'FitBoxToText','off','FontSize',32,'EdgeColor','none','FontWeight','bold')
annotation(gcf,'textbox',...
    [0.54266666 0.950697674418605 0.0213958333333333 0.0437209302325584],...
    'String',{'B'},...
    'FitBoxToText','off','FontSize',32,'EdgeColor','none','FontWeight','bold')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulate Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [choice,RT, allFix] = simulateMASC(f)
 [choice, RT, allFix] = switchModelCall(f.settings,f.parameters,f.attValues,f.nTrials,f.nSubj,f.model);


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Best Parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions for computing stats.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% H1 & H2: Choice Consistency & Number of Fixations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function [statsH1, statsH2] = computeH1_H2(CC_l,RT_l)
statsH1 = table(); statsH2 = table();
namesH1H2 = {'easy-med','easy-hard','med-hard'};
comps = [1 2; 1 3; 2 3];
for ii = 1:3
    [~,p12,~,stats] = ttest(CC_l(:,comps(ii,1)),CC_l(:,comps(ii,2)),'tail','right');
    statsH1.t(ii) = stats.tstat;
    statsH1.df(ii) = stats.df;
    statsH1.p(ii) = p12;
    statsH1.d(ii) = computeCohen_d(CC_l(:,comps(ii,1)), CC_l(:,comps(ii,2)),'paired');


    [ ~,p12,~,stats] = ttest(RT_l(:,comps(ii,1)),RT_l(:,comps(ii,2)),'tail','left');
    statsH2.t(ii) = stats.tstat;
    statsH2.df(ii) = stats.df;
    statsH2.p(ii) = p12;
    statsH2.d(ii) = computeCohen_d(RT_l(:,comps(ii,1)), RT_l(:,comps(ii,2)),'paired');
end

statsH1.Properties.RowNames = namesH1H2;
statsH2.Properties.RowNames = namesH1H2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% H3: Dwell Time & Choice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function statsH3 = computeH3(nSubj,nTrials, attValues, w, choice,propFix1)
betasH3 = zeros(nSubj,3);
for s = 1:nSubj
    valDiff = zeros(nTrials,1); %control for value difference in this test
    for t = 1:nTrials
        optValues = attValues(:,:,t,s)*w(:,s);
        valDiff(t) = optValues(1)-mean(optValues(2:end));
    end
    warning('off')
    lastwarn('')
    betasH3(s,:) = glmfit([propFix1(s,:)',valDiff],choice(:,s)==1,'binomial');%,'LikelihoodPenalty','jeffreys-prior');
    %warningMessage = lastwarn;
    % if ~isempty(warningMessage)
    %betas(s,:) = nan;
    %  bad(s) = 1;
    %  else
    %      bad(s) = 0;
    % end
end
statsH3 = table();
[~,p,~,stats] = ttest(betasH3(:,2),zeros(nSubj,1),'tail','right');
statsH3.t = stats.tstat;
statsH3.df = stats.df;
statsH3.p = p;
statsH3.d = computeCohen_d(betasH3(:,2),zeros(nSubj,1),'paired');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% H4: Fixations & Value Difference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function statsH4 = computeH4(betas)
namesH4 = {'Intercept','FixNum','ValueDiff','Interact'};

[~,p1,~,s1] = ttest(betas,zeros(size(betas)),'tail','right');
d1 = mean(betas)./std(betas);

statsH4 = table();
for ii = 1:4
    statsH4.tstat(ii) = s1.tstat(ii);
    statsH4.df(ii) = s1.df(ii);
    statsH4.p(ii) = p1(ii);
    statsH4.d(ii) = d1(ii);
end
statsH4.Properties.RowNames = namesH4;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% H5: Attraction Search Effect
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function statsH5 = computeH5(ASE)
statsH5 = table();
[~,p,~,stats] = ttest(ASE(:,1),ASE(:,2),'tail','right');
statsH5.t = stats.tstat;
statsH5.df = stats.df;
statsH5.p = p;
statsH5.d = computeCohen_d(ASE(:,1),ASE(:,2),'paired');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% H6: Fixations over Time by Weights
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function statsH6 = computeH6(nSubj, nTrials,allFix,w,matrixOAP)
%stats
[~,mH6] = size(matrixOAP);
maxFix = 20;
betasH6 = zeros(nSubj,2);
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
    betasH6(s,:) = glmfit(zFixNum,glmHighestFix(glmFixNum>0),'binomial');
    betasH6(s,:) = glmfit(glmFixNum(glmFixNum>0)-1,glmHighestFix(glmFixNum>0),'binomial');
end
chanceLevelIntercept = log((1/mH6)/(1-(1/mH6))); %chance level of the logistic intercept for m attributes

namesH6 = {'Intercept','Coefficient'};
statsH6 = table();
[~,p1,~,stats] = ttest(betasH6(:,1),zeros(nSubj,1)+chanceLevelIntercept,'tail','right');
statsH6.t(1) = stats.tstat;
statsH6.df(1) = stats.df;
statsH6.p(1) = p1;
statsH6.d(1) = computeCohen_d(betasH6(:,1),zeros(nSubj,1)+chanceLevelIntercept,'paired');

[~,p2,~,stats] = ttest(betasH6(:,2),zeros(nSubj,1),'tail','left');
statsH6.t(2) = stats.tstat;
statsH6.df(2) = stats.df;
statsH6.p(2) = p2;
statsH6.d(2) = computeCohen_d(betasH6(:,2),zeros(nSubj,1),'paired');

statsH6.Properties.RowNames = namesH6;
%significantEffects.test6 = (p1<.05)&(p2<.05);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% H7: Payne Index & Weight Difference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function statsH7 = computeH7(distWeight,PayneIndex)
statsH7 = table();
[r,p] = corr(distWeight(~isnan(PayneIndex)),PayneIndex(~isnan(PayneIndex)),'tail','left');
statsH7.r = r;
statsH7.p = p;
%significantEffects.test7 = p<.05;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute Cohen's d
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = computeCohen_d(x1, x2, varargin)
%
% call: d = computeCohen_d(x1, x2, varargin)
%
% EFFECT SIZE of the difference between the two
% means of two samples, x1 and x2 (that are vectors),
% computed as "Cohen's d".
%
% If x1 and x2 can be either two independent or paired
% samples, and should be treated accordingly:
%
%   d = computeCohen_d(x1, x2, 'independent');  [default]
%   d = computeCohen_d(x1, x2, 'paired');
%
% Note: according to Cohen and Sawilowsky:
%
%      d = 0.01  --> very small effect size
%      d = 0.20  --> small effect size
%      d = 0.50  --> medium effect size
%      d = 0.80  --> large effect size
%      d = 1.20  --> very large effect size
%      d = 2.00  --> huge effect size
%
%
% Ruggero G. Bettinardi (RGB)
% Cellular & System Neurobiology, CRG
% -------------------------------------------------------------------------------------------
%
% Code History:
%
% 25 Jan 2017, RGB: Function is created

if nargin < 3, testType = 'independent';
else           testType = varargin{1};
end

% basic quantities:
n1       = numel(x1);
n2       = numel(x2);
mean_x1  = nanmean(x1);
mean_x2  = nanmean(x2);
var_x1   = nanvar(x1);
var_x2   = nanvar(x2);
meanDiff = (mean_x1 - mean_x2);

% select type of test:
isIndependent = strcmp(testType, 'independent');
isPaired      = strcmp(testType, 'paired');

% compute 'd' accordingly:
if isIndependent

    sv1      = ((n1-1)*var_x1);
    sv2      = ((n2-1)*var_x2);
    numer    =  sv1 + sv2;
    denom    = (n1 + n2 - 2);
    pooledSD =  sqrt(numer / denom); % pooled Standard Deviation
    s        = pooledSD;             % re-name
    d        =  meanDiff / s;        % Cohen's d (for independent samples)

elseif isPaired

    haveNotSameLength = ~isequal( numel(x1), numel(x2) );
    if haveNotSameLength, error('In a paired test, x1 and x2 have to be of same length!'), end

    deltas   = x1 - x2;         % differences
    sdDeltas = nanstd(deltas);  % standard deviation of the diffferences
    s        = sdDeltas;        % re-name
    d        =  meanDiff / s;   % Cohen's d (paired version)

end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



end

