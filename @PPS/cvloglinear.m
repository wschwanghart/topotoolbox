function [ypredcv,stats] = cvloglinear(P,mdl,varargin)

%CVLOGLINEAR Cross-validate a loglinear point process model
%
% Syntax
%
%     ypredcv = cvloglinear(P,mdl)
%     ypredcv = cvloglinear(P,mdl,'pn',pv,...)
%     [ypredcv,stats] = ...
%
% Description
%
%     cvloglinear cross-validates a loglinear model obtained with
%     PPS/fitloglinear. The function subdivides the stream network stored
%     in the PPS object in segments with more or less equal length (note
%     that the network is always split at confluences and thus shorter
%     reaches likely exist). Subsequently, the segments are randomly
%     partitioned into five groups, each of which has approximately the 
%     same number of observations. The first group is then used as test data 
%     while the other groups are used as training data. This step is
%     repeated for each group.
%
% Input arguments
%
%     P           PPS object
%     mdl         model fit with fitloglinear
%     seglength   length of segments into which the stream network
%                 underlying P is partitioned.
%     kfold       number of disjoint subsamples in the k-fold cross-
%                 validation. Default is 5.
%
% Output arguments
%
%     ypredcv     node-attribute list of cross-validated predictions of
%                 point density.
%     stats       AUC values
%
% Example
%
%
% 
% See also: PPS, PPS/fitloglinear
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 24. April, 2023


p = inputParser;
p.FunctionName = 'cvloglinear';

addParameter(p,'seglength',500*P.S.cellsize);
addParameter(p,'kfold',5)
addParameter(p,'repeat',1)
parse(p,varargin{:});


% Create segments
CS = STREAMobj2cell(P.S,'segments',p.Results.seglength);
nseg = numel(CS);

% Preallocate output
ypredcv = repmat(getnal(P.S),1,p.Results.repeat);

for iter = 1:p.Results.repeat

% partition object
c = cvpartition(nseg,'KFold',p.Results.kfold);

for r = 1:p.Results.kfold
    
    idxtrain = training(c,r);
    idxtest  = test(c,r);
    
    if nnz(idxtest) == 1
        STEST = CS{idxtest};
    else
        STEST  = union(CS{idxtest});
    end

    if nnz(idxtrain) == 1
        STRAIN = CS{idxtrain};
    else
        STRAIN = union(CS{idxtrain});
    end
    
    warning off
    PTEST  = PPS(STEST,'pp',P.S.IXgrid(P.PP));
    PTRAIN = PPS(STRAIN,'pp',P.S.IXgrid(P.PP));
    warning on
    nvar   = size(mdl.Variables,2)-1;
    covariates_test  = repmat(getnal(PTEST.S),1,nvar);
    covariates_train = repmat(getnal(PTRAIN.S),1,nvar);

    for r2 = 1:nvar
        covariates_test(:,r2)  = nal2nal(PTEST.S, P.S,table2array(mdl.Variables(:,r2)));
        covariates_train(:,r2) = nal2nal(PTRAIN.S,P.S,table2array(mdl.Variables(:,r2)));
    end

    mdltrain = fitloglinear(PTRAIN,covariates_train,'model',mdl.Formula);
    
    % Use trained model to predict
    ypredtest = predict(mdltrain,covariates_test);
    ypredcv(:,iter) = nal2nal(P.S,PTEST.S,ypredtest,ypredcv(:,iter));

end
end
    
d   = distance(P.S,'node_to_node');
d   = mean(d);
ypredcv = ypredcv./d; %.cellsize;

for iter = 1:p.Results.repeat
    [~,~,stats.auc(iter)] = roc(P,ypredcv(:,iter),'plot',false,'perfcurve',false);
end
    


