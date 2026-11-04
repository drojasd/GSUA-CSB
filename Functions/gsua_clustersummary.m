function S = gsua_clustersummary(T,clusterInfo)
%GSUA_CLUSTERSUMMARY Per-cluster, per-parameter summary with an identifiability verdict
%
%   S = gsua_clustersummary(T,clusterInfo)
%
%   Builds a table with one row per (cluster, parameter) from the clusterInfo returned by GSUA_IA
%   or GSUA_DIA. For each parameter, in each cluster, it reports where the group CENTER sits and how
%   WIDE the plausible range is, and a suggested identifiability verdict -- the summary used in the
%   epidemiology hands-on of the Semillero course.
%
%   Two kinds of interval are reported, and they answer different questions:
%     * MedianCI / MeanCI describe the CENTER of the group (how well pinned the central value is):
%       the distribution-free confidence interval of the median (GSUA_MEDIANCI) and the Student-t
%       confidence interval of the mean. They shrink as the cluster grows.
%     * The 2.5-97.5 PERCENTILE interval [P2_5, P97_5] describes the SPREAD of plausible parameter
%       values across the ensemble. Its width relative to the parameter's ORIGINAL range
%       (PctWidthRel) is what the verdict is based on -- it does NOT shrink just because more
%       estimates were collected.
%
%   Verdict (on PctWidthRel): < 0.25 -> 'identifiable'; > 0.50 -> 'not identifiable'; in between
%   -> 'dubious'.
%
%   Inputs:
%     T           <-- summary table (supplies the parameter names via its RowNames)
%     clusterInfo <-- second output of GSUA_IA/GSUA_DIA (uses .NumClusters, .Est, .OrigRange)
%
%   Output:
%     S <-- table with columns Cluster, Parameter, n, Median, MedianCI_low, MedianCI_high, Mean,
%           MeanCI_low, MeanCI_high, P2_5, P97_5, PctWidthRel, Verdict
%
%   Example:
%     [T,ci] = gsua_ia(T,Est,false,false,true,true);
%     S = gsua_clustersummary(T,ci);
%     disp(S)
%
%   See also GSUA_IA, GSUA_DIA, GSUA_MEDIANCI, GSUA_BOOTSTRAP.

if ~isfield(clusterInfo,'OrigRange') || ~isfield(clusterInfo,'Est')
    error('gsua_clustersummary:BadClusterInfo', ...
        'clusterInfo must come from GSUA_IA/GSUA_DIA (needs the .Est and .OrigRange fields).');
end

names = T.Properties.RowNames;
origRange = clusterInfo.OrigRange;
rangeW = origRange(:,2)-origRange(:,1);
Np = numel(names);
Kc = clusterInfo.NumClusters;

Cluster=[]; Parameter={}; nCol=[]; Median=[]; MedianCI_low=[]; MedianCI_high=[];
Mean=[]; MeanCI_low=[]; MeanCI_high=[]; P2_5=[]; P97_5=[]; PctWidthRel=[]; Verdict={};

for c=1:Kc
    est = clusterInfo.Est{c};          % Np x n_c
    nC = size(est,2);
    [mlo,mhi] = gsua_medianCI(est,0.05);
    for i=1:Np
        e = est(i,:);
        med = median(e);
        mn = mean(e);
        sd = std(e);
        if nC>1
            moe = tinv(0.975,nC-1)*sd/sqrt(nC);
        else
            moe = 0;
        end
        pr = prctile(e,[2.5 97.5]);
        if rangeW(i)>0
            pw = (pr(2)-pr(1))/rangeW(i);
        else
            pw = 0;                     % fixed parameter: no spread to assess
        end
        if pw<0.25
            v = 'identifiable';
        elseif pw>0.5
            v = 'not identifiable';
        else
            v = 'dubious';
        end
        Cluster(end+1,1)=c; %#ok<*AGROW>
        Parameter{end+1,1}=names{i};
        nCol(end+1,1)=nC;
        Median(end+1,1)=med;
        MedianCI_low(end+1,1)=mlo(i);
        MedianCI_high(end+1,1)=mhi(i);
        Mean(end+1,1)=mn;
        MeanCI_low(end+1,1)=mn-moe;
        MeanCI_high(end+1,1)=mn+moe;
        P2_5(end+1,1)=pr(1);
        P97_5(end+1,1)=pr(2);
        PctWidthRel(end+1,1)=pw;
        Verdict{end+1,1}=v;
    end
end

S = table(Cluster,Parameter,nCol,Median,MedianCI_low,MedianCI_high,Mean,MeanCI_low,MeanCI_high, ...
    P2_5,P97_5,PctWidthRel,Verdict, ...
    'VariableNames',{'Cluster','Parameter','n','Median','MedianCI_low','MedianCI_high', ...
    'Mean','MeanCI_low','MeanCI_high','P2_5','P97_5','PctWidthRel','Verdict'});
end
