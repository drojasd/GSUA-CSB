function [T,clusterInfo] = gsua_dia(T,T_est,varargin)
%GSUA_DIA Practical identifiability diagnostics (headless, no plots)
%
%   T = gsua_dia(T,T_est)
%   T = gsua_dia(T,T_est,outlier)
%   [T,clusterInfo] = gsua_dia(T,T_est,outlier,cluster)
%   [T,clusterInfo] = gsua_dia(T,T_est,outlier,cluster,maxK,silThreshold)
%   [T,clusterInfo] = gsua_dia(...,'cost',cost,'CostMethod','gap')
%
%   Computes an updated confidence range and an identifiability index
%   for each parameter from a set of repeated estimations, without
%   producing any figures. This is the non-plotting counterpart of
%   GSUA_IA: same range/index computation (parametric normal CI on the
%   median), useful for scripted/batch pipelines where the diagnostic
%   plots of GSUA_IA are not needed.
%
%   Multiple global minima (cluster = true):
%   Same spectral-clustering detection as GSUA_IA, without plots: the
%   normalized estimates are clustered for candidate counts 2:maxK using
%   spectralcluster, the split with the best mean silhouette score
%   (evalclusters) is kept only if that score is >= silThreshold,
%   otherwise the estimates are treated as one basin. Requires Statistics
%   and Machine Learning Toolbox. See GSUA_IA for the full explanation.
%
%   Fit-quality filtering ('cost'):
%   A repeated-estimation run that converged to a much worse cost than
%   the best one contributes essentially arbitrary parameter values --
%   left in, it can make a genuinely well-identified parameter look
%   poorly identified purely because an optimizer run failed. Passing
%   'cost' (e.g. the res output of gsua_pe, one entry per column of
%   T_est) screens those runs out before anything else is computed. See
%   GSUA_COSTCUTOFF for the two available methods ('rtol'/'gap') and the
%   minKeep safety floor.
%
%   Clustering runs BEFORE outlier removal, on purpose: global
%   Mahalanobis-distance outlier detection assumes one unimodal
%   population, so running it first would treat a genuine second basin
%   as "outliers" relative to the pooled mean/covariance across both
%   basins and could delete it before clustering ever finds it. When a
%   genuine multi-cluster split is found, T.Range/T.index and the
%   correlation matrix are computed from the dominant (largest)
%   cluster's points only -- pooling separated basins into one interval
%   or correlation isn't a meaningful summary -- and outlier removal, if
%   requested, is scoped to that single cluster, where the unimodal
%   assumption behind Mahalanobis distance is actually valid.
%
%   Note on outlier=true: the bundled DetectMultVarOutliers utility
%   (Add Funcs/AntonSemechko-Multivariate-Outliers-.../) defaults to
%   assuming up to roughly half the sample could be contaminated, which
%   makes it considerably more aggressive than a plain Mahalanobis cutoff
%   on data that is actually clean and tightly clustered -- confirmed to
%   remove the large majority of a low-noise synthetic single cluster in
%   testing, independent of this function's cluster/outlier ordering.
%   This is a pre-existing characteristic of that third-party utility,
%   not something introduced here; inspect T.Est after the call if
%   outlier=true results look surprising.
%
%   Inputs:
%     T            <-- summary table from gsua_dataprep
%     T_est        <-- Np x N matrix, one column per repeated parameter
%                      estimation (e.g. Table.Estlsqc from gsua_pe)
%     outlier      <-- (optional, logical) remove multivariate outliers
%                      before computing the index. Default: false
%     cluster      <-- (optional, logical) run spectral clustering to
%                      detect multiple global minima. Default: false
%     maxK         <-- (optional) largest cluster count to test. Default: 5
%     silThreshold <-- (optional) minimum mean silhouette score required
%                      to accept a multi-cluster split. Default: 0.5
%     'cost'         <-- (optional, paired) 1xN fit cost, one per column
%                        of T_est. Default: [] (no filtering)
%     'CostMethod'   <-- (optional, paired) 'rtol' (default) or 'gap'
%     'CostRtol'     <-- (optional, paired) Default: 0.1
%     'CostAtol'     <-- (optional, paired) Default: 1e-8
%     'CostGapRatio' <-- (optional, paired) Default: 3
%     'MinKeep'      <-- (optional, paired) Default: 3
%
%   Outputs:
%     T           <-- same table with T.Range replaced by the new
%                     confidence interval, T.Nominal set to the median
%                     estimate, T.index set to the per-parameter
%                     identifiability index (0 = well identified, 1 =
%                     poorly identified), and T.Est set to the estimates
%                     actually used to compute the above (after cost
%                     filtering, cluster restriction, and outlier
%                     removal have all been applied -- may be a subset
%                     of the input T_est)
%     clusterInfo <-- struct with fields NumClusters, Idx, Centers
%                     (candidate global minima, one row per cluster),
%                     ClusterSize, Silhouette, and Dominant (index of the
%                     largest cluster, the one T.Range/T.index reflect
%                     when NumClusters > 1). Single-basin default
%                     (NumClusters = 1, Dominant = 1) when cluster is
%                     false or no split beat silThreshold.
%
%   See also GSUA_IA, GSUA_PE, GSUA_MEDIANCI, GSUA_COSTCUTOFF.

p = inputParser;
addRequired(p,'T');
addRequired(p,'T_est');
addOptional(p,'outlier',false);
addOptional(p,'cluster',false);
addOptional(p,'maxK',5);
addOptional(p,'silThreshold',0.5);
addParameter(p,'cost',[]);
addParameter(p,'CostMethod','rtol');
addParameter(p,'CostRtol',0.1);
addParameter(p,'CostAtol',1e-8);
addParameter(p,'CostGapRatio',3);
addParameter(p,'MinKeep',3);
parse(p,T,T_est,varargin{:});

outlier=p.Results.outlier;
cluster=p.Results.cluster;
maxK=p.Results.maxK;
silThreshold=p.Results.silThreshold;
cost=p.Results.cost;

if ~isempty(cost)
    keepCost=gsua_costcutoff(cost,p.Results.CostMethod,p.Results.CostRtol, ...
        p.Results.CostAtol,p.Results.CostGapRatio,p.Results.MinKeep);
    nBadFit=sum(~keepCost);
    if nBadFit>0
        disp([num2str(nBadFit) ' run(s) dropped by the fit-quality filter (converged to a much worse cost than the best one)'])
    end
    T_est=T_est(:,keepCost);
end

Normalized=zeros(size(T,1),size(T_est,2));
for i=1:size(T,1)
    Normalized(i,:)=(T_est(i,:)-T.Range(i,1))/(T.Range(i,2)-T.Range(i,1));
end

clusterInfo=struct('NumClusters',1,'Idx',ones(1,size(T_est,2)),'Centers',median(T_est,2)', ...
    'ClusterSize',size(T_est,2),'Silhouette',[],'Dominant',1);
statT_est=T_est; % point set stats below are computed from; narrows to the dominant cluster if found
if cluster
    Nruns=size(Normalized,2);
    X=Normalized';
    kmax=min(maxK,Nruns-1);
    if kmax>=2
        try
            eva=evalclusters(X,@(x,k) spectralcluster(x,k),'silhouette','KList',2:kmax);
            bestK=eva.OptimalK;
            bestSil=max(eva.CriterionValues,[],'omitnan');
        catch ME
            warning(ME.message)
            bestK=1;
            bestSil=0;
        end
        if bestK>1 && bestSil>=silThreshold
            cIdx=spectralcluster(X,bestK)';
            centers=zeros(bestK,size(T_est,1));
            clusterSize=zeros(1,bestK);
            for k=1:bestK
                centers(k,:)=median(T_est(:,cIdx==k),2)';
                clusterSize(k)=sum(cIdx==k);
            end
            [~,dominant]=max(clusterSize);
            clusterInfo=struct('NumClusters',bestK,'Idx',cIdx,'Centers',centers, ...
                'ClusterSize',clusterSize,'Silhouette',bestSil,'Dominant',dominant);
            disp(['Detected ' num2str(bestK) ' candidate global minima (mean silhouette ' num2str(bestSil,3) ')'])
            disp(['Cluster sizes (of ' num2str(Nruns) ' runs): ' num2str(clusterSize)])
            statT_est=T_est(:,cIdx==dominant);
        else
            disp(['No confident multi-cluster split found (best mean silhouette ' num2str(bestSil,3) ' < ' num2str(silThreshold) ')'])
            disp('Treating estimations as a single basin (1 candidate global minimum)')
        end
    else
        warning('Not enough estimation runs to test clustering (need at least 3 for maxK>=2)')
    end
end

if outlier
    disp('Removing outliers...')
    [~,~,RD,chi_crt]=DetectMultVarOutliers(statT_est');
    id_in=RD<chi_crt(4);
    statT_est=statT_est(:,id_in);
    disp(num2str(sum(id_in))+" outliers were removed")
end

RHO = corr(statT_est');

x=statT_est;
med=median(x,2);
N=size(x,2);
desv=std(x,[],2);

lb = med-1.96*sqrt(pi/2)*desv/sqrt(N);
up = med+1.96*sqrt(pi/2)*desv/sqrt(N);

if (any(lb(lb<T.Range(:,1))))||(any(up(up>T.Range(:,2))))
    lb(lb<T.Range(:,1))= T.Range(lb<T.Range(:,1),1);
    up(up>T.Range(:,2))= T.Range(up>T.Range(:,2),2);
end


boxin=(up-lb)./(T.Range(:,2)-T.Range(:,1));
len=length(boxin);
corrin=sum(abs(RHO))/len;
extrin=sum(abs(RHO)>0.5)/len;
ind=(2*boxin'+corrin+extrin)/4;

T.Est=statT_est; % the point set the stats above were actually computed from (post cost/cluster/outlier filtering)
T.Range=[lb,up];
T.index=ind';
T.Nominal=med;

if cluster
    try
        T=addprop(T,{'ClusterInfo'},{'table'});
    catch
    end
    try
        T.Properties.CustomProperties.ClusterInfo=clusterInfo;
    catch
    end
end

end
