function [T,clusterInfo] = gsua_ia(T,T_est,varargin)
%GSUA_IA Practical identifiability analysis with diagnostic plots
%
%   T = gsua_ia(T,T_est)
%   T = gsua_ia(T,T_est,correction,outlier,show)
%   [T,clusterInfo] = gsua_ia(T,T_est,correction,outlier,show,cluster)
%   [T,clusterInfo] = gsua_ia(T,T_est,correction,outlier,show,cluster,maxK,silThreshold)
%   [T,clusterInfo] = gsua_ia(...,'cost',cost,'CostMethod','gap')
%
%   Analyzes practical parameter identifiability from a set of repeated
%   estimations: normalized CDFs and boxplots per parameter, a
%   parameter-correlation heatmap, and (optionally) an identifiability
%   index scatter/graph highlighting strongly correlated parameter
%   pairs. Updates T.Range to a distribution-free confidence interval
%   of the median (see GSUA_MEDIANCI) and reports an identifiability
%   index per parameter combining interval width, mean correlation, and
%   strong-correlation count. This is the plotting counterpart of
%   GSUA_DIA.
%
%   Multiple global minima (cluster = true):
%   Repeated estimations from a non-convex problem can land in distinct
%   basins of attraction ("multiple global minima") instead of scattering
%   around one point. When cluster is true, gsua_ia runs spectral
%   clustering (spectralcluster) on the normalized estimates for a range
%   of candidate cluster counts 2:maxK, and picks the count with the best
%   mean silhouette score (evalclusters). If the best mean silhouette is
%   below silThreshold, the estimates are treated as a single basin
%   (NumClusters = 1) — a low silhouette score means the "best" split
%   found is not actually well separated, so the sample is more likely
%   unimodal than confidently multimodal. Requires Statistics and Machine
%   Learning Toolbox.
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
%   genuine multi-cluster split is found, every plot and statistic below
%   (ECDFs, boxplots, the correlation heatmap, T.Range/T.index) reflects
%   the dominant (largest) cluster's points only -- pooling separated
%   basins into one interval or correlation isn't meaningful -- except
%   the "Candidate global minima" plot, which by design overlays all
%   clusters on the full (cost-filtered) sample. Outlier removal, if
%   requested, is scoped to the dominant cluster only, where the
%   unimodal assumption behind Mahalanobis distance is actually valid.
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
%     correction   <-- (optional, logical) clip the new range to the
%                      original T.Range when it would otherwise expand
%                      beyond it. Default: false
%     outlier      <-- (optional, logical) remove multivariate outliers
%                      before computing statistics. Default: false
%     show         <-- (optional, logical) plot the identifiability index
%                      scatter, correlation graph, and (if cluster is
%                      true) the clustering diagnostic plots. Default: false
%     cluster      <-- (optional, logical) run spectral clustering to
%                      detect multiple global minima among the repeated
%                      estimations. Default: false
%     maxK         <-- (optional) largest cluster count to test. Default: 5
%     silThreshold <-- (optional) minimum mean silhouette score required
%                      to accept a multi-cluster split over a single
%                      basin. Default: 0.5
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
%     clusterInfo <-- struct with fields NumClusters, Idx (1xN cluster
%                     label per estimation run), Centers (NumClusters x
%                     Np candidate global minima, one row per cluster,
%                     each the per-cluster median of T_est), ClusterSize
%                     (estimation count per cluster), Silhouette (mean
%                     silhouette score of the chosen split), and Dominant
%                     (index of the largest cluster, the one T.Range/
%                     T.index reflect when NumClusters > 1). Single-basin
%                     default (NumClusters = 1, Dominant = 1) when
%                     cluster is false or no split beat silThreshold.
%
%   Example:
%     [T,ci] = gsua_ia(T,Table.Estlsqc,false,false,true,true);
%     if ci.NumClusters > 1
%         disp('Candidate global minima (one row per cluster):')
%         disp(ci.Centers)
%     end
%
%   See also GSUA_DIA, GSUA_PE, GSUA_MEDIANCI, GSUA_COSTCUTOFF.

p = inputParser;
addRequired(p,'T');
addRequired(p,'T_est');
addOptional(p,'correction',false);
addOptional(p,'outlier',false);
addOptional(p,'show',false);
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

correction=p.Results.correction;
outlier=p.Results.outlier;
show=p.Results.show;
cluster=p.Results.cluster;
maxK=p.Results.maxK;
silThreshold=p.Results.silThreshold;
cost=p.Results.cost;

% --- Fit-quality filtering (BEFORE anything else, including clustering/outliers) ---
if ~isempty(cost)
    keepCost=gsua_costcutoff(cost,p.Results.CostMethod,p.Results.CostRtol, ...
        p.Results.CostAtol,p.Results.CostGapRatio,p.Results.MinKeep);
    nBadFit=sum(~keepCost);
    if nBadFit>0
        disp([num2str(nBadFit) ' run(s) dropped by the fit-quality filter (converged to a much worse cost than the best one)'])
    end
    T_est=T_est(:,keepCost);
end

D1 = floor(sqrt(size(T,1))); % Number of rows of subplot
D2 = D1+ceil((size(T,1)-D1^2)/D1);

NormalizedAll=zeros(size(T,1),size(T_est,2));
for i=1:size(T,1)
    NormalizedAll(i,:)=(T_est(i,:)-T.Range(i,1))/(T.Range(i,2)-T.Range(i,1));
end

% --- Clustering runs BEFORE outlier removal -- see header comment for why ---
clusterInfo=struct('NumClusters',1,'Idx',ones(1,size(T_est,2)),'Centers',median(T_est,2)', ...
    'ClusterSize',size(T_est,2),'Silhouette',[],'Dominant',1);
statT_est=T_est; % point set stats/plots below are computed from; narrows to the dominant cluster if found
X=NormalizedAll'; % one row per estimation run, one column per parameter (used below and by cluster plots)
if cluster
    Nruns=size(NormalizedAll,2);
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

% --- Outlier removal AFTER clustering, scoped to the dominant cluster ---
if outlier
    disp('Removing outliers...')
    [~,~,RD,chi_crt]=DetectMultVarOutliers(statT_est');
    id_in=RD<chi_crt(4);
    statT_est=statT_est(:,id_in);
    disp(num2str(sum(id_in))+" outliers were removed")
end

% --- Everything below is computed from statT_est (final filtered set) ---
Normalized=zeros(size(T,1),size(statT_est,2));
for i=1:size(T,1)
    Normalized(i,:)=(statT_est(i,:)-T.Range(i,1))/(T.Range(i,2)-T.Range(i,1));
end
nnominal=(statT_est(:,1)-T.Range(:,1))./(T.Range(:,2)-T.Range(:,1));

RHO = corr(statT_est');

med=median(statT_est,2);
[lb, up] = gsua_medianCI(statT_est, 0.05);

if (any(lb(lb<T.Range(:,1))))||(any(up(up>T.Range(:,2))))
    action=correction;
    if action
        lb(lb<T.Range(:,1))= T.Range(lb<T.Range(:,1),1);
        up(up>T.Range(:,2))= T.Range(up>T.Range(:,2),2);
    end
end

boxin=(up-lb)./(T.Range(:,2)-T.Range(:,1));
len=length(boxin);
corrin=sum(abs(RHO))/len;
extrin=sum(abs(RHO)>0.5)/len;
ind=(2*boxin'+corrin+extrin)/4;
[~,idx]=sort(ind,'descend');

% --- Plots ---
figure('Name','CDF Range')
for i=1:size(T,1)
subplot(D1,D2,i)
ecdf(statT_est(i,:),'Bounds','on')
title(T.Properties.RowNames{i})
end

figure('Name','Normalized Practical identifiability')
clf
boxplot(Normalized','Labels',T.Properties.RowNames,'LabelOrientation','inline')
h = findobj(gca, 'type', 'text');
set(h, 'Interpreter', 'tex','FontWeight','bold');

B = filloutliers(Normalized','center','median');
indicator=max(B,[],1)-min(B,[],1);
[~,index]=sort(indicator,'descend');
nn=B(:,index);
names=T.Properties.RowNames(index);
figure('Name','Sorted parameter range')
clf
boxplot(nn,'Labels',names,'LabelOrientation','inline')
h = findobj(gca, 'type', 'text');
set(h, 'Interpreter', 'tex','FontWeight','bold');
hold on
nnominal_sorted=nnominal(index);
plot(nnominal_sorted,'.','MarkerSize',20,'Color','black')

figure('Name','Correlation')
imagesc(RHO)
xticks(1:size(T,1))
xticklabels(T.Properties.RowNames)
yticks(1:size(T,1))
yticklabels(T.Properties.RowNames)
clim([-1 1])
colormap jet
colorbar

if show

    figure('Name','Identifiability Index')
    scatter3(boxin(idx(1)),corrin(idx(1)),extrin(idx(1)),60,ind(idx(1)),'filled','DisplayName',T.Properties.RowNames{idx(1)})
    hold on
    for i=idx(2:end)
        scatter3(boxin(i),corrin(i),extrin(i),60,ind(i),'filled','DisplayName',T.Properties.RowNames{i})
    end
    colormap jet
    legend('NumColumns',2','Orientation','vertical')
    xlabel('Interval Index')
    ylabel('Correlation Index')
    zlabel('Strong correlation Index')


    figure('Name','Identifiability graph')
    G=graph((abs(RHO)>0.5)-eye(len),T.Properties.RowNames);
    G.Nodes.Weights=ind';
    G.Nodes.NodeColors = ind';
    p2=plot(G);
    p2.NodeCData = G.Nodes.NodeColors;
    p2.NodeFontSize=12;
    p2.MarkerSize=7;
    colormap jet
    a = colorbar;
    a.Label.String = 'Identifiability Index';
    title('Strong correlations')

    if cluster && clusterInfo.NumClusters>1
        figure('Name','Multiple global minima - clustering')
        clf
        [~,score]=pca(X);
        subplot(1,2,1)
        gscatter(score(:,1),score(:,2),clusterInfo.Idx)
        xlabel('PC1'); ylabel('PC2')
        title(['Estimation runs grouped into ' num2str(clusterInfo.NumClusters) ' clusters'])
        subplot(1,2,2)
        silhouette(X,clusterInfo.Idx)
        title(['Mean silhouette = ' num2str(clusterInfo.Silhouette,3)])

        figure('Name','Candidate global minima')
        clf
        boxplot(NormalizedAll','Labels',T.Properties.RowNames,'LabelOrientation','inline')
        hold on
        markers={'o','s','^','d','v','p','h'};
        for k=1:clusterInfo.NumClusters
            centerNorm=(clusterInfo.Centers(k,:)'-T.Range(:,1))./(T.Range(:,2)-T.Range(:,1));
            plot(centerNorm,markers{mod(k-1,numel(markers))+1},'MarkerSize',10,'LineWidth',2,...
                'DisplayName',['Cluster ' num2str(k) ' (n=' num2str(clusterInfo.ClusterSize(k)) ')'])
        end
        legend('show')
        title('Per-cluster candidate global minima (normalized parameters)')
    end
end

T.Est=statT_est;
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
