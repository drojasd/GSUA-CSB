function [T,clusterInfo] = gsua_ia(T,T_est,varargin)
%GSUA_IA Practical identifiability analysis with diagnostic plots
%
%   T = gsua_ia(T,T_est)
%   T = gsua_ia(T,T_est,correction,outlier,show)
%   [T,clusterInfo] = gsua_ia(T,T_est,correction,outlier,show,cluster)
%   [T,clusterInfo] = gsua_ia(T,T_est,correction,outlier,show,cluster,maxK,silThreshold)
%   [T,clusterInfo] = gsua_ia(...,'cost',cost,'CostMethod','gap')
%   [T,clusterInfo] = gsua_ia(...,'IndexMethod','spread','truth',truth,'xdata',xdata,'ydata',ydata)
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
%   (NumClusters = 1) -- a low silhouette score means the "best" split
%   found is not actually well separated, so the sample is more likely
%   unimodal than confidently multimodal. Requires Statistics and Machine
%   Learning Toolbox.
%
%   Multistart vs bootstrap ensemble:
%   Feeding gsua_ia a MULTISTART ensemble (many fits of the SAME data)
%   measures whether the optimizer found the global optimum, not parameter
%   uncertainty -- it collapses to the best optimum. For honest
%   uncertainty intervals and correlations, feed a PARAMETRIC-BOOTSTRAP
%   ensemble instead (GSUA_BOOTSTRAP: refits of synthetic datasets drawn
%   from the best fit's own observation noise), whose spread reflects how
%   much the data noise moves each estimate.
%
%   Fit-quality filtering ('cost'):
%   A repeated-estimation run that converged to a much worse cost than
%   the best one contributes essentially arbitrary parameter values --
%   left in, it can make a genuinely well-identified parameter look
%   poorly identified purely because an optimizer run failed. Passing
%   'cost' (e.g. the res output of gsua_pe, one entry per column of
%   T_est) screens those runs out before anything else is computed. See
%   GSUA_COSTCUTOFF for the two available methods ('rtol'/'gap') and the
%   minKeep safety floor. The surviving columns are reported in
%   clusterInfo.Kept (see below).
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
%   the "Candidate global minima" plot and the per-cluster fields of
%   clusterInfo, which cover every cluster. Outlier removal, if
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
%     'IndexMethod'  <-- (optional, paired) 'spread' (default) or 'ci'. How the
%                        identifiability index measures interval width. 'spread'
%                        uses the 2.5-97.5 percentile dispersion of the estimates
%                        relative to the ORIGINAL range, so a parameter scattered
%                        across its whole prior is not reported as identified even
%                        when its median is well pinned (e.g. an initial-condition
%                        factor the data barely constrains). 'ci' is the legacy
%                        behaviour: median-CI width over the original range, which
%                        shrinks as the sample grows and can rate a flat parameter
%                        as identified. The correlation terms are unchanged.
%     'truth'        <-- (optional, paired) Np x 1 reference ("true") parameter
%                        values (e.g. for synthetic teaching data). When given and
%                        show is true, they are marked on the CDFs, the boxplots and
%                        the per-cluster scatter. Default: [] (none)
%     'xdata'        <-- (optional, paired) domain for the per-cluster curve plot
%                        (#7); required together with 'ydata' to draw it.
%     'ydata'        <-- (optional, paired) observed data for the per-cluster curve
%                        plot: inputs x numel(xdata). With show and >1 cluster, each
%                        cluster's center curve is overlaid on the data.
%
%   Outputs:
%     T           <-- same table with T.Range replaced by the new
%                     confidence interval, T.Nominal set to the median
%                     estimate, T.index set to the per-parameter
%                     identifiability index (0 = well identified, 1 =
%                     poorly identified; see IndexMethod), and T.Est set
%                     to the estimates actually used (after cost filtering,
%                     cluster restriction and outlier removal -- may be a
%                     subset of the input T_est). Reflects the dominant
%                     cluster when NumClusters > 1.
%     clusterInfo <-- struct describing the clustering and per-cluster
%                     statistics:
%                       NumClusters   number of candidate global minima
%                       Idx           1xN cluster label per ORIGINAL T_est
%                                     column; NaN for columns dropped by
%                                     the cost filter
%                       Kept          1xN logical over the ORIGINAL T_est
%                                     columns: which survived the cost
%                                     filter (all true if no 'cost' given)
%                       Centers       NumClusters x Np per-cluster medians
%                                     (candidate global minima)
%                       ClusterSize   estimation count per cluster
%                       Silhouette    mean silhouette of the chosen split
%                       SilhouetteByK mean silhouette for each tested K
%                       SilhouetteK   the K values SilhouetteByK refers to
%                       Dominant      index of the largest cluster (the one
%                                     T.Range/T.index reflect)
%                       Range         1xNumClusters cell; Range{c} is the
%                                     Np x 2 median CI of cluster c
%                       Index         1xNumClusters cell; Index{c} is the
%                                     Np x 1 identifiability index of cluster c
%                       Est           1xNumClusters cell; Est{c} is the
%                                     Np x n_c estimates of cluster c
%                       Median        1xNumClusters cell; Median{c} is the
%                                     Np x 1 median of cluster c
%                     Single-basin default (NumClusters = 1) when cluster
%                     is false or no split beat silThreshold; the per-cluster
%                     cells then hold one entry equal to the dominant stats.
%
%   Example:
%     [T,ci] = gsua_ia(T,Table.Estlsqc,false,false,true,true);
%     if ci.NumClusters > 1
%         disp('Candidate global minima (one row per cluster):')
%         disp(ci.Centers)
%         disp('Per-cluster identifiability index of cluster 2:')
%         disp(ci.Index{2})
%     end
%
%   See also GSUA_DIA, GSUA_PE, GSUA_BOOTSTRAP, GSUA_MEDIANCI, GSUA_COSTCUTOFF, GSUA_CLUSTERSUMMARY.

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
addParameter(p,'IndexMethod','spread');
addParameter(p,'truth',[]);
addParameter(p,'xdata',[]);
addParameter(p,'ydata',[]);
parse(p,T,T_est,varargin{:});

correction=p.Results.correction;
outlier=p.Results.outlier;
show=p.Results.show;
cluster=p.Results.cluster;
maxK=p.Results.maxK;
silThreshold=p.Results.silThreshold;
cost=p.Results.cost;
indexMethod=validatestring(p.Results.IndexMethod,{'spread','ci'});
truth=p.Results.truth;
curveX=p.Results.xdata;
curveY=p.Results.ydata;

% --- Fit-quality filtering (BEFORE anything else, including clustering/outliers) ---
N0=size(T_est,2);        % original number of runs (clusterInfo.Idx/.Kept index these)
origRange=T.Range;       % original bounds: the index denominator (T.Range is overwritten below)
keepCost=true(1,N0);
if ~isempty(cost)
    keepCost=gsua_costcutoff(cost,p.Results.CostMethod,p.Results.CostRtol, ...
        p.Results.CostAtol,p.Results.CostGapRatio,p.Results.MinKeep);
    keepCost=logical(keepCost(:)');
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
    NormalizedAll(i,:)=(T_est(i,:)-origRange(i,1))/(origRange(i,2)-origRange(i,1));
end

% --- Clustering runs BEFORE outlier removal -- see header comment for why ---
% cIdxF holds cluster labels over the COST-FILTERED columns (what the plots and per-cluster stats
% below use); clusterInfo.Idx (built further down) maps these back onto the original N columns.
cIdxF=ones(1,size(T_est,2));
centers=median(T_est,2)';
clusterSize=size(T_est,2);
bestSil=[];
silByK=[];
silK=[];
dominant=1;
X=NormalizedAll'; % one row per estimation run (used below and by cluster plots)
if cluster
    Nruns=size(NormalizedAll,2);
    kmax=min(maxK,Nruns-1);
    if kmax>=2
        silK=2:kmax;
        try
            eva=evalclusters(X,@(x,k) spectralcluster(x,k),'silhouette','KList',silK);
            silByK=eva.CriterionValues(:)';
            bestK=eva.OptimalK;
            bestSil=max(eva.CriterionValues,[],'omitnan');
        catch ME
            warning(ME.message)
            bestK=1;
            bestSil=0;
            silByK=[];
            silK=[];
        end
        if bestK>1 && bestSil>=silThreshold
            cIdxF=spectralcluster(X,bestK)';
            centers=zeros(bestK,size(T_est,1));
            clusterSize=zeros(1,bestK);
            for k=1:bestK
                centers(k,:)=median(T_est(:,cIdxF==k),2)';
                clusterSize(k)=sum(cIdxF==k);
            end
            [~,dominant]=max(clusterSize);
            disp(['Detected ' num2str(bestK) ' candidate global minima (mean silhouette ' num2str(bestSil,3) ')'])
            disp(['Cluster sizes (of ' num2str(Nruns) ' runs): ' num2str(clusterSize)])
        else
            disp(['No confident multi-cluster split found (best mean silhouette ' num2str(bestSil,3) ' < ' num2str(silThreshold) ')'])
            disp('Treating estimations as a single basin (1 candidate global minimum)')
        end
    else
        warning('Not enough estimation runs to test clustering (need at least 3 for maxK>=2)')
    end
end

% --- Dominant cluster point set, then outlier removal scoped to it ---
statT_est=T_est(:,cIdxF==dominant);
if outlier
    disp('Removing outliers...')
    [~,~,RD,chi_crt]=DetectMultVarOutliers(statT_est');
    id_in=RD<chi_crt(4);
    statT_est=statT_est(:,id_in);
    disp(num2str(sum(id_in))+" outliers were removed")
end

% --- Dominant-cluster range + identifiability index (shared helper GSUA_IASTATS) ---
Normalized=zeros(size(T,1),size(statT_est,2));
for i=1:size(T,1)
    Normalized(i,:)=(statT_est(i,:)-origRange(i,1))/(origRange(i,2)-origRange(i,1));
end
nnominal=(statT_est(:,1)-origRange(:,1))./(origRange(:,2)-origRange(:,1));

[lb,up,ind,med,RHO,boxin,corrin,extrin]=gsua_iastats(statT_est,origRange,indexMethod,'orderstat',correction);
len=size(statT_est,1);
[~,idx]=sort(ind,'descend');

% --- Per-cluster statistics (#1): same Range/Index/Est/Median the dominant cluster gets, for every
% cluster, so a multi-basin result exposes each basin's own interval and index. ---
Kc=max(cIdxF);
clRange=cell(1,Kc); clIndex=cell(1,Kc); clEst=cell(1,Kc); clMedian=cell(1,Kc);
for c=1:Kc
    est_c=T_est(:,cIdxF==c);
    clEst{c}=est_c;
    clMedian{c}=median(est_c,2);
    [lbc,upc,indc]=gsua_iastats(est_c,origRange,indexMethod,'orderstat',correction);
    clRange{c}=[lbc,upc];
    clIndex{c}=indc';
end

% --- clusterInfo (#1 per-cluster, #2 Kept/Idx over original columns, #6 silhouette-by-K data) ---
IdxFull=nan(1,N0);
IdxFull(keepCost)=cIdxF;
clusterInfo=struct('NumClusters',Kc,'Idx',IdxFull,'Kept',keepCost, ...
    'Centers',centers,'ClusterSize',clusterSize,'Silhouette',bestSil,'Dominant',dominant);
clusterInfo.SilhouetteByK=silByK;
clusterInfo.SilhouetteK=silK;
clusterInfo.Range=clRange;
clusterInfo.Index=clIndex;
clusterInfo.Est=clEst;
clusterInfo.Median=clMedian;
clusterInfo.OrigRange=origRange;   % original bounds (GSUA_CLUSTERSUMMARY's percentile-width denominator)

% --- Finalize T (before plotting, so a plotting error cannot cost the caller the results) ---
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

% --- Plots (guarded: a plotting failure warns but never discards the computed T/clusterInfo) ---
try
    % Parameter names may be TeX or LaTeX; gsua_labels converts them and picks the interpreter.
    [iaNames, iaInterp] = gsua_labels(T.Properties.RowNames);

    if ~isempty(truth)
        truth=truth(:);
        truthNorm=(truth-origRange(:,1))./(origRange(:,2)-origRange(:,1));
    else
        truthNorm=[];
    end

    figure('Name','CDF Range')
    for i=1:size(T,1)
        subplot(D1,D2,i)
        ecdf(statT_est(i,:),'Bounds','on')
        title(iaNames{i},'Interpreter',iaInterp)
        if ~isempty(truth)
            xline(truth(i),'r--','LineWidth',1.2);
        end
    end

    figure('Name','Normalized Practical identifiability')
    clf
    boxplot(Normalized','Labels',iaNames,'LabelOrientation','inline')
    h = findobj(gca, 'type', 'text');
    set(h, 'Interpreter', iaInterp,'FontWeight','bold');
    if ~isempty(truthNorm)
        hold on
        plot(1:numel(truthNorm),truthNorm,'r.','MarkerSize',20)
        hold off
    end

    B = filloutliers(Normalized','center','median');
    indicator=max(B,[],1)-min(B,[],1);
    [~,index]=sort(indicator,'descend');
    nn=B(:,index);
    names=iaNames(index);
    figure('Name','Sorted parameter range')
    clf
    boxplot(nn,'Labels',names,'LabelOrientation','inline')
    h = findobj(gca, 'type', 'text');
    set(h, 'Interpreter', iaInterp,'FontWeight','bold');
    hold on
    nnominal_sorted=nnominal(index);
    plot(nnominal_sorted,'.','MarkerSize',20,'Color','black')
    if ~isempty(truthNorm)
        plot(truthNorm(index),'rx','MarkerSize',10,'LineWidth',1.5)
    end
    hold off

    figure('Name','Correlation')
    imagesc(RHO)
    xticks(1:size(T,1))
    xticklabels(iaNames)
    yticks(1:size(T,1))
    yticklabels(iaNames)
    set(gca,'TickLabelInterpreter',iaInterp)
    clim([-1 1])
    colormap jet
    colorbar

    if show

        figure('Name','Identifiability Index')
        scatter3(boxin(idx(1)),corrin(idx(1)),extrin(idx(1)),60,ind(idx(1)),'filled','DisplayName',iaNames{idx(1)})
        hold on
        for i=idx(2:end)
            scatter3(boxin(i),corrin(i),extrin(i),60,ind(i),'filled','DisplayName',iaNames{i})
        end
        colormap jet
        legend('NumColumns',2','Orientation','vertical','Interpreter',iaInterp)
        xlabel('Interval Index')
        ylabel('Correlation Index')
        zlabel('Strong correlation Index')


        figure('Name','Identifiability graph')
        G=graph((abs(RHO)>0.5)-eye(len),iaNames);
        G.Nodes.Weights=ind';
        G.Nodes.NodeColors = ind';
        p2=plot(G,'Interpreter',iaInterp);
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
            gscatter(score(:,1),score(:,2),cIdxF)
            xlabel('PC1'); ylabel('PC2')
            title(['Estimation runs grouped into ' num2str(clusterInfo.NumClusters) ' clusters'])
            subplot(1,2,2)
            silhouette(X,cIdxF)
            title(['Mean silhouette = ' num2str(clusterInfo.Silhouette,3)])

            % #6: mean silhouette vs candidate K, with the acceptance threshold.
            if ~isempty(silByK)
                figure('Name','Silhouette vs K')
                clf
                plot(silK,silByK,'-o','LineWidth',1.5,'MarkerFaceColor','auto')
                hold on
                yline(silThreshold,'r--','silThreshold','LabelHorizontalAlignment','left');
                xlabel('Number of clusters K'); ylabel('Mean silhouette')
                title('Cluster-count selection'); grid on
                xticks(silK)
            end

            % #7a: pairwise scatter matrix of the cost-filtered estimates, colored by cluster.
            figure('Name','Per-cluster pairwise scatter')
            clf
            try
                gplotmatrix(T_est',[],cIdxF(:),[],[],[],'on','hist',iaNames)
            catch
                gplotmatrix(T_est',[],cIdxF(:))
            end

            % #7b: each cluster's center curve vs the data (requires xdata/ydata).
            if ~isempty(curveX) && ~isempty(curveY)
                figure('Name','Per-cluster curves vs data')
                clf
                nOut=size(curveY,1);
                cmap=lines(clusterInfo.NumClusters);
                for o=1:nOut
                    subplot(nOut,1,o)
                    plot(curveX,curveY(o,:),'k.','MarkerSize',8,'DisplayName','data'); hold on
                    for k=1:clusterInfo.NumClusters
                        yk=gsua_deval(clusterInfo.Centers(k,:),T,curveX);
                        plot(curveX,yk(o,:),'LineWidth',1.5,'Color',cmap(k,:), ...
                            'DisplayName',['Cluster ' num2str(k)])
                    end
                    hold off; legend('show'); xlabel('x'); ylabel(['output ' num2str(o)])
                end
            end

            figure('Name','Candidate global minima')
            clf
            boxplot(NormalizedAll','Labels',iaNames,'LabelOrientation','inline')
            set(findobj(gca,'type','text'),'Interpreter',iaInterp)
            hold on
            markers={'o','s','^','d','v','p','h'};
            for k=1:clusterInfo.NumClusters
                centerNorm=(clusterInfo.Centers(k,:)'-origRange(:,1))./(origRange(:,2)-origRange(:,1));
                plot(centerNorm,markers{mod(k-1,numel(markers))+1},'MarkerSize',10,'LineWidth',2,...
                    'DisplayName',['Cluster ' num2str(k) ' (n=' num2str(clusterInfo.ClusterSize(k)) ')'])
            end
            if ~isempty(truthNorm)
                plot(truthNorm,'rx','MarkerSize',12,'LineWidth',2,'DisplayName','truth')
            end
            legend('show')
            title('Per-cluster candidate global minima (normalized parameters)')
        end
    end
catch plotME
    warning('gsua_ia:PlottingFailed','Results computed; a diagnostic plot failed: %s', plotME.message);
end

end
