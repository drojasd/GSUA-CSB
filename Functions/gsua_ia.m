function [T,clusterInfo] = gsua_ia(T,T_est,correction,outlier,show,cluster,maxK,silThreshold)
%GSUA_IA Practical identifiability analysis with diagnostic plots
%
%   T = gsua_ia(T,T_est)
%   T = gsua_ia(T,T_est,correction,outlier,show)
%   [T,clusterInfo] = gsua_ia(T,T_est,correction,outlier,show,cluster)
%   [T,clusterInfo] = gsua_ia(T,T_est,correction,outlier,show,cluster,maxK,silThreshold)
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
%   Inputs:
%     T            <-- summary table from gsua_dataprep
%     T_est        <-- Np x N matrix, one column per repeated parameter
%                      estimation (e.g. Table.Estlsqc from gsua_pe)
%     correction   <-- (optional, logical) clip the new range to the
%                      original T.Range when it would otherwise expand
%                      beyond it. Default: false
%     outlier      <-- (optional, logical) remove multivariate outliers
%                      from T_est before computing statistics. Default: false
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
%
%   Outputs:
%     T           <-- same table with T.Range replaced by the new
%                     confidence interval, T.Nominal set to the median
%                     estimate, and T.index set to the per-parameter
%                     identifiability index (0 = well identified, 1 =
%                     poorly identified)
%     clusterInfo <-- struct with fields NumClusters, Idx (1xN cluster
%                     label per estimation run), Centers (NumClusters x
%                     Np candidate global minima, one row per cluster,
%                     each the per-cluster median of T_est), ClusterSize
%                     (estimation count per cluster), and Silhouette
%                     (mean silhouette score of the chosen split). Empty
%                     struct with NumClusters = 1 when cluster is false.
%
%   Example:
%     [T,ci] = gsua_ia(T,Table.Estlsqc,false,false,true,true);
%     if ci.NumClusters > 1
%         disp('Candidate global minima (one row per cluster):')
%         disp(ci.Centers)
%     end
%
%   See also GSUA_DIA, GSUA_PE, GSUA_MEDIANCI.

if nargin<8
    silThreshold=0.5;
end
if nargin<7
    maxK=5;
end
if nargin<6
    cluster=false;
end
if nargin<5
    show=false;
end
if nargin<4
    outlier=false;
end
if nargin <3
    correction=false;
end
if outlier
    disp('Removing outliers...')
    [~,~,RD,chi_crt]=DetectMultVarOutliers(T_est');
    id_in=RD<chi_crt(4);
    T_est=T_est(:,id_in);
    disp(num2str(sum(id_in))+" outliers were removed")
end
T.Est=T_est;
T.Nominal=T.Est(:,1);
D1 = floor(sqrt(size(T,1))); % Number of rows of subplot
D2 = D1+ceil((size(T,1)-D1^2)/D1);
figure('Name','CDF Range')
for i=1:size(T,1)
subplot(D1,D2,i)
ecdf(T.Est(i,:),'Bounds','on')
title(T.Properties.RowNames{i})
end

Normalized=zeros(size(T,1),size(T.Est,2));
for i=1:size(T,1)
    Normalized(i,:)=(T.Est(i,:)-T.Range(i,1))/(T.Range(i,2)-T.Range(i,1));
end
nnominal=(T.Nominal-T.Range(:,1))./(T.Range(:,2)-T.Range(:,1));
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
nnominal=nnominal(index);
plot(nnominal,'.','MarkerSize',20,'Color','black')

figure('Name','Correlation')
RHO = corr(T.Est');
imagesc(RHO)
xticks(1:size(T,1))
xticklabels(T.Properties.RowNames)
%xtickangle(60)
yticks(1:size(T,1))
yticklabels(T.Properties.RowNames)
clim([-1 1])
colormap jet
colorbar
% figure('Name','NPI without outliers')
% clf
% boxplot(B,'Labels',T.Properties.RowNames,'LabelOrientation','inline')
% 

x=T.Est;
med=median(x,2);
% N=size(x,2);
% desv=std(x,[],2);

% lb = med-1.96*sqrt(pi/2)*desv/sqrt(N);
% up = med+1.96*sqrt(pi/2)*desv/sqrt(N);

[lb, up] = gsua_medianCI(x, 0.05);

if (any(lb(lb<T.Range(:,1))))||(any(up(up>T.Range(:,2))))
    %action=input('Replace out ranges by estimation ranges? (true,false): ');
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

clusterInfo=struct('NumClusters',1,'Idx',ones(1,size(T_est,2)),'Centers',med','ClusterSize',size(T_est,2),'Silhouette',[]);
if cluster
    Nruns=size(Normalized,2);
    X=Normalized'; % one row per estimation run, one column per parameter
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
            clusterInfo=struct('NumClusters',bestK,'Idx',cIdx,'Centers',centers,'ClusterSize',clusterSize,'Silhouette',bestSil);
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

if show

    figure('Name','Identifiability Index')
    scatter3(boxin(idx(1)),corrin(idx(1)),extrin(idx(1)),60,ind(idx(1)),'filled','DisplayName',T.Properties.RowNames{idx(1)})
    hold on
    for i=idx(2:end)
        scatter3(boxin(i),corrin(i),extrin(i),60,ind(i),'filled','DisplayName',T.Properties.RowNames{i})
    end
    colormap jet
    % b=colorbar('Location','eastoutside');
    % b.Label.String = 'Identifiability Index';
    legend('NumColumns',2','Orientation','vertical')
    xlabel('Interval Index')
    ylabel('Correlation Index')
    zlabel('Strong correlation Index')


    figure('Name','Identifiability graph')
    G=graph((abs(RHO)>0.5)-eye(len),T.Properties.RowNames);
    G.Nodes.Weights=ind';
    G.Nodes.NodeColors = ind';
    p=plot(G);
    p.NodeCData = G.Nodes.NodeColors;
    p.NodeFontSize=12;
    p.MarkerSize=7;
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
        boxplot(Normalized','Labels',T.Properties.RowNames,'LabelOrientation','inline')
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


% indicator=[max(B,[],1)-min(B,[],1)];
% [~,index]=sort(indicator);
% nn=B(:,index);
% names=T.Properties.RowNames(index);
% val=median(nn,1);
% nn=[nn-val];
% nnominal=nnominal(index)'-val;
% figure('Name','Sorted parameter range')
% boxplot(nn,'Labels',names,'LabelOrientation','inline')
% hold on
% plot(nnominal,'Marker','*','MarkerSize',2)
end