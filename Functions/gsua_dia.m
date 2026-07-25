function [T,clusterInfo] = gsua_dia(T,T_est,outlier,cluster,maxK,silThreshold)
%GSUA_DIA Practical identifiability diagnostics (headless, no plots)
%
%   T = gsua_dia(T,T_est)
%   T = gsua_dia(T,T_est,outlier)
%   [T,clusterInfo] = gsua_dia(T,T_est,outlier,cluster)
%   [T,clusterInfo] = gsua_dia(T,T_est,outlier,cluster,maxK,silThreshold)
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
%   Inputs:
%     T            <-- summary table from gsua_dataprep
%     T_est        <-- Np x N matrix, one column per repeated parameter
%                      estimation (e.g. Table.Estlsqc from gsua_pe)
%     outlier      <-- (optional, logical) remove multivariate outliers
%                      from T_est before computing the index. Default: false
%     cluster      <-- (optional, logical) run spectral clustering to
%                      detect multiple global minima. Default: false
%     maxK         <-- (optional) largest cluster count to test. Default: 5
%     silThreshold <-- (optional) minimum mean silhouette score required
%                      to accept a multi-cluster split. Default: 0.5
%
%   Outputs:
%     T           <-- same table with T.Range replaced by the new
%                     confidence interval, T.Nominal set to the median
%                     estimate, and T.index set to the per-parameter
%                     identifiability index (0 = well identified, 1 =
%                     poorly identified)
%     clusterInfo <-- struct with fields NumClusters, Idx, Centers
%                     (candidate global minima, one row per cluster),
%                     ClusterSize, and Silhouette. Empty split
%                     (NumClusters = 1) when cluster is false.
%
%   See also GSUA_IA, GSUA_PE, GSUA_MEDIANCI.

if nargin<6
    silThreshold=0.5;
end
if nargin<5
    maxK=5;
end
if nargin<4
    cluster=false;
end
if nargin <3
    outlier=false;
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


Normalized=zeros(size(T,1),size(T.Est,2));
for i=1:size(T,1)
    Normalized(i,:)=(T.Est(i,:)-T.Range(i,1))/(T.Range(i,2)-T.Range(i,1));
end

RHO = corr(T.Est');

x=T.Est;
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

clusterInfo=struct('NumClusters',1,'Idx',ones(1,size(T_est,2)),'Centers',med','ClusterSize',size(T_est,2),'Silhouette',[]);
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