function T = gsua_ia(T,T_est)
T.Estlsqc=T_est;
T.Nominal=T.Estlsqc(:,1);
D1 = floor(sqrt(size(T,1))); % Number of rows of subplot
D2 = D1+ceil((size(T,1)-D1^2)/D1);
figure('Name','CDF Range')
for i=1:size(T,1)
subplot(D1,D2,i)
ecdf(T.Estlsqc(i,:),'Bounds','on')
title(T.Properties.RowNames{i})
end

Normalized=zeros(size(T,1),size(T.Estlsqc,2));
for i=1:size(T,1)
    Normalized(i,:)=(T.Estlsqc(i,:)-T.Range(i,1))/(T.Range(i,2)-T.Range(i,1));
end
nnominal=(T.Nominal-T.Range(:,1))./(T.Range(:,2)-T.Range(:,1));
figure('Name','Normalized Practical identifiability')
clf
boxplot(Normalized','Labels',T.Properties.RowNames,'LabelOrientation','inline')


B = filloutliers(Normalized','center','median');
indicator=[max(B,[],1)-min(B,[],1)];
[~,index]=sort(indicator,'descend');
nn=B(:,index);
names=T.Properties.RowNames(index);
figure('Name','Sorted parameter range')
boxplot(nn,'Labels',names,'LabelOrientation','inline')
hold on
nnominal=nnominal(index);
plot(nnominal,'.','MarkerSize',20,'Color','black')


% figure('Name','NPI without outliers')
% clf
% boxplot(B,'Labels',T.Properties.RowNames,'LabelOrientation','inline')
% 
 nRange=zeros(size(T,1),2);
for i=1:size(nRange,1)
    nRange(i,:)=quantile(T.Estlsqc(i,:), [0.05,0.95]);
end
T.Range=nRange;

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