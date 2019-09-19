function [ParU,New_range,J_test,Y_test,sup]= gsua_ucis(Par,N,varargin)
% Function for uncertainty-based confidence intervals
%
% T= gsua_ucis(T,N)
% Parameters:
% T <-- summary table from gsua_oatr function
% N <-- Number of samples per cycle
% Outputs:
% T <-- summary table with new confidence intervals in T.Range
% Additional features:
% It is possible to request the following additional positional outputs
% New_range <-- array with range modification record.
% J_test <-- MSE output of last iteration
% Y_test <-- output of last iteration
% sup    <-- record of good scalars
% Also, it is possible to apply the next paired features
% 'ynom',ynom <-- to especify a new output objective
% 'reps',k <--to perform k cycles. Default: 100.
% 'recort',r <-- to set another recort criteria (r). Default: 0.5.
% 'select',s <-- to set another select criteria (s). Default: 0.5.
% 'lim',l <-- to set another distance criteria (l). Default: 0.3.
% 'stop',p <-- to set the confidence of the interval in (p). Default: 0.95.
% 'parallel',false <-- to avoid the speed up of parallel computing toolbox.
% [T,New_range,J_test,Y_test,sup]= gsua_ucis(T,N,'lim',0.5,'select',0.4)

p=inputParser;

defaultReps=100;
defaultRecort=0.5; 
defaultSelect=0.5;
defaultLim=0.3;
defaultStop=0.95;
defaultynom=[];
defaultParallel=true;

addRequired(p,'Par',@istable);
addRequired(p,'N',@isnumeric);
addParameter(p,'ynom',defaultynom,@isnumeric);
addParameter(p,'reps',defaultReps,@isnumeric);
addParameter(p,'recort',defaultRecort,@isnumeric);
addParameter(p,'select',defaultSelect,@isnumeric);
addParameter(p,'lim',defaultLim,@isnumeric);
addParameter(p,'stop',defaultStop,@isnumeric);
addParameter(p,'parallel',defaultParallel,@islogical);

parse(p,Par,N,varargin{:})
Par=p.Results.Par;
try
    TP=Par.Properties.CustomProperties;
catch
    TP=load('ATable.mat');
    TP=TP.Table2;
end
ParU=Par;
N=p.Results.N;
y_Nom=p.Results.ynom;
reps=p.Results.reps;
recort=p.Results.recort;
select=p.Results.select;
lim=p.Results.lim;
stop=p.Results.stop;
parallel=p.Results.parallel;



kind=TP.Kind;

if isempty(y_Nom)
    if strcmp(kind,'mat')
        domain=TP.Domain;
        xdata=linspace(domain(1),domain(2),(domain(2)-domain(1)));
    else
        time=TP.tout;
        xdata=min(time):max(time);
    end
    y_Nom=gsua_deval(Par.Nominal',Par,xdata);
else
    xdata=linspace(0,size(y_Nom)-1,size(y_Nom));
end

if nargin<9
    stop=0.95;
end

Np=size(Par,1);%Number of parameters
try
    fixed=TP.Fixed;
    if isempty(fixed)
        fixed=false(1,Np);
    end
catch
    fixed=false(1,Np);
end
bins_div=divisors(N);
bins_div=bins_div(bins_div>=10);
if Np>10
    bins_div=bins_div(bins_div<=ceil(N/Np));
else
    bins_div=bins_div(bins_div<=ceil(N/10));
end

or_nBins=bins_div(1);
or_select=select;

%Par to be modified
New_range=cell(reps,Np);%Ranges per iteration
D1 = floor(sqrt(Np));%Number of subplots
D2 = D1+ceil((Np-D1^2)/D1);%Number of subplots
sup=zeros(reps,1);
noms=Par.Nominal';
Range=Par.Range;

%Ranges expansion


limit=sum((y_Nom-(y_Nom*(1+lim))).^2,2);
Best_esc=[0,noms];%Parameters of the best curve

disp('------------------------------------')
disp("Ranges calculation is being launched")    
disp('------------------------------------') 

for i=1:reps
    intern_counter=1;
    intern_counter2=1;
    [M_test] = gsua_dmatrix(ParU,N);
    Y_test=gsua_pardeval(M_test,Par,xdata,parallel,[(i-1)*N,reps*N]);
    J_test=sum((Y_test-y_Nom).^2,2);
    BFM_esc=[J_test,M_test,Y_test];
    BFM_esc=sortrows(BFM_esc); 
    nBins=or_nBins;
    select=or_select;
    
    figure(2)
    clf
    try
        sup(i)=find(BFM_esc(:,1)<limit, 1, 'last' );  
        gsua_plot('UncertaintyAnalysis',Par,BFM_esc(:,Np+2:end),xdata,y_Nom,sup(i));
    catch
        sup(i)=0;
        gsua_plot('UncertaintyAnalysis',Par,BFM_esc(:,Np+2:end),xdata,y_Nom);
    end
    drawnow
    
    disp('------------------------------------')
    disp([num2str(sup(i)) ' Scalars has been detected in the especified range'])
    disp('------------------------------------')

      
    if sup(i)>N*stop %criterio de parada
        disp('At least the 95% of curves within the specified range')
        Par.Range=ParU.Range;
        return
    end
            
    if sup(i)>ceil(0.75*N)
        method=2;
    else
        method=1;
    end
    
    TRUE=1;
    while TRUE==1
        Values=ones(1,nBins+1);%Values of beans
        
        if method==2 
            BFM_esc_best=BFM_esc(1:sup(i),:);
        else      
            BFM_esc_best=BFM_esc(1:(ceil(N*select)),:);    
        end

        figure(3)
        clf
        for w=1:Np
            if ~(fixed(w))
            subplot(D1,D2,w)
            bins=histogram(BFM_esc_best(:,w+1),nBins);
            if i>1
                bins.BinEdges = New_range{i-1,w}(1):(New_range{i-1,w}(2)-New_range{i-1,w}(1))/nBins:New_range{i-1,w}(2);
            else
                bins.BinEdges = Par.Range(w,1):(Par.Range(w,2)-Par.Range(w,1))/nBins:Par.Range(w,2);
            end
            hold on
            plot(Par.Nominal(w),0,'o','MarkerSize',10,'MarkerEdgeColor','r')
            xlabel(Par.Properties.RowNames(w))


            Values(1:end-1)=bins.Values;
            Values(end)=Values(end-1);

            New_hist=bins.BinEdges(Values>(max(Values)*recort));
            New_range{i,w}=[min(New_hist);max(New_hist)];


            if New_range{i,w}(1)>=Best_esc(w+1)
                New_range{i,w}(1)=Best_esc(w+1)-(bins.BinEdges(end)-bins.BinEdges(1))*0.1;%Edge scape with 5%
                New_range{i,w}(2)=New_range{i,w}(2)-(bins.BinEdges(end)-bins.BinEdges(1))*0.1;

            elseif New_range{i,w}(2)<=Best_esc(w+1)
                New_range{i,w}(2)=Best_esc(w+1)+(bins.BinEdges(end)-bins.BinEdges(1))*0.1;%Edge scape with 5%
                New_range{i,w}(1)=New_range{i,w}(1)+(bins.BinEdges(end)-bins.BinEdges(1))*0.1;
            end
            end
        end

        if i>1 
            cuts=sum(sum([New_range{i,:}]==[New_range{i-1,:}]))+sum(fixed)*2;
        else
            cuts=sum(sum([New_range{i,:}]==Par.Range(~fixed,:)'))+sum(fixed)*2;
        end  
        if cuts<Np*2
            disp([num2str(Np*2-cuts) '--> Reductions done in this iteration'])
        end

        if cuts==Np*2

            try
                nBins=bins_div(intern_counter+1);
            catch
                %nBins=bins_div(end);
                method=1;
                if sup(i)>select
                    select=sup(i)/N*0.9^intern_counter2;%% solucionar
                else
                    select=select/N*0.9^intern_counter2;
                end
                disp('------------------------------')
                disp(["Changing Select to " num2str(select)])
                disp('------------------------------')
                intern_counter2=intern_counter2+1;
                intern_counter=0;
                nBins=bins_div(1);
            end
            disp('------------------------------')
            disp(["Changing nBins to " num2str(nBins)])
            disp('------------------------------')
        else
            TRUE=0;
        end
        intern_counter=intern_counter+1;
    end
    
    for w=1:Np
        if ~fixed(w)
%         disp(['Initial range of (' num2str(w) '): ' Par.Properties.RowNames(w) '--> ' num2str(Range(w,:))])
%         disp(['Actual range --> ' num2str(New_range{i,w}')])
        ParU.Range(w,:)=[New_range{i,w}]';  
        end
    end
    
    disp('------------------------------')
    disp(["Running method " num2str(method)])
    disp('------------------------------')
    Normalized=zeros(2,Np);
    for h=1:Np
        if ~fixed(h)
        Normalized(:,h)=[(ParU.Range(h,1)-Range(h,1))/(Range(h,2)-Range(h,1));...
            (ParU.Range(h,2)-Range(h,1))/(Range(h,2)-Range(h,1))];
        end
    end
    figure(1)
    clf
    boxplot(Normalized(:,~fixed),'Labels',Par.Properties.RowNames(~fixed)')
    title(['Iteration # ' num2str(i) ', Scalars: ' num2str(sup(i))])
    drawnow
    
    figure(4)
    clf
    plot(1:i,sup(1:i)/N)
    title('# Scalars vs iteration')
    xlabel('Iteration')
    ylabel('Scalars (%)')
    save('Par_Final.mat','ParU');
end
end