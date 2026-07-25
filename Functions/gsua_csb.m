function [ParU,New_range,J_test,Y_test,sup,condition]= gsua_csb(Par,N,varargin)
% Function for uncertainty-based confidence intervals
%
% T= gsua_csb(T,N)
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
% [T,New_range,J_test,Y_test,sup]= gsua_csb(T,N,'lim',0.5,'select',0.4)
%
% See also GSUA_OATR, GSUA_OATR2, GSUA_DMATRIX, GSUA_PARDEVAL.

p=inputParser;

defaultReps=100;
defaultRecort=0.75; 
defaultSelect=0.5;
defaultLim=0.3;
defaultStop=0.95;
defaultynom=[];
defaultParallel=false;
defaultCorrect=true;
defaultBreaking=1;
defaultStretch=3;

addRequired(p,'Par',@istable);
addRequired(p,'N',@isnumeric);
addParameter(p,'ynom',defaultynom,@isnumeric);
addParameter(p,'reps',defaultReps,@isnumeric);
addParameter(p,'select',defaultSelect,@isnumeric);
addParameter(p,'lim',defaultLim,@isnumeric);
addParameter(p,'stop',defaultStop,@isnumeric);
addParameter(p,'parallel',defaultParallel,@islogical);
addParameter(p,'correct',defaultCorrect,@islogical);
addParameter(p,'protect',defaultParallel,@islogical);
addParameter(p,'breaking',defaultBreaking,@isnumeric);
addParameter(p,'show',defaultParallel,@islogical);
addParameter(p,'switch',defaultRecort,@isnumeric);
addParameter(p,'stretch',defaultStretch,@isnumeric);


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
select=p.Results.select;
lim=p.Results.lim;
stop=p.Results.stop;
parallel=p.Results.parallel;
correct=p.Results.correct;
protect=p.Results.protect;
breaking=p.Results.breaking;
show=p.Results.show;
recort=p.Results.switch;
stretch=round(p.Results.stretch);
condition=0;

kind=TP.Kind;

if isempty(y_Nom)
    if ismember(kind,[2 3 4 5])
        domain=TP.Domain;
        xdata=linspace(domain(1),domain(2),(domain(2)-domain(1)));
    else
        time=TP.Tout;
        xdata=min(time):max(time);
    end
    y_Nom=gsua_deval(Par.Nominal',Par,xdata);
else
    xdata=linspace(0,length(y_Nom)-1,length(y_Nom));
end

Np=size(Par,1);%Number of parameters
if N<Np*10
    N=Np*10;
    fprintf('Changing N to %u \n',N) 
end
try
    fixed=TP.Fixed;
    if isempty(fixed)
        fixed=false(1,Np);
    end
catch
    fixed=false(1,Np);
end

or_select=select;

%Par to be modified
New_range=zeros(Np,2,reps+1);%Ranges per iteration
D1 = floor(sqrt(Np));%Number of subplots
D2 = D1+ceil((Np-D1^2)/D1);%Number of subplots
sup=zeros(reps,1);
noms=Par.Nominal';
Range=Par.Range;
New_range(:,:,1)=ParU.Range;

%Ranges expansion

Best_esc=[0,noms];%Parameters of the best curve

if breaking==1
    disp('------------------------------------')
    disp("Ranges calculation is being launched")    
    disp('------------------------------------') 
end

for i=1:reps
    intern_counter=1;
    [M_test] = gsua_dmatrix(ParU,N);
    if breaking==1
        Y_test=gsua_pardeval(M_test,Par,xdata,parallel,[(i-1)*N,reps*N]);
    else
         Y_test=gsua_pardeval(M_test,Par,xdata,parallel);
    end
    J_test = gsua_costfMulti(Y_test,y_Nom,lim, parallel);
    BFM_esc=[J_test',M_test];
    [BFM_esc,order]=sortrows(BFM_esc); 
    select=or_select;

    try
        sup(i)=find(BFM_esc(:,1) < 1, 1, 'last' );  
    catch
        sup(i)=0;
    end
    drawnow
    
    if breaking==1
        disp('------------------------------------')
        disp([num2str(sup(i)) ' Scalars have been detected in the especified range'])
        disp('------------------------------------')
    end

      
    if sup(i)>N*stop %criterio de parada
        disp('At least the 95% of curves within the specified range')
        Par.Range=ParU.Range;
        return
    end
            
    if sup(i)>ceil(recort*N)
        method=2;
    else
        method=1;
    end
    
    TRUE=1;
    while TRUE==1

        if method==2 
            BFM_esc_best=BFM_esc(1:sup(i),:);
        else      
            BFM_esc_best=BFM_esc(1:(ceil(N*select)),:);    
        end
        for w=1:Np
            if ~(fixed(w))

            data = BFM_esc_best(:,w+1);
            maxreturn=New_range(w,2,i);
            minreturn=New_range(w,1,i);
            inter = [0,1];

%                     
            data = normalize(data,'range');
            dist = makedist('uniform');


            flag = kstest(data,'CDF',dist,'Alpha',0.1);
            %sskew=[];
            sp = [];
            mcounter =1;
            kurto = kurtosis(data);

            while flag~=0 && kurto>1.79
                skew = skewness(data);
                linter = inter(2)-inter(1);
                if abs(skew)<0.1
                    interL = inter(1)+linter*0.05;
                    interU = inter(2)-linter*0.05;
                    inter = [interL,interU];
                elseif skew>0
                    interU = inter(2)-linter*0.05;
                    interL = inter(1);
                    inter = [interL,interU];
                else
                    interL = inter(1)+linter*0.05;
                    interU = inter(2);
                    inter = [interL,interU];
                end
                data = data(data>inter(1));
                data = data(data<inter(2));

                [flag,p] = kstest(normalize(data,'range'),'CDF',dist,'Alpha',0.1);
                kurto = kurtosis(normalize(data,'range'));
                sp=[sp,p];

                if mcounter>2
                    if abs((sp(end-1)-sp(end))/sp(end-1))<2
                        break
                    end
                    if sp(end-1)>sp(end)
                        break
                    end
                end
                mcounter = mcounter+1;
            end

            New_range(w,:,i+1)=inter*(maxreturn-minreturn)+minreturn;

            if protect
                if New_range(w,1,i+1)>=Best_esc(w+1)
                    New_range(w,1,i+1)=Best_esc(w+1)-(New_range(w,2,i+1)-New_range(w,1,i+1))*0.1;%Edge scape with 5%
                    New_range(w,2,i+1)=New_range(w,2,i+1)-(New_range(w,2,i+1)-New_range(w,1,i+1))*0.1;
                    if correct
                        if New_range(w,1,i+1) < Range(w,1)
                            New_range(w,1,i+1) = Range(w,1);
                        end
                    end

                elseif New_range(w,2,i+1)<=Best_esc(w+1)
                    New_range(w,2,i+1)=Best_esc(w+1)+(New_range(w,2,i+1)-New_range(w,1,i+1))*0.1;%Edge scape with 5%
                    New_range(w,1,i+1)=New_range(w,1,i+1)+(New_range(w,2,i+1)-New_range(w,1,i+1))*0.1;
                    if correct
                        if New_range(w,2,i+1) > Range(w,2)
                            New_range(w,2,i+1) = Range(w,2);
                        end
                    end
                end
            end
            end
        end
            

        cuts=sum(sum(New_range(:,:,i)==New_range(:,:,i+1)))+sum(fixed)*2;

        if cuts<Np*2
            if breaking==1
                disp([num2str(Np*2-cuts) '--> Reductions done in this iteration'])
            end
        end

        if cuts==Np*2

            method=1;
            if sup(i)>select
                select=sup(i)/N*0.9^intern_counter;%% solucionar
            else
                select=select*0.9^intern_counter;
            end

            if select<0.1 || select*N<20
                if breaking<=stretch
                    N2 = round(N*2);
                    [AuxPar,~,~,~,~,cond]=gsua_csb(ParU,N2,'ynom',y_Nom,'reps',1,'select',or_select,...
                        'lim',lim,'stop',stop,'parallel',parallel,'protect',protect,'correct',correct,...
                        'switch',recort,'breaking',breaking+1,'stretch',stretch);
                    if cond==100
                        condition=cond;
                        if breaking==1
                            condition=-1*i;
                            disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
                            disp("Failure on iteration "+num2str(i))
                            disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
                        end
                        return
                    else
                        condition = condition+cond+1;
                        New_range(:,:,i+1) = AuxPar.Range;
                        TRUE = 0;
                    end
                else
                    condition=100;
                    return
                end
            else
                disp('------------------------------')
                disp(["Changing Select to " num2str(select)])
                disp('------------------------------')
            end
         
        else
            TRUE=0;
        end
        intern_counter=intern_counter+1;
    end
    
%     for w=1:Np
%         if ~fixed(w)
% %         disp(['Initial range of (' num2str(w) '): ' Par.Properties.RowNames(w) '--> ' num2str(Range(w,:))])
% %         disp(['Actual range --> ' num2str(New_range{i,w}')])
%         ParU.Range(w,:)=New_range(:,:,i+1);  
%         end
%     end
    ParU.Range=New_range(:,:,i+1);  
    if breaking == 1
        disp('------------------------------')
        disp(["Running method " num2str(method)])
        disp('------------------------------')
    end

if show    
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
    
    figure(2)
    clf
    for w=1:Np
        subplot(D1,D2,w)
        histogram(BFM_esc_best(:,w+1))
        xlabel(Par.Properties.RowNames(w));
        hold on
        xl1 = xline(New_range(w,1,i+1));
        xl2 = xline(New_range(w,2,i+1));
        xl1.LineWidth = 2;
        xl2.LineWidth = 2;
        xl1.Color = [1 0 0];
        xl2.Color = [1 0 0];
    end
    
    figure(3)
    clf
    if sup(i)>0
        gsua_plot('UncertaintyAnalysis',Par,Y_test(order,:,:),xdata,y_Nom,sup(i));
    else
        gsua_plot('UncertaintyAnalysis',Par,Y_test(order,:,:),xdata,y_Nom);
    end
    
    figure(4)
    clf
    plot(sup(1:i)/N)
    title('Progress')
%     save('Par_Final.mat','ParU');
end
end
end