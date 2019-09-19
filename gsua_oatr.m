function ParT=gsua_oatr(ParT,varargin)
% Function for once-at-time ranges expantion
%
% T=gsua_oatr(T)
% Parameters:
% T <-- summary table from gsua_dataprep function. The objective output is
% given by T.Nominal.
% Outputs:
% T <-- summary table with new factor intervals in T.Range
% Additional features:
% To define limit for range expantion use the paired feature 'lim'.
% Default:0.3.
% To speed up range calculus, use the paired feature 'parallel',true.
% T=gsua_oatr(T,'lim',0.2,'parallel',true)

p=inputParser;
defaultLim=0.3;
defaultParallel=false;
addRequired(p,'ParT',@istable);
addParameter(p,'lim',defaultLim,@isnumeric);
addParameter(p,'parallel',defaultParallel,@islogical);

parse(p,ParT,varargin{:})
ParT=p.Results.ParT;
lim=p.Results.lim;
parallel=p.Results.parallel;
names=ParT.Properties.RowNames;
Range=ParT.Range;
noms=ParT.Nominal';
try
    TP=ParT.Properties.CustomProperties;
catch
    TP=load('ATable.mat');
    TP=TP.Table2;
end
kind=TP.Kind;
if strcmp('mat',kind)
domain=TP.Domain;
xdata=linspace(domain(1),domain(2),(domain(2)-domain(1)));
else
    time=TP.output;
    xdata=min(time):max(time);
end
y_Nom=gsua_deval(ParT.Nominal',ParT,xdata);
fixed=TP.Fixed;

disp("Expansion-reduction Method OAT is being launched")
disp('------------------------------------')
limit=sum((y_Nom-(y_Nom*(1+lim))).^2,2);
Np=size(ParT,1);
if isempty(fixed)
    fixed=false(1,Np);
end

Tf=mat2cell(ParT.Range,ones(1,1,Np));

%%
if parallel
    parforArg=Inf;
else
    parforArg=0;
end

parfor (i=1:Np,parforArg)
    if ~fixed(i)
    tnom=noms;
    intern_counter=0;
    aumenter=0.7;
    counter=zeros(1,2);
    flag=1;
    saver_down=1;
    saver_up=1;
    disp(['Calculus of ' names{i} ' Starting']) 
    while flag==1 || flag==2
        try
            tnom(i)=noms(i)*aumenter;
            y_temp =  gsua_deval(tnom,ParT,xdata);
         
        catch
            warning(['Parameter ' names{i} ' has caused a error'])
            if flag==1
                aumenter=(aumenter+saver_down)/2;
                %Tf{i}(flag)=noms(i)*saver_down;
            else
                aumenter=(aumenter+saver_up)/2;
                %Tf{i}(flag)=noms(i)*saver_up;
            end
            counter(1)=3;
            continue
%             aumenter=1.5;
%             counter=zeros(1,2);
%             flag=flag+1;
%             saver_down=1;
%             saver_up=1;
%             continue
        end
        J_temp=sum((y_Nom-y_temp).^2,2);
        if J_temp>limit*1.1
            act=1;
        elseif J_temp<limit
            act=2;
        else
            act=3;
        end
       
        switch act
            case 1
                if counter(1)==1 || counter(1)==0 
                    saver_down=aumenter;
                    counter(2)=counter(2)+1;
                    aumenter=aumenter^(1/(counter(2)+1));
                    counter(1)=1;

                else
                    saver_down=aumenter;
                    aumenter=(saver_up+aumenter)/2;
                    counter(1)=3;
                end
            case 2
                if counter(1)==2 || counter(1)==0
                    saver_up=aumenter;
                    counter(2)=counter(2)+1;
                    aumenter=aumenter^(counter(2)+1);
                    counter(1)=2;
                else
                    saver_up=aumenter;
                    aumenter=(saver_down+aumenter)/2;
                    counter(1)=3;
                end

            otherwise
                Tf{i}(flag)=noms(i)*aumenter;
                aumenter=1.5;
                counter=zeros(1,2);
                flag=flag+1;
                saver_down=1;
                saver_up=1;
    
        end
    
        intern_counter=intern_counter+1;
        if intern_counter>100
            disp(['Activating intern_counter for ' names{i}])
            Tf{i}(flag)=noms(i)*aumenter;
            aumenter=1.5;
            counter=zeros(1,2);
            flag=flag+1;
            saver_down=1;
            saver_up=1;
        end
    end
    disp(['Range for ' names{i} ' Done!'])
    disp(['Initial range of ' names{i} '(' num2str(i) ')--> ' num2str(Range(i,:))])
    disp(['Actual range --> ' num2str(Tf{i})])
    end
end
for i=1:Np
    if Tf{i}(1)>Tf{i}(2)
        Tf{i}=[Tf{i}(2),Tf{i}(1)];
    end
end
Normalized=zeros(2,Np);
    for h=1:Np
        Normalized(:,h)=[(Tf{h}(1)-Range(h,1))/(Range(h,2)-Range(h,1));...
            (Tf{h}(2)-Range(h,1))/(Range(h,2)-Range(h,1))];
    end
    ParT.Range=cell2mat(Tf);   
figure(1)
    clf
    boxplot(Normalized(:,~fixed),'Labels',names(~fixed)')
    
    title('Iteration #0')
    drawnow
end