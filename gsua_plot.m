function []=gsua_plot(plot_type,p1,p2,p3,p4,p5)
% Visualization function for GSUA toolbox
%
% []=gsua_plot(plot_type,p1,p2,p3,p4,p5)
% Check gsua_userguide for function features (section 4.2)
% open gsua_userguide

colormap jet
Par = p1;
Np = size(Par,1);
try
   TP=Par.Properties.CustomProperties; 
catch 
    TP=load('ATable.mat');
    TP=TP.Table2;
end
fixed=TP.Fixed;
if isempty(fixed)
    fixed=false(1,Np);
end
out=TP.output;
Vars=TP.Vars;
try
    SensMethod=TP.SensMethod;
catch
    SensMethod='no one';
end
Np = size(Par,1)-sum(fixed);
if any(fixed)
    names=Par.Properties.RowNames(~fixed);
else
    names=Par.Properties.RowNames;
end
switch plot_type
    case 'UncertaintyAnalysis'
        Y = p2; t = p3; nout=size(Y,3);
        D1 = floor(sqrt(nout)); % Number of rows of subplot
        D2 = D1+ceil((nout-D1^2)/D1); % Number of columns of subplot
        n=size(Y,1);
        width=round(max(1,log10(n)/log10(30)));
        for j=1:nout
            subplot(D1,D2,j)
        switch nargin
            case 4
                if size(t,2)>1
                plot(t,Y(:,:,j),'b');
                else
                    boxplot(Y(:,:,j));
                end
            case 5
                y_nom = p4(j,:);
                if size(t,2)>1
                
                h = plot(t,Y(:,:,j),'b',t,y_nom,'r');
                set(h(size(Y,1)+1),'linewidth',width);
                else
                    boxplot(Y(:,:,j));
                    hold on
                    plot(1,y_nom,'o','linewidth',width);
                    legend('Nominal output')
                end
            case 6
                y_nom = p4;
                sup=p5;
                h1=plot(t,Y(sup+1:end,:),'Color',[ 0 0.75 1]);
                hold on
                h2=plot(t,Y(1:sup,:),'Color',[0 0.2 0.8]);
                h3=plot(t,y_nom,'Color',[1,0,0]);
                legend([h1(1),h2(end-1),h3(end)],{'Distant outputs','Right outputs','Nominal output'})
                set(h3,'linewidth',2);
            otherwise
                disp('Give the right number of 4 or 5 function inputs')
                disp('sens_plot(''UncertaintyAnalysis'',Par,Y,t)')
                disp('sens_plot(''UncertaintyAnalysis'',Par,Y,t,y_nom)')
                return
        end
        try
            title(Vars(out(j)))
        catch
            title(strcat('Output',num2str(j)))
        end
        xlabel ('Time')
        end
        suptitle(strcat('N= ',num2str(size(Y,1))));
    case 'FractionalSensitivityArea'
        if nargin ~=4
            disp('Give the right number of 4 function inputs')
            disp('sens_plot(''FractionalSensitivityArea'',Par,S,t)')
            return
        end
        S = p2; t = p3';
        Nd = size(t,1);
        t0 = 2; % The sensitivity in t = 0 is not computed
        while sum(isnan(S(:,t0))) > 0
            t0 = t0+1;
        end
        area(t(t0:Nd),S(:,t0:Nd)')
        xlabel('Time')
        ylabel('Si = Vi/V')
        title({['Time-dependent fractional sensitivity indices using ' SensMethod ' method']; ' '},'Color','r')
        legend(Par.Properties.RowNames{1:Np},'Location','BestOutside')
        colormap('colorcube')
        
    case 'FractionalSensitivityPlots'
        if nargin ~=4
            disp('Give the right number of 4 function inputs')
            disp('sens_plot(''FractionalSensitivityPlots'',Par,S,t)')
            return
        end
        S = p2; t = p3';
        Nd = size(t,1);
        D1 = floor(sqrt(Np)); % Number of rows of subplot
        D2 = D1+ceil((Np-D1^2)/D1); % Number of columns of subplot
        t0 = 2;
        while sum(isnan(S(:,t0)))>0
            t0 = t0+1;
        end
        for i = 1:Np
            subplot(D1,D2,i)
            plot(t(t0:Nd),S(i,t0:Nd)')
            xlabel('Time')
            ylabel('STi = Vi/V')
            title(Par.Properties.RowNames{i})
        end
        h = title(axes,{['Time-dependent fractional sensitivity indices using ' SensMethod ' method'];' '},'Color','r');
        set(gca,'visible','off')
        set(h,'visible','on')
        
    case 'TotalSensitivityArea'
        if nargin ~=4
            disp('Give the right number of 4 function inputs')
            disp('sens_plot(''TotalSensitivityArea'',Par,ST,t)')
            return
        end
        ST = p2; t = p3';
        Nd = size(t,1);
        t0 = 1;
        while sum(isnan(ST(:,t0)))>0
            t0 = t0+1;
        end
        [~,index]=sort(nansum(ST,2),'descend');
        ST=ST(index,:);
        names=names(index);
        for i = t0:Nd
            ST(:,i)=ST(:,i)./sum(ST(:,i));
        end
        area(t(t0:Nd),ST(:,t0:Nd)','facecolor','flat')
        xlabel('Time')
        ylabel('STi = VTi/V   STi=STi/sum(STi)')
        title({['Normalized vectorial sensitivity indices (' SensMethod ')']; ' '},'Color','r')
        legend(names{1:Np},'Location','BestOutside')
        colormap('colorcube')
        
    case 'TotalSensitivityPlots'
        if nargin ~=4
            disp('Give the right number of 4 function inputs')
            disp('sens_plot(''TotalSensitivityArea'',Par,ST, t)')
            return
        end
        ST = p2(~fixed,:); t = p3';
        D1 = floor(sqrt(Np)); % Number of rows of subplot
        D2 = D1+ceil((Np-D1^2)/D1); % Number of columns of subplot
        for i = 1:Np
            subplot(D1,D2,i)
            plot(t,ST(i,:)')
            xlabel('Time')
            ylabel('SiT = ViT/V')
            title(Par.Properties.RowNames{i})
        end
        h = title(axes,{['Time-dependent total sensitivity indices using ' SensMethod ' method'];' '},'Color','r');
        set(gca,'visible','off')
        set(h,'visible','on')
        
    case 'ScatterOutput'  % Scatter plots of Y vs Parameters in the time instant tref
        Y = p2; M = p3; t = p4; 
        if nargin ~= 6
            tref=1;
            c=1;
        else
            tref = p5(1);
            [c] = find(abs(t-tref)<=min(abs(t-tref))); % Find the column that correspons to a time instant
        end
        D1 = floor(sqrt(Np));
        D2 = D1+ceil((Np-D1^2)/D1);

        for i = 1:Np
            subplot(D1,D2,i)
            plot(M(:,i),Y(:,c),'.')
            xlabel(Par.Properties.RowNames{i})
            ylabel('Y')
        end
        h = title(axes,{['Scatterplots of Y vs Parameters in t = ' num2str(tref)];' '},'Color','r');
        set(gca,'visible','off')
        set(h,'visible','on')
        
    case 'ScatterParameter'  % Scatter plots of every pair of parameters
        if nargin ~=4
            disp('Give the right number of 4 function inputs')
            disp('sens_plot(''ScatterParameter'',Par,SampleMethod,M)')
            return
        end
        SampleMethod = p2; M = p3(:,~fixed);
        Table = array2table(M,'VariableNames',genvarname(Par.Properties.RowNames(~fixed)));
        sdo.scatterPlot(Table)
        h = title(axes,{['Scatterplot of pair of parameters with ' SampleMethod ' method'];' '},'Color','r');
        set(gca,'visible','off')
        set(h,'visible','on')
        
    case 'Pie'
        switch nargin
            case 5 % Time-dependent pie charts of sensitivity indices
                S = p2; t = p3; tref = p4;
                Nplot = size(tref,2);
                D1 = floor(sqrt(Nplot));
                D2 = D1+ceil((Nplot-D1^2)/D1);
                for j = 1:Nplot
                    % Find the column that correspons to a time instant
                    [c] = find( abs(t-tref(j)) <= min(abs(t-tref(j))) );
                    labels = cell(Np,1);
                    for i = 1:Np
                        if S(i,c)<=0
                            S(i,c) = 1e-6;
                        end
                        if S(i,c)*100/max([sum(S(:,c)),1]) <= 2
                            labels{i}='';
                        else
                            labels{i}=[Par.Properties.RowNames{i} ' (' num2str(S(i,c)*100/max([sum(S(:,c)),1]),3) '%)'];
                        end
                    end
                    subplot(D1,D2,j)
                    h1 = pie(S(:,c),labels);
                    hText = findobj(h1,'Type','text'); % text object handles
                    set(hText,'FontSize',7);
                    title({['Si ( t=' num2str(t(c)) ' )'];' '},'color','red')
                end % Pie chat
                legend(Par.Properties.RowNames{1:Np},'Location','BestOutside')
                h = title(axes,{['Time-dependent pie charts of sensitivity indices using ' SensMethod ' method'];''},'Color','r');
                set(gca,'visible','off')
                set(h,'visible','on')
            case {4,3} % The pie chart of sensitivity indices for MSE function
                if nargin==3
                    lim=90;
                else
                    lim=p3;
                end
                [S,index]= sort(p2,'descend');    
                labels = cell(Np,1);
                explode=zeros(Np,1);
                names=names(index);
                S(S<=0)=1e-6;
                S=S.*100/max([sum(S),1]);  
                others=Np-find(cumsum(S)<90,1,'last')-1;
                for i = 1:Np
                    if sum(S(1:i))-S(i)<lim
                        labels{i}=[names{i} ' (' num2str(S(i),3) '%)'];
                    else
                        explode(i:end)=1;
                        labels{i}=['\Sigma (' num2str(others) ' Factors) = '  num2str(sum(S(i:end)),3) '%'];
                        labels(i+1:end)=repmat(cellstr(''),size(labels(i+1:end)));

                        break      
                    end
                end
                h1 = pie(S,explode,labels);
                hText = findobj(h1,'Type','text'); % text object handles
                set(hText,'FontSize',8);
                title({['Escalar sensitivity indices (' SensMethod ')'];' '},'Color','r')
                legend(names{1:Np},'Location','BestOutside')
            otherwise
                disp('Give the right number of 3 or 5 function inputs')
                disp('sens_plot(''Pie'',Par,S,t,tref)')
                disp('sens_plot(''Pie'',Par,S)')
                return
        end
        colormap('colorcube')
        
    case 'Bar'
        switch nargin
            case 5 % Time-dependent bar charts of sensitivity indices
                [S,index] = sort(p2,'ascend'); t = p3; tref = p4;
                Nplot = size(tref,2);
                D1 = floor(sqrt(Nplot));
                D2 = D1+ceil((Nplot-D1^2)/D1);
                for j = 1:Nplot
                    % Find the column that correspons to a time instant
                    [~,c]  =find(abs(t-tref(j))<=min(abs(t-tref(j))));
                    if sum(S(:,c))<1
                        S1=[S(:,c); 1-sum(S(:,c))];
                    else
                        S1=[S(:,c); 0.00001];
                    end
                    subplot(D1,D2,j)
                    barh(S1)
                    set(gca,'ytick',1:Np,'yticklabel',Par.Properties.RowNames(index),'FontSize',6)
                    title(['Si ( t=' num2str(t(c)) ' )'])
                    for i = 1:Np
                        if S1(i) > 0.3
                            alignment = 'right';
                            colorl = 'white';
                        else
                            alignment = 'left';
                            colorl = 'blue';
                        end
                        text(S1(i),i,['  ' num2str(S1(i),3)],'HorizontalAlignment',alignment,'color',colorl,'FontSize',6)
                    end
                end
                h = title(axes,{['Time-dependent bar charts of sensitivity indices using ' SensMethod ' method'];' '},'Color','r');
                set(gca,'visible','off')
                set(h,'visible','on')
            case 3 % Bar chart of sensitivity indices for MSE function
                Par = p1; [S,index] = sort(p2,'ascend');
                %S=normalize(S,'scale',sum(S));
                Np = size(Par,1);
                barh(S);
                set(gca,'ytick',1:Np,'yticklabel',Par.Properties.RowNames(index),'FontSize',10,'FontWeight','Bold')
                title({['Bar chart of sensitivity indices for MSE function using ' SensMethod ' method'];' '},'Color','r')
                for i = 1:Np
                    if S(i) > 0.3
                        alignment = 'right';
                        colorl = 'white';
                    else
                        alignment = 'left';
                        colorl = 'blue';
                    end
                    text(S(i),i,['    ' num2str(S(i),3)],'HorizontalAlignment',alignment,'color',colorl,'FontSize',8)
                end
                    

            otherwise
                disp('Give the right number of 3 or 5 function inputs')
                disp('sens_plot(''Bar'',Par,S)')
                disp('sens_plot(''Bar'',Par,S,t,tref)')
                return
        end
        colormap('colorcube')
     
    case 'MC_filtering'    
         Par=p1; M=p2; Jsup=p3; Jinf=p4;
         Np=size(Par,1);
         D1 = floor(sqrt(Np)); % Number of rows of subplot
         D2 = D1+ceil((Np-D1^2)/D1); % Number of columns of subplot

         for k = 1:Np
             subplot(D1,D2,k);
             cdfplot(Jinf{k})
             hold on
             cdfplot(Jsup{k})
             cdfplot(M(:,k))
             xlabel('Value')
             ylabel('CDF')
             title(Par.Properties.RowNames{k})
             hold off
         end
         legend({'Low Values','High Values','Prior'})
         h = title(axes,{[SensMethod ' method'];' '},'Color','r');
         set(gca,'visible','off')
         set(h,'visible','on')
    otherwise
        disp('Select the right plot type')
        return
end
end