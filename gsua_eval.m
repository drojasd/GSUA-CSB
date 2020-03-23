function y = gsua_eval(par,Table,xdata,ydata)
% Function for few model evaluations
% 
% Y=gsua_eval(values,Table)
% Parameters:
% values <-- array of NpxN (number of factors x number of simulations)
% Table  <-- summary table from gsua_dataprep function
% Outputs:
% Y <-- array with model output
% Additional features:
% If you provide an array (xdata) that belongs to the model domain, then
% the model output is interpolated to match xdata. Also you can provide a
% previous model output or experimental data as an array (ydata) that must
% coincide with the xdata array:
% Y=gsua_eval(values,Table,xdata,ydata)
if nargin<4
    ydata=[];
end
par=par';
try
    Table2=Table.Properties.CustomProperties;
catch
    Table2=load('ATable.mat');
    Table2=Table2.Table2;
end
kind=Table2.Kind; 


if kind==1
    [y,~]=sens_montecarlo(Table,par,xdata,false,[0,size(par,1)]);
else
    fun=Table2.Solver;
    ni=Table2.NumVars;
    out=Table2.output;
    if nargin ==2
        try
        xdata=Table2.Domain(1):Table2.Domain(2);
        catch
            xdata=1;
        end
    end
nout=size(out,2);
switch kind
    case 3
        y=zeros(size(par,1),size(xdata,2),nout);
        if size(par,1)>40
            parfor i=1:size(par,1)
                y_temp2=deval(fun(par(i,1:ni),par(i,ni+1:end)),xdata,out);
                for j=1:nout
                    y(i,:,j)=y_temp2(j,:);
                end
                clc
                disp(strcat('Sim ',num2str(i),' Done'))
            end
        else
            for i=1:size(par,1)
                y_temp=deval(fun(par(i,1:ni),par(i,ni+1:end)),xdata,out);
                for j=1:nout
                    y(i,:,j)=y_temp(j,:);
                end
            end
        end
    case 2
        y=fun(par);
    case 4
        y=gsua_intrp(fun(par(:,1:ni),par(:,ni+1:end)),xdata,out);
    case 5
        y=gsua_intrp(fun(par),xdata,out);
    case 6
        y=zeros(size(par,1),1);
        for i=1:size(par,1)
            y(i)=fun(par(i,:));
        end
end
end
if ~isempty(ydata)
    gsua_plot('UncertaintyAnalysis',Table,y,xdata,ydata)
else
    gsua_plot('UncertaintyAnalysis',Table,y,xdata)
end

if size(y,1)==1 && (kind==3)
    y=y_temp;
end

end