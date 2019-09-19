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
if ~strcmp('mat',Table.Properties.CustomProperties.Kind)
    [y,~]=sens_montecarlo(Table,par,xdata,false,[0,size(par,1)]);
    Sname=false;
else
    fun=Table.Properties.CustomProperties.Solver;

try
    ni=Table.Properties.CustomProperties.NumVars;
    out=Table.Properties.CustomProperties.output;
    Sname=strcmp(Table.Properties.CustomProperties.Sname,'ode45');
    if nargin ==2
        try
        xdata=Table.Properties.CustomProperties.Domain(1):Table.Properties.CustomProperties.Domain(2);
        catch
            xdata=1;
        end
    end
catch
    T2=load('ATable.mat');
    T2=T2.Table2;
    ni=T2.NumVars;
    out=T2.output;
    Sname=strcmp(T2.Sname,'ode45');
    if nargin ==2
        xdata=T2.Domain(1):T2.Domain(2);
    end
end

nout=size(out,2);
if size(xdata,2)>1 && Sname
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
else
    if Sname
        y=fun(par);
    else
        y=gsua_intrp(fun(par(:,1:ni),par(:,ni+1:end)),xdata,out);
    end
end
end
if ~isempty(ydata)
    gsua_plot('UncertaintyAnalysis',Table,y,xdata,ydata)
else
    gsua_plot('UncertaintyAnalysis',Table,y,xdata)
end

if size(y,1)==1 && size(xdata,2)>1 && Sname
    y=y_temp;
end

end