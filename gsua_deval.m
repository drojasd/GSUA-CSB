function y = gsua_deval(par,Table2,xdata)
try
    Table=Table2.Properties.CustomProperties;
catch
    Table=load('ATable.mat');
    Table=Table.Table2;
end
kind=Table.Kind;   
  
if ~strcmp(kind,'mat')
    [y,~]=sens_montecarlo(Table2,par,xdata,false,[0,size(par,1)]);
else
ni=Table.NumVars;
fun=Table.Solver; 
out=Table.output;
Sname=strcmp(Table.Sname,'ode45');
try
    if size(xdata,2)>1 && Sname
        y=deval(fun(par(:,1:ni),par(:,ni+1:end)),xdata,out);
    else
        if Sname
            y=fun(par);
        else
            y=gsua_intrp(fun(par(:,1:ni),par(:,ni+1:end)),xdata,out);
        end
    end
catch ME
    warning(ME)
    y=inf(1,max(size(xdata)));
end
end
end