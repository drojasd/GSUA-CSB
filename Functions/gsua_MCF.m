function [Jsup,Jinf] = gsua_MCF(T,M,Y,y_exp,tIn,tOut)
%GSUA_MCF Monte-Carlo filtering plots (behavioral vs non-behavioral parameter sets)
%
%   [Jsup,Jinf] = gsua_MCF(T,M,Y,y_exp)
%   [Jsup,Jinf] = gsua_MCF(T,M,Y,y_exp,tIn,tOut)
%
%   Splits a design/sample matrix M into "low" and "high" subsets based
%   on whether each run's scalar output is below/above the experimental
%   reference y_exp (or, if tIn/tOut are given, based on the output
%   value at a single time instant tIn), then plots the empirical CDF of
%   each parameter for the full sample and both subsets. Called by
%   GSUA_UA after an uncertainty analysis to visualize which parameter
%   regions drive high/low model output.
%
%   Inputs:
%     T     <-- summary table from gsua_dataprep (used for
%               CustomProperties.Fixed and parameter names)
%     M     <-- N x Np design/sample matrix (from gsua_dmatrix)
%     Y     <-- model output for each row of M (scalar per row, or
%               N x Nd if dimension>1 and tIn/tOut are not given)
%     y_exp <-- experimental/reference output to compare against
%     tIn   <-- (optional) single time instant to filter on
%     tOut  <-- (optional) time vector matching the columns of Y, used
%               to locate tIn
%
%   Outputs:
%     Jsup <-- 1 x Np cell array, parameter values for runs above y_exp
%     Jinf <-- 1 x Np cell array, parameter values for runs below y_exp
%
%   See also GSUA_UA, GSUA_DMATRIX.

try
     fixed=T.Properties.CustomProperties.Fixed;
 catch
     TP=load('ATable.mat');
     TP=TP.Table2;
     fixed=TP.Fixed;
end

if isempty(fixed)
 fixed=zeros(1,length(M(1,:)));
end

M=M(:,~fixed);
names=T.Properties.RowNames(~fixed);

[N,Np]=size(M);


% 
D1 = floor(sqrt(Np)); % Number of rows of subplot
D2 = D1+ceil((Np-D1^2)/D1); % Number of columns of subplot
Jinf=cell(Np,1);
Jsup=cell(Np,1);

dimension=size(Y,1);
try
    if nargin<5

        Jexp=sum(y_exp);
        if dimension>1
            J=sum(Y,2);
        else
            J=Y;
        end
        t_in=0;

        for i=1:Np
            Jinf{i}=M(J<Jexp,i);
            Jsup{i}=M(J>Jexp,i);
            if isempty(Jinf{i})
              Jinf{i}=min(M(:,i));  
            end
            if isempty(Jsup{i})
              Jsup{i}=max(M(:,i));  
            end
        end
    else
        t_in=find(tOut==tIn);
        Jexp=(y_exp(t_in));
        J=(Y(:,t_in));
        for i=1:Np
            Jinf{i}=M(J<Jexp,i);
            Jsup{i}=M(J>Jexp,i);
        end
        if isempty(Jinf{i})
          Jinf{i}=min(M(:,i));  
        end
        if isempty(Jsup{i})
          Jsup{i}=max(M(:,i));  
        end
    end
        
for k = 1:Np
             subplot(D1,D2,k);
             cdfplot(M(:,k))
             hold on
             cdfplot(Jinf{k})
             cdfplot(Jsup{k})
             xlabel('Value')
             ylabel('CDF')
             title(names{k})
             hold off
end
         legend({'Prior','Low Values','High Values'})
         if t_in==0
             
            h = title(axes,{['Montecarlo Filtering for escalar Y with N: ' num2str(N)];' '},'Color','r');
         else
            h = title(axes,{['Montecarlo Filtering with N: ' num2str(N) ' in t: ' num2str(tIn)];' '},'Color','r');
         end
         set(gca,'visible','off')
         set(h,'visible','on')
catch ME
    disp(ME.message)
    return
end

end