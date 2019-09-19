function [M,T] = gsua_dmatrix(Table,N,varargin)
% Function for design of experiments (factor space sampling)
%
% [M,T2]=gsua_dmatrix(T,N)
% Parameters:
% T <-- summary table from gsua_dataprep
% N <-- number of samples
% Outpus:
% M  <-- design matrix of NxNp for later routines
% T2 <-- summary table with fixed parameters actualized
% Additional features:
% You can choose a method for factor space sampling between uniform
% distribution design and latin hypercube design (default). To switch
% between methods use the paired feature 'Method' and 'Uniform' or
% 'LatinHypercube'.
% Also you can visualize the sampling result using the paired feature
% 'Show', 'on'
% M=gsua_dmatrix(T,N,'Method','Uniform','Show','on')
p=inputParser;
defaultShow='off';
validShow={'off' 'on'};
checkShow = @(x) any(validatestring(x,validShow));
defaultMethod='LatinHypercube';
validMethod={'LatinHypercube','Uniform'};
checkMethod = @(x) any(validatestring(x,validMethod));

addRequired(p,'Table');
addRequired(p,'N',@isnumeric);
addParameter(p,'Method',defaultMethod,checkMethod);
addParameter(p,'Show',defaultShow,checkShow);

parse(p,Table,N,varargin{:})
T=p.Results.Table;
N=p.Results.N;
method=p.Results.Method;
show=p.Results.Show;
Range=T.Range';
Np=size(T,1);
try
    Table2=T.Properties.CustomProperties;
catch
    TP=load('ATable.mat');
    Table2=TP.Table2;
end

if strcmp(Table2.rMethod,{'normal'})
    disp('Given ranges represent a normal distribution')
    M = zeros(N,Np);
    for k = 1:Np
        pdfun = makedist('normal','mu',Range(1,k),'sigma',Range(2,k));
        M(:,k) = random(pdfun,1,N); % Normal distribution between two values given by the range of parameters
    end
else
    switch method
        case 'LatinHypercube'
            M = lhsdesign(N,Np).*( ones(N,1)*(Range(2,:)-Range(1,:)) ) + ones(N,1)*Range(1,:);
        case 'Uniform'
            M = zeros(N,Np);
            for k = 1:Np
                pdfun = makedist('Uniform','lower',Range(1,k),'upper',Range(2,k));
                M(:,k) = random(pdfun,1,N); % uniform distribution between two values given by the range of parameters
            end
     end
end
try
    T.Properties.CustomProperties.Fixed=Range(1,:)==Range(2,:);
catch 
    Table2.Fixed=Range(1,:)==Range(2,:);
    save('ATable','Table2');
end

if strcmp(show,'on')
    figure
    clf
    gsua_plot('ScatterParameter',T,method,M)
end
end