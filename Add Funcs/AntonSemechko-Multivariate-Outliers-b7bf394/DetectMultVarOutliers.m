function [mu,S,RD,chi_crt]=DetectMultVarOutliers(X,n_out,id_in,vis)
% Detect outliers contained in a normally distributed multivariate dataset
% using procedure described in [1].
%
% INPUT:
%   - X       : N-by-d data matrix, where N>=(d+1) is the number of
%               samples/observations and d>1 is their dimensionality.
%               All entries in X must be real numbers.
%   - n_out   : optional argument specifying an upper limit on the 
%               number of outliers contained in X; 1 <= n_out <= (N-d-1).
%               Default setting is n_out = N - h, where h=floor((N+d+1)/2); 
%               same as in [1]. Lower values of n_out produce more 
%               accurate estimates of robust Mahalanobis distance 
%               (see RD below) at the expense of robustness.
%   - id_in   : optional N-by-1 binary vector indicating the subset of
%               observations in X that are known a priori to be valid 
%               (i.e., inliers). id_in=false(N,1) is the default 
%               setting; meaning that no prior knowledge about
%               the dataset is available. If there are at least d+1 
%               valid (and unique) observations, they will be used to 
%               obtain an initial estimate of the location and dispersion 
%               of the Gaussian distribution from which the samples were 
%               generated. Otherwise, these model parameters will be
%               initialized with a procedure similar to the one described 
%               in [1]; the difference being that the geometric median 
%               is used here to obtain a robust estimate of central 
%               tendency instead of coordinate-wise medians; the former is
%               rotation invariant where as the latter isn't.
%   - vis     : set vis=false to suppress visualization of the scatter 
%               plot of robust distances. vis=true is the default setting.                
%
% OUTPUT:
%   - mu      : 1-by-d vector specifying robust estimate of central
%               tendency of X (i.e., mean)
%   - S       : d-by-d matrix specifying robust estimate of dispersion
%               of X (i.e., covariance)
%   - RD      : N-by-1 array of robust Mahalanobis distances for 
%               observations in X. See [1] for details.
%   - chi_crt : 1-by-3 vector of critical values of inv_chi^2(1-alpha/2,d) 
%               for significance values (i.e., alpha) 0.2, 0.1, 0.05, and
%               0.01. For example, outlyingness of the i-th observation at
%               0.05 significance level can be tested as RD(i)>=chi_crt(3).
%               As stated in [1], these tests should be taken with a grain 
%               of salt and inspected visually (see 'vis' option) since RD 
%               follows chi^2 distribution only approximately. 
%
% REFERENCES:
% [1] Hadi, A.S., (1992), 'Identifying multiple outliers in multivariate 
%     data', Journal of the Royal Statistical Society. Series B 
%     (Methodological), Vol. 54, No. 3, pp. 761-771
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%


% Basic error checking
% -------------------------------------------------------------------------
if nargin<1 || isempty(X)
    error('Insufficient number of input arguments')
end

if ~isnumeric(X) || ~ismatrix(X)|| sum(isnan(X(:)) | isinf(X(:)))>0
    error('Incorrect entry for 1st input argument (X)')
end

[N,d]=size(X);
if N<(d+1)
    error('Too few samples (%u) relative to dimensionality of the data (%u)',N,d)
end

Y=unique(X,'rows','stable');
N2=size(Y,1);
if N2~=N && N2<(d+1) 
    error('Dataset contains only %u unique observations, which are too few relative dimensionality of the data (%u)',N2,d)
end
clear Y N2

if nargin<2 || isempty(n_out)
    n_out=N-floor((N+d+1)/2);
elseif ~isnumeric(n_out) || numel(n_out)~=1 || n_out<1 || n_out>(N-d-1) || n_out~=round(n_out)
    error('Incorrect entry for 2nd input argument (k)')
end

if nargin<3 || isempty(id_in)
    id_in=false(N,1);
elseif ~islogical(id_in) || numel(id_in)~=N
    error('Incorrect entry for 3rd input argument (id_val)')
end
id_in=id_in(:);

flag_in=false;
if sum(id_in)>0, flag_in=true; end

if nargin<4 || isempty(vis)
    vis=false;
elseif ~(isnumeric(vis) || islogical(vis)) || numel(vis)~=1
    error('Incorrect entry for 4th input argument (vis)')
end


% Obtain initial estimate of location and dispersion
% -------------------------------------------------------------------------
X_in=unique(X(id_in,:),'rows');
n_in=size(X_in,1);
if n_in<(d+1) 
    
    % Geometric median
    mu=GeometricMedian(X);
    
    % Covariance
    S=bsxfun(@minus,X,mu);
    S=(S'*S)/(N-1);        
    
    % Mahalanobis distance
    MD2=MahDistSquared(X,mu,S);
    
    % Subset of floor((N+d+1)/2) samples with smallest MD
    if flag_in, MD2(id_in)=0; end  % bias selection towards inliers; if they are known
    [~,id_srt]=sort(MD2,'ascend');
    ho=floor((N+d+1)/2);
    X_sub=X(id_srt(1:ho),:);
    
    % Initial estimates of location and dispersion
    mu=mean(X_sub,1);
    S=bsxfun(@minus,X_sub,mu);
    S=(S'*S)/(ho-1);
    
    n_in=d+1;
else
    mu=mean(X(id_in,:),1);
    S=bsxfun(@minus,X(id_in,:),mu);
    S=(S'*S)/(sum(id_in)-1);
end

% Refine estimates of location and dispersion; unless N_in>=h, in which 
% case can proceed to last step 
% -------------------------------------------------------------------------
h=N-n_out;
if n_in<h
    
    % Find h inliers
    Sd=det(S); 
    chk_det=Sd>eps;
    
    i=n_in-1;
    while i<h
        
        i=i+1;
        
        % Mahalanobis distance 
        MD2=MahDistSquared(X,mu,S);
        
        % Subset of i samples with smallest MD
        if flag_in, MD2(id_in)=0; end
        [~,id_srt]=sort(MD2,'ascend');
        X_sub=X(id_srt(1:i),:);
        
        % Updated estimates of location and dispersion
        mu_i=mean(X_sub,1);
        Si=bsxfun(@minus,X_sub,mu_i);
        Si=(Si'*Si)/(i-1); 
        
        if chk_det
            Sdi=det(Si);
            if Sdi<eps
                j=i;
                while Sdi<eps && j<h 
                    j=j+1;
                    X_sub=X(id_srt(1:j),:);
                    mu_i=mean(X_sub,1);
                    Si=bsxfun(@minus,X_sub,mu_i);
                    Si=(Si'*Si)/(j-1);
                    Sdi=det(Si);
                end
                i=j; %#ok<*FXSET>
            end
            chk_det=false;
        end        
        mu=mu_i;
        S=Si;        
        
    end
    
    n_in=h; 
     
    id_in=false(N,1);
    id_in(id_srt(1:h))=true;
end


% Test of outlyingness based on robust distances
% -------------------------------------------------------------------------

% Mahalanobis distances
MD=MahDistSquared(X,mu,S);
MD=sqrt(MD);
MD_srt=sort(MD,'ascend');

% Robust distances
m_j=MD_srt(n_in);
c_npr=(1+n_in/(N-d))^2;
c_b=(c_npr*m_j)/chi2inv(0.5,d);
RD=MD/sqrt(c_b);

% Critical values of chi^2 with d DOF
p_crt=[0.2 0.1 0.05 0.01];
chi_crt=sqrt(chi2inv(1-p_crt/2,d));


% Visualize robust Mahalanobis distances
% -------------------------------------------------------------------------
if ~vis, return; end

hf=figure('color','w');
drawnow
pause(0.1)
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
j_hf=get(hf,'JavaFrame'); 
j_hf.setMaximized(true);
warning('on','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
drawnow
pause(0.1)

spl_id=1:N;
h_a=plot(spl_id(id_in),RD(id_in)','ok','MarkerSize',6,'MarkerFaceColor','w'); % basis subset
hold on
h_b=plot(spl_id(~id_in),RD(~id_in)','ok','MarkerSize',6,'MarkerFaceColor','k'); % complement of the basis subset

if N<=10
    dx=0.25;
elseif N<=100
    dx=0.5;
elseif N<=250
    dx=1;
elseif N<=500
    dx=2;
elseif N<=1E3
    dx=4;
else
    dx=10;
end

h1=plot([1-dx N+dx],chi_crt(1)*[1 1],':','Color',[0.75 0 0.75],'LineWidth',1);
h2=plot([1-dx N+dx],chi_crt(2)*[1 1],':','Color',[0.75 0 0],'LineWidth',1);
h3=plot([1-dx N+dx],chi_crt(3)*[1 1],':','Color',[0 0.75 0],'LineWidth',1);
h4=plot([1-dx N+dx],chi_crt(4)*[1 1],':','Color',[0 0 0.75],'LineWidth',1);

xl=sprintf('sample id (1 to %u)',N);
xlabel(xl,'FontSize',25,'FontWeight','bold','Color','k')
ylabel('Robust Mahalanobis Distance','FontSize',25,'FontWeight','bold','Color','k')
set(gca,'FontSize',20,'XLim',[0 N+1],'XColor','k','YColor','k')

xt=get(gca,'XTick');
xt(1)=1;
xt=unique(xt);
set(gca,'XLim',[1-dx N+dx],'XTick',xt,'TickDir','out')

h=legend([h_a h_b h1 h2 h3 h4],...
         {'basis subset' 'basis complement' '$$\alpha=0.2$$' '$$\alpha=0.1$$' '$$\alpha=0.05$$' '$$\alpha=0.01$$'},...
          'Location','EastOutside',...
          'Orientation','vertical',...
          'Interpreter','latex');

p=get(h,'Position');
x=1-1.02*p(3); p(1)=x; set(h,'Position',p)

drawnow
pause(0.1)

p=get(gca,'Position'); w=0.98*(x-p(1)); p(3)=w; set(gca,'Position',p)

if nargout<1
    clear mu S RD chi_crt
end


% =========================================================================
function MD2=MahDistSquared(X,mu,S)

% Residuals
dX=bsxfun(@minus,X,mu);
d=size(dX,2);

% Covariance
if d==1
    dX=dX/(std(dX)+1E-8);
else
    
    % SVD of empirical covariance 
    [U,L,~]=svd(S);
    L=diag(L);
    L=sqrt(L(:)'); 

    % Regularization
    if min(L)<=1E-8
        id_cut=find(L>1E-8,1,'last');
        if ~isempty(id_cut)
            L=max(L,L(id_cut));
        end
    end
    
    % Project residuals on U and scale
    dX=dX*U;
    dX=bsxfun(@rdivide,dX,L);
    
end

% Mahalanobis distance (squared)
MD2=sum(dX.^2,2);

