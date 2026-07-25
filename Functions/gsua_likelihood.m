function Tfinal = gsua_likelihood(T,xdata,ydata,alpha,step,margin,tolerance1,tolerance2,limit,reps,show,parallel,saver,pars)
%GSUA_LIKELIHOOD Profile-likelihood confidence intervals for each parameter
%
%   Tfinal = gsua_likelihood(T,xdata,ydata,alpha,step,margin,tolerance1,tolerance2,limit,reps,show,parallel,saver,pars)
%
%   Computes a profile-likelihood-based confidence interval for each
%   parameter in T: starting from the nominal estimate, each parameter
%   is stepped away (independently, up and down) while T is
%   re-estimated (via GSUA_PE) with that parameter fixed at each trial
%   value, until the resulting likelihood-ratio statistic crosses the
%   chi-squared threshold implied by alpha. This is an alternative to
%   the sampling-based ranges used by GSUA_OATR/GSUA_CSB, based directly
%   on the likelihood-ratio test.
%
%   Inputs:
%     T          <-- summary table from gsua_dataprep (T.Nominal is the
%                    reference parameter set)
%     xdata      <-- domain at which the model is evaluated
%     ydata      <-- experimental/reference output at xdata (if empty,
%                    the nominal model output is used instead)
%     alpha      <-- confidence level for the chi-squared threshold
%                    (e.g. 0.95)
%     step       <-- initial relative step size used to perturb each
%                    parameter away from its nominal value
%     margin     <-- relative standard deviation used in the Gaussian
%                    likelihood weighting (see gsua_likecost)
%     tolerance1 <-- convergence tolerance on the likelihood-ratio
%                    statistic (bisection stopping criterion)
%     tolerance2 <-- convergence tolerance on the step size (relative to
%                    the initial nominal-to-bound distance)
%     limit      <-- maximum number of bisection iterations per
%                    parameter bound
%     reps       <-- number of estimation runs per GSUA_PE call while
%                    profiling
%     show       <-- (logical) plot each parameter's likelihood profile
%                    and corresponding model fit as it is computed
%     parallel   <-- (logical) let the inner GSUA_PE estimations use
%                    parallel (fmincon UseParallel) optimization
%     saver      <-- (logical) append each parameter's result to
%                    'likelihoodCI.mat' as it is computed
%     pars       <-- indices of the parameters to profile (empty = all)
%
%   Outputs:
%     Tfinal <-- copy of T with Tfinal.Range replaced by the
%                profile-likelihood interval bounds [lower,upper] per
%                parameter
%
%   See also GSUA_PE, GSUA_LIKECOST, GSUA_OATR, GSUA_CSB.
    margin2=margin-1;
    if isempty(pars)
        npars = size(T,1);
        pars=1:npars;
    else
        npars=length(pars);
    end
    lambda = chi2inv(alpha,1)/2;
    %lambda = chi2inv(alpha,npars)/2;
    if isempty(ydata)
        ydata=gsua_deval(T.Nominal',T,xdata);
    end
    ymodel = gsua_deval(T.Nominal',T,xdata);
    desv = (ydata*margin2).^2;
    cost = sum(log(2*pi*desv) +(ydata-ymodel).^2./desv,'all','omitnan')/2;
    thresh = cost+lambda;
    
    Tfinal=T;
    init = [0,npars*2];
    gsua_timer(0,init);

%    if show
%         %D1 = floor(sqrt(npars)); % Number of rows of subplot
%         %D2 = D1+ceil((npars-D1^2)/D1); % Number of columns of subplot
%         f1 = figure('Name','Likelihood profiles');
%    end
    Taux=zeros(npars,2); 
%    if parallel
%         show = false;
%         parforArg=Inf;
%         parallel2=false;
%     else
%         parforArg=0;
%    end

   if parallel
       opt = optimoptions('fmincon','UseParallel',true);
   else
       opt = optimoptions('fmincon','UseParallel',false);
   end

    %parfor (i=1:npars,parforArg)
    for i=pars
        
        for j=1:2
            counter2=1;
            Testim = T;
            stepPar = abs(Testim.Nominal(i)-Testim.Range(i,j))*step;
            if stepPar == 0
                stepPar=abs(Testim.Nominal(i)*step);
            end
            if j==1
                Testim.Nominal(i)=Testim.Nominal(i)-stepPar;
            else
                Testim.Nominal(i)=Testim.Nominal(i)+stepPar;
            end
            par = [T.Nominal(i), Testim.Nominal(i)];
            biflag=false;
            obj = [];

            while true
                
                Testim.Range(i,:)=[Testim.Nominal(i), Testim.Nominal(i)];
                [Tprueba,res] = gsua_pe(Testim,xdata,ydata,'margin',-margin,'ipoint',repmat(Testim.Nominal,1,reps)','solver','fmincon','N',reps,'save',false,'timer',false,'opt',opt);
                obj =  [obj, thresh-res(1)];

                if show
                    f1=figure(i);
                    %set(f1, 'WindowStyle', 'Docked');
                    f1.Name=append('Likelihood profile for ',T.Properties.RowNames{i});
                    subplot(1,2,j)
                    scatter(par(2:end-1),obj(1:end-1),'filled')
                    hold on
                    scatter(par(end),obj(end),'filled')
                    hold off
                    xline(T.Nominal(i))
                    xlim([min(par),max(par)])
                    yline(tolerance1)
                    yline(-tolerance1)
                    if j==1
                        title(append('Lower ',T.Properties.RowNames{i}))
                    else
                        title(append('Upper ',T.Properties.RowNames{i}))
                    end
                    drawnow
                    figure(i+1)
                    gsua_eval(Tprueba.Estfmincon,Tprueba,xdata,ydata);
                    drawnow
                end

               if Testim.Nominal(i)>T.Range(i,2) || Testim.Nominal(i)<T.Range(i,1)
                    Taux(i,j)=T.Range(i,j);
                    break
               end

                if abs(obj(end))<tolerance1 %terminación por tolerancia
                     Taux(i,j)=Tprueba.Nominal(i);
                     %Tfinal.Range(i,j)=Tprueba.Nominal(i);
                     break
                end

                if ~biflag
                    if obj(end)>0
                        if j==1
                            Testim.Nominal(i)=Testim.Nominal(i)-stepPar;
                        else
                            Testim.Nominal(i)=Testim.Nominal(i)+stepPar;
                        end
                        par = [par, Testim.Nominal(i)];
                        
                    else
                        Testim.Nominal(i)=(par(end)+par(end-1))/2;
                        par = [par, Testim.Nominal(i)];
                        biflag=true;
                        bp = [par(end-2),par(end-1)];%[a,b]
                    end
                else
                    counter2=counter2+1;
                    if obj(end)*obj(end-1)>0
                        bp = [bp(1),par(end)];
                    else
                        bp = [par(end-1),par(end)];
                    end
                    Testim.Nominal(i)=(bp(1)+bp(2))/2;
                    par = [par, Testim.Nominal(i)];
                end

                if counter2>limit %terminación por límite de interacciones
                     Taux(i,j)=Tprueba.Nominal(i);
                     break
                end
                if abs(par(end)-par(end-1))<tolerance2*abs(T.Nominal(i)-T.Range(i,j))%acciones en caso de iteraciones atascadas
                     Taux(i,j)=Tprueba.Nominal(i);
                     break
                end
            end
            if saver
                save_for_parfor('likelihoodCI.mat',1,Taux);
            end
            clc
            gsua_timer(1,[(2*i-(2-j)),init(2)]);
        end
    end
    Tfinal.Range=Taux;

end

function save_for_parfor(fname,numvars,varargin)
%SAVE_FOR_PARFOR Save parfor-loop variables to a .mat file by their caller-side names
%
%   save_for_parfor(fname,numvars,var1,...)
%
%   Local helper for gsua_likelihood: the built-in SAVE cannot be called
%   directly with a variable inside a parfor loop body, so this
%   recovers each input's original name via INPUTNAME and saves it
%   under that name instead.
%
%   Inputs:
%     fname   <-- output .mat file name
%     numvars <-- number of trailing variables to save
%     varargin <-- the variables to save, passed through so INPUTNAME
%                  can recover their caller-side names
    for i = 1:numvars
       eval([inputname(i+2),' = varargin{i};']);  
    end
    save('-mat',fname,inputname(3));
    for i = 2:numvars    
        save('-mat',fname,inputname(i+2),'-append');
    end
end
   