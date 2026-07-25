clc
fprintf(2,'HELP:\n')
disp('     The models must to be configured as specified in userguide (open gsua_userguide) ')
disp('     Parameters name is an array as follows: {''p1'',''p2'',...,''pk''}')
disp('     Parameters Range is a numerical array as follows [inf1,sup1; inf2,sup2; ... ; infk,supK] ')
disp('     Odes, if required, is an array with the symbolic system of equations  ')
disp('     Vars, if required, is an array with the symbolic states of the system ')
disp('     output, if required, is a scalar with the number of the output state of the system')
disp('     domain, if required, is a numerical array of 1x2 for model domain as follows: [inf,sup].')
disp('     This tool needs the Statistics Toolbox')
disp(' ')
% Selection of routine
routine = input('Select toolbox routine [PE (1), UA(2), SA(3), UCI(4)]: ');
user = input('Do you want to load and run an example (1) or a user-defined model (2)?: ');

% Selection of Simulink model and parameters
if user==1
disp('Simulink examples:''SIR'', ''Pendulum'', ''PID''')
disp('Symbolic Matlab example: ''sSIR''')
model_case = input('Type the name of the example: ');
switch model_case
    case 'SIR'
        model = 'sens_example_sir_sim';
        ParIn = {'beta','alpha','I0'};
        p=30;
        Ranges=[0.15 p; 0.45 p; 0.1 p];   
        y_exp = [];
        kind='simulink';
    case 'Pendulum'
        ParIn= {'m','l','g','f'};
        Ranges=[2.4 20;1.3 20;9.78,1;0.69,20];
        model = 'sens_example_pendulum_sim';
        y_exp = [];
        kind='simulink';
    case 'PID'
        ParIn={'b2','a1','a2'};
        Ranges=[3.2 20;2 20;4 20];
        model = 'sens_example_PID_sim';
        y_exp = [];
        kind='simulink';
    case 'sSIR'
        syms S(t) I(t) R(t) P(t) beta gamma P %symbolic state variables and parameters
        P=S+I+R;
        %Defining the system of differential equations
        ode1 = diff(S) == -S*I*beta/P;
        ode2 = diff(I) == S*I*beta/P - gamma*I;
        ode3 = diff(R) == gamma*I;
        %Array with the system
        odes=[ode1; ode2; ode3];
        vars = [R S I];
        domain=[0 264];
        modelName='sSIR';
        Ranges=[0 0; 100 3000; 1 500; 0 1; 0 1];
        kind='symbolic';
        time=linspace(domain(1),domain(2),1+domain(2)-domain(1));
    otherwise
       disp('Wrong example selection')
       return
end
if strcmp(kind,'simulink')    
    [T,time] = gsua_dataprep(model,Ranges,ParIn,'rMethod','percent');
else
    [T,solver] = gsua_dataprep(odes,vars,domain,modelName,'range',Ranges,'output',3);
end
else
    if user==2
        kind=input('Do work with simulink (1) or symbolic Matlab (2)?: ');
        if kind==1
            warning('Remember that your model must be well configured (Check userguide in: Working with simulink models')
            model = input('Give model name (archive): ');
            ParIn = input('Give the parameters name: ');
            Ranges = input('Give the parameters range: ');
            y_exp = input('Give the name of experimental time-response vector ('''' or set empty [] for using nominal time response): ');
            [T,time] = gsua_dataprep(model,Ranges,ParIn,'rMethod','percent');
        else
            warning('Remember that your model must be well configured (Check userguide in: Working with symbolic Matlab models')
            modelName = input('Give a name for the model: ');
            odes = input('Give the system of equations: ');
            vars = input('Give the system variables: ');
            domain = input('Give the system domain: ');
            [~] = gsua_dataprep(odes,vars,domain,modelName);
            Ranges = input('Give Range for parameters: ');
            output = input('Give the objective output [1,2,..,model_order]: ');
            [T,solver] = gsua_dataprep(odes,vars,domain,modelName,'range',Ranges,'output',output);
            time=linspace(domain(1),domain(2),1+domain(2)-domain(1));
        end
    else
        disp('Wrong option')
        return
    end
end
switch routine
    case 1
        y_est=input('Give an experimental output for optimization. If empty [], the toolbox will provide a random output: ');
        if isempty(y_est)
            y_est=gsua_deval(gsua_dmatrix(T,1),T,time);  
            noise= normrnd(0,0.3,size(time));
            y_est=y_est+y_est.*noise;%create a noisy output
        end
        solver=input('Select a solver for optimization: lsqcurvefit (1), genetic algorithm (2), particle swarm (3) ');
        parallel=input('Do use parallel? Yes: true, No: false ');
        switch solver %options for solvers
            case 1
                solver='lsqc';
                opt=optimoptions('lsqcurvefit','UseParallel',parallel,'Display','iter');
            case 2
                solver='ga';
                opt=optimoptions('ga','UseParallel',parallel,'Display','iter');
            case 3
                solver='particle';
                opt=optimoptions('particleswarm','UseParallel',parallel,'Display','iter');
        end
        N=input('Give number of optimization processes (if >9 an identifiability analysis is also performed): ');
        [T,res]=gsua_pe(T,time,y_est,'solver',solver,'Show','on','opt',opt,'N',N);
        if N>=10
            T_newRange=gsua_ia(T,T{:,3})
        end
    case {2,3}
        % Selection of sample method
        SampleMethod = input('Select the sample method (''Uniform'', ''LatinHypercube''): ');
        % Selection of sample size
        N = input('Sample size: ');
        % Selection or not of parallel computing
        parallel=input('Do use parallel computing? (Yes:true, No: false): ');
        %Table to data management
        M=gsua_dmatrix(T,N,'Method',SampleMethod,'Show','on');
        switch routine
            case 2
                Y=gsua_ua(M,T,'parallel',parallel);
            case 3
                % Selection of sensitivity method
                SensMethod = input('Select the sensitivity method (''brute-force'', ''Sobol'', ''Jansen'', ''Saltelli'', ''OAT'', ''Distance''): ');
                [T,J,Y] = gsua_sa(M,T,'SensMethod',SensMethod,'ynom',y_exp,'parallel',parallel);
                % Plots of relevant information of sensitivity analysis
                figure(2)
                gsua_plot('UncertaintyAnalysis',T,Y,time)
                figure(3)
                subplot(2,2,1), gsua_plot('FractionalSensitivityArea',T,SensMethod,T.STi_vec,time)
                subplot(2,2,2), gsua_plot('TotalSensitivityArea',T,SensMethod,T.STi_vec,time)
                subplot(2,2,3), gsua_plot('Pie',T,SensMethod,T.STi)
                subplot(2,2,4), gsua_plot('Bar',T,SensMethod,T.STi)
                figure(4)
                gsua_plot('ScatterOutput',T,Y,M,time,8)
                figure(5)
                gsua_plot('Pie',T, SensMethod,T.STi_vec,time,[1,3,5,7,10])
                figure(6)
                gsua_plot('Bar',T, SensMethod,T.STi_vec,time,[1,3,5,7,10])
                figure(7)
                gsua_plot('FractionalSensitivityPlots',T,SensMethod,T.Si_vec,time)
        end  
    case 4
        nominal=input('Give a set of parameters for UCI interval calculus (Give empty [] to calculus over a random set) ');
        if ~isempty(nominal)
            T.Nominal=nominal;
        else
            T.Nominal=gsua_dmatrix(T,1)';
        end
        parallel=input('Do use parallel computing? (Yes:true, No: false): ');
        T_oat=gsua_oatr(T);
        T_uci=gsua_uci(T_oat,100,'parallel',parallel)
    otherwise
        disp('Wrong routine selection')
        return
end
