syms S(t) I(t) R(t) P(t) P beta gamma %symbolic state variables and parameters
P=S+I+R;
%Defining the system of differential equations
ode1 = diff(S) == -S*I*beta/P;
ode2 = diff(I) == S*I*beta/P - gamma*I;
ode3 = diff(R) == gamma*I;
%Array with the system
odes=[ode1; ode2; ode3];
vars = [R S I];
[T,~] = gsua_dataprep(odes,vars,[0 264],'SIR','range',[0 0; 100 3000; 1 500; 0 1; 0 1],'output',3);