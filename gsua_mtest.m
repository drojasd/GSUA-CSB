function [] = gsua_mtest(name,T,simvals,xdata,gname,pause)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
close all
T.Nominal=T.Estlsqc(:,1);
T.gaming=T.Nominal;
prop=T.Properties.CustomProperties;
disp('Getting Nominal curve ...')
y_nom=gsua_deval(T.Nominal',prop,xdata);
disp('Nominal curve done')
nvals=size(simvals,1);
for i=1:nvals
    T{name,'gaming'}=simvals(i,:)';
    fig=figure(1);
    gsua_eval(prop.Solver,T.gaming,T,xdata,y_nom);
    suptitle(num2str(i))
    im{i} = frame2im(getframe(fig));
    disp(strcat(num2str(i/nvals*100),'%'))
end
if nvals>1
    gname=strcat(gname,'.gif');
    n=1;
    [A,map] = rgb2ind(im{n},8);
    imwrite(A,map,gname,'gif','LoopCount',Inf,'DelayTime',pause);
    for n = 2:nvals
        [A,map] = rgb2ind(im{n},8);
        imwrite(A,map,gname,'gif','WriteMode','append','DelayTime',pause);
    end
end
disp('All work Done!') 
end
