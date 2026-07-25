function [] = gsua_mtest(name,T,simvals,xdata,gname,pause)
%GSUA_MTEST Animate model output as one parameter is swept, saved as a GIF
%
%   gsua_mtest(name,T,simvals,xdata,gname,pause)
%
%   Diagnostic/presentation utility: repeatedly overrides a single
%   parameter (T{name,'gaming'}) with each value in simvals, evaluates
%   the model against its nominal curve using GSUA_EVAL, captures each
%   frame, and (when simvals has more than one row) assembles the frames
%   into an animated GIF. Useful for visually exploring model
%   sensitivity to one parameter at a time.
%
%   Inputs:
%     name    <-- row name (char) of the parameter to sweep
%     T       <-- summary table from gsua_dataprep, with an
%                 Estlsqc/estimation column providing the nominal
%                 parameter set
%     simvals <-- Nframes x Np matrix of parameter sets to evaluate, one
%                 row per animation frame
%     xdata   <-- domain at which to evaluate the model
%     gname   <-- base file name (without extension) for the output GIF
%     pause   <-- delay (seconds) between GIF frames
%
%   Outputs: none (writes '<gname>.gif' to the current folder when
%   simvals has more than one row)
%
%   See also GSUA_EVAL, GSUA_DEVAL.
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
    sgtitle(num2str(i))
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
