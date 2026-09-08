function buildExamplePages()
%BUILDEXAMPLEPAGES Export the MATLAB example Live Scripts to docs/examples/*.html.
%
%   buildExamplePages()
%
%   Run from the repository root. Executes each example and writes its rendered HTML
%   into docs/examples/, which GitHub Pages serves.
%
%   It also patches one thing the exporter gets wrong for a web page. MATLAB writes
%   every figure as <img ... style="width: 100%">, so a 560x337 px plot is stretched
%   to whatever the container is -- on a wide screen that upscales it two or three
%   times, which is why the axis labels and legends look enormous and slightly
%   blurred. The patch caps figures at their natural width and centres them; the
%   small inline images MATLAB uses for LaTeX are left alone, since they carry a
%   vertical-align style instead and must stay on the text baseline.
%
%   This file lives in docs/, which the .prj excludes, so it is not shipped in the
%   packaged toolbox.
%
%   See also EXPORT, BUILDTOOLBOX.

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
outDir = fullfile(here,'examples');
if ~isfolder(outDir), mkdir(outDir); end

examples = { 'pk_user_defined',  'pk-user-defined-matlab.html'
             'sir_symbolic',     'sir-symbolic-matlab.html' };

addpath(genpath(root));
for k = 1:size(examples,1)
    src = fullfile(root,'Examples',[examples{k,1} '.m']);
    dst = fullfile(outDir, examples{k,2});
    fprintf('building %s ...\n', examples{k,2});
    t0 = tic;
    if isfile(dst), delete(dst); end
    export(src, dst, 'Run', true);
    capFigureWidth(dst);
    fprintf('   done in %.1fs (%.0f KB)\n', toc(t0), dir(dst).bytes/1024);
end
end

function capFigureWidth(htmlFile)
%CAPFIGUREWIDTH Stop exported figures from being upscaled to the container width.
txt = fileread(htmlFile);
css = [ ...
    '<style>' ...
    'img[style*="width: 100%"]{max-width:620px!important;width:auto!important;' ...
    'height:auto;display:block;margin:18px auto;}' ...
    'body{max-width:900px;margin:0 auto;padding:0 20px;}' ...
    '</style>'];
% Insert immediately before </head> so it overrides the exporter's own rules.
idx = strfind(txt,'</head>');
if isempty(idx)
    warning('buildExamplePages:NoHead','No </head> in %s; figure width not capped.',htmlFile);
    return
end
txt = [txt(1:idx(1)-1) css txt(idx(1):end)];
fid = fopen(htmlFile,'w','n','UTF-8');
fwrite(fid, unicode2native(txt,'UTF-8'));
fclose(fid);
end
