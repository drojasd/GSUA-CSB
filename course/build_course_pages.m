function build_course_pages(modulos)
%BUILD_COURSE_PAGES Render the course's MATLAB demos to docs/course/*.html.
%
%   build_course_pages()
%   build_course_pages({'00','01'})   % solo esos módulos
%
%   Runs each module's MATLAB demo with PUBLISH and writes its rendered HTML
%   (sectioned code, console output and figures) into docs/course/, which
%   GitHub Pages serves. The demos are plain .m scripts with %% sections, so
%   PUBLISH is the right tool (EXPORT only accepts Live Scripts).
%
%   Figure images are emitted next to the HTML as <demo>_NN.png and referenced
%   relatively, so both the .html and its .png files belong in docs/course/.
%
%   This file lives in course/, which toolbox.ignore excludes, so it is not
%   shipped in the packaged toolbox.
%
%   See also PUBLISH, BUILDEXAMPLEPAGES.

here   = fileparts(mfilename('fullpath'));      % .../course
root   = fileparts(here);                        % repo root
outDir = fullfile(root, 'docs', 'course');
if ~isfolder(outDir), mkdir(outDir); end

demos = { ...
    '00_Escribir_Modelos',             'demo_00_escribir_modelos', 'modulo-00-matlab.html'
    '01_Introduccion_Modelos',         'demo_01_introduccion',  'modulo-01-matlab.html'
    '02_Muestreo',                     'demo_02_muestreo',      'modulo-02-matlab.html'
    '03_Incertidumbre_MonteCarlo',     'demo_03_incertidumbre', 'modulo-03-matlab.html'
    '04_Analisis_Sensibilidad',        'demo_04_sensibilidad',  'modulo-04-matlab.html'
    '05_Estimacion_Identificabilidad', 'demo_05_estimacion',    'modulo-05-matlab.html' };

if nargin >= 1
    demos = demos(startsWith(demos(:,1), string(modulos)), :);
end

opts = struct('format','html', 'outputDir',outDir, 'showCode',true, ...
              'evalCode',true, 'maxWidth',620, 'catchError',false);

for k = 1:size(demos,1)
    src = fullfile(here, demos{k,1}, 'matlab', [demos{k,2} '.m']);
    dst = fullfile(outDir, demos{k,3});
    fprintf('building %s ...\n', demos{k,3});
    t0 = tic;
    addpath(fullfile(here, demos{k,1}, 'matlab'));   % publish necesita el archivo en el path
    htmlOut = publish(src, opts);           % outDir/<demo>.html (+ figure PNGs)
    if isfile(dst), delete(dst); end
    movefile(htmlOut, dst);                  % rename to modulo-0X-matlab.html
    close all;
    fprintf('   done in %.1fs (%.0f KB)\n', toc(t0), dir(dst).bytes/1024);
end
end
