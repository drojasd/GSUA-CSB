function outFile = buildToolbox(version)
%BUILDTOOLBOX Package GSUA-CSB into release/<version>/GSUA-CSB.mltbx.
%
%   outFile = buildToolbox()          % version from GSUA-UCI.prj
%   outFile = buildToolbox("1.12")    % explicit version
%
%   WHY THIS EXISTS. matlab.addons.toolbox.ToolboxOptions('GSUA-UCI.prj') does not apply the
%   <fileset.rootfiles> whitelist or the <param.exclude.filters> blacklist that the .prj
%   declares -- it hands back every file under the project folder. Packaging straight from it
%   therefore ships whatever happens to be on disk. That is not hypothetical: v1.10 shipped at
%   15.5 MB carrying releases 1.7/1.8/1.9, the whole python/ port with __pycache__ and PEtab
%   test data, the MATLAB test suite, the duplicate "- copia" image folder, and a local editor
%   config. Because each release lands in release/, the next build then swallowed the previous
%   release too, so the artifact roughly doubled every version (2.5 MB -> 15.5 MB -> 31 MB).
%
%   This function applies what the .prj already says, reading both lists FROM the .prj rather
%   than restating them here, so a layout change stays correct without editing this file.
%
%   It also trims ToolboxMatlabPath to the shipped folders. Left alone, that metadata names
%   directories the package no longer contains and MATLAB emits an addpath warning per missing
%   folder on every user's install.
%
%   See also MATLAB.ADDONS.TOOLBOX.PACKAGETOOLBOX, MATLAB.ADDONS.INSTALL.

prj = 'GSUA-UCI.prj';
if ~isfile(prj)
    error('buildToolbox:NoProject','%s not found; run from the toolbox root.',prj);
end
txt = fileread(prj);

if nargin < 1 || isempty(version)
    version = regexp(txt,'<param\.version>(.*?)</param\.version>','tokens','once');
    if isempty(version)
        error('buildToolbox:NoVersion','Could not read <param.version> from %s.',prj);
    end
    version = string(version{1});
end
version = string(version);

opts = matlab.addons.toolbox.ToolboxOptions(prj);
root = string(opts.ToolboxFolder);

[keepRoots, fileRoots] = local_rootfiles(txt);
excludes = local_excludeRegexps(txt);

files = string(opts.ToolboxFiles);
keepF = local_select(erase(files,root), keepRoots, fileRoots, excludes);
opts.ToolboxFiles = cellstr(files(keepF));

% The toolbox root itself is dropped from the path metadata deliberately: it holds no
% functions, and MATLAB rewrites that absolute entry into a literal "D__\..." subfolder of the
% install directory, warning about it on every install.
mp = string(opts.ToolboxMatlabPath);
relM = erase(mp, root);
keepM = strlength(relM) > 0 & local_select(relM + filesep, keepRoots, fileRoots, excludes);
opts.ToolboxMatlabPath = cellstr(mp(keepM));

opts.ToolboxVersion = version;   % set explicitly: the .prj parse reports a stale value

% Validate BEFORE touching anything on disk, so a rejected build leaves the previous
% artifact intact rather than deleting it on the way to failing.
local_assertClean(erase(string(opts.ToolboxFiles), root));

releaseDir = fullfile('release', char(version));
if ~isfolder(releaseDir), mkdir(releaseDir); end
outFile = fullfile(releaseDir, 'GSUA-CSB.mltbx');
if isfile(outFile), delete(outFile); end
opts.OutputFile = outFile;

matlab.addons.toolbox.packageToolbox(opts);

info = dir(outFile);
fprintf('Packaged  : %s (%.1f KB)\n', outFile, info.bytes/1024);
fprintf('Version   : %s\n', opts.ToolboxVersion);
fprintf('Identifier: %s   (stable across versions -- never change)\n', opts.Identifier);
fprintf('Files     : %u of %u offered by the .prj parse\n', numel(opts.ToolboxFiles), numel(files));
fprintf('Path      : %u folder(s)\n', numel(opts.ToolboxMatlabPath));
end

function [dirRoots, fileRoots] = local_rootfiles(txt)
%LOCAL_ROOTFILES The <fileset.rootfiles> whitelist: what the .prj says the toolbox consists of.
block = regexp(txt,'<fileset\.rootfiles>(.*?)</fileset\.rootfiles>','tokens','once');
if isempty(block)
    error('buildToolbox:NoRootfiles','Could not read <fileset.rootfiles> from the .prj.');
end
entries = regexp(block{1},'<file>(.*?)</file>','tokens');
entries = string(cellfun(@(c) c{1}, entries, 'UniformOutput', false));
rel = strtrim(erase(entries,'${PROJECT_ROOT}'));
rel = regexprep(rel,'^[\\/]+','');
rel = strrep(rel,'/',filesep);

isDir = arrayfun(@(r) isfolder(fullfile(pwd,r)), rel);
dirRoots  = rel(isDir);
fileRoots = rel(~isDir);
end

function res = local_excludeRegexps(txt)
%LOCAL_EXCLUDEREGEXPS The <param.exclude.filters> blacklist, as anchored regexps.
block = regexp(txt,'<param\.exclude\.filters>(.*?)</param\.exclude\.filters>','tokens','once');
res = string.empty(1,0);
if isempty(block), return, end
lines = strtrim(string(splitlines(block{1})));
lines = lines(strlength(lines) > 0 & ~startsWith(lines,'%'));
for L = lines(:)'
    res(end+1) = local_globToRegexp(L); %#ok<AGROW>
end
end

function rx = local_globToRegexp(pat)
%LOCAL_GLOBTOREGEXP Translate one .prj exclude pattern into an anchored regexp.
%   "**/x" matches x at any depth; "*" matches within a single path segment; a bare name
%   matches that file or anything beneath that folder.
p = strrep(strtrim(pat),'/','\');
anyDepth = startsWith(p,'**\');
if anyDepth, p = extractAfter(p,3); end
q = regexptranslate('escape',p);
q = strrep(q,'\*','[^\\]*');          % '*' escaped by the line above, now made a segment wildcard
if anyDepth
    rx = "^(.*\\)?" + q + "(\\.*)?$";
else
    rx = "^" + q + "(\\.*)?$";
end
end

function tf = local_select(rel, dirRoots, fileRoots, excludes)
%LOCAL_SELECT Keep paths under a declared root, minus anything the exclude filters name.
rel = rel(:);
tf = false(size(rel));
for r = fileRoots(:)'
    tf = tf | rel == r | rel == (r + filesep);
end
for r = dirRoots(:)'
    % Match a whole path segment: a bare prefix test would let "doc" swallow "docs\".
    tf = tf | startsWith(rel, r + filesep);
end
for rx = excludes(:)'
    tf = tf & cellfun(@isempty, regexpi(rel, rx, 'once'));
end
end

function local_assertClean(rel)
%LOCAL_ASSERTCLEAN Fail the build loudly if a known contaminant survived the filter.
sep = string(filesep);
bad = ["release"+sep, "python"+sep, ".claude", "__pycache__", ".pytest_cache", ...
       "tests"+sep, "docs"+sep, " - copia", ".MATLABDriveTag", ".mltbx"];
for b = bad
    hit = rel(contains(rel, b));
    if ~isempty(hit)
        error('buildToolbox:Contaminated', ...
            'Refusing to package: %u file(s) matching "%s" survived the filter, e.g. %s', ...
            numel(hit), b, hit(1));
    end
end
end
