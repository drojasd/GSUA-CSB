function [labels, interp] = gsua_labels(names)
%GSUA_LABELS Prepare parameter or output names for display and pick the text interpreter
%
%   [labels, interp] = gsua_labels(names)
%
%   Figures show names with MATLAB's default 'tex' interpreter, which understands TeX
%   markup such as '\beta_h' or '\mu_{w,l}' but not LaTeX: a name written as '$\phi$'
%   keeps its dollar signs, and commands such as \hat, \mathrm or \frac are printed
%   literally (with a "valid interpreter syntax" warning). Names copied from a LaTeX
%   document therefore rendered badly.
%
%   If any name looks like LaTeX (contains '$' or a LaTeX-only command), the whole set
%   switches to the 'latex' interpreter -- a legend or a tick-label set has a single
%   interpreter, so the names must agree -- and each name is converted: already in
%   $...$ is kept, TeX-style math (backslash, '_' or '^') is wrapped in $...$, and a
%   plain name is left as text. Otherwise the names come back unchanged with 'tex',
%   which is exactly the previous behavior.
%
%   Inputs:
%     names  <-- char, string array, or cell array of names
%
%   Outputs:
%     labels <-- the names ready to display: a char vector for a char input, a cell
%                array of the same size otherwise
%     interp <-- 'tex' or 'latex', to pass as the 'Interpreter' (or
%                'TickLabelInterpreter') of whatever displays the labels
%
%   See also GSUA_PLOT, GSUA_MCF, GSUA_IA.

wasChar = ischar(names);
c = cellstr(names);

% LaTeX commands that MATLAB's 'tex' interpreter does not know.
latexOnly = ['\\(hat|bar|tilde|dot|ddot|vec|check|breve|acute|grave|widehat|' ...
    'widetilde|overline|underline|mathrm|mathbf|mathit|mathsf|mathtt|mathcal|' ...
    'mathbb|boldsymbol|text|textrm|textit|textbf|frac|dfrac|tfrac|sqrt|' ...
    'operatorname|left|right)(?![a-zA-Z])'];
looksLatex = contains(c, '$') | ~cellfun(@isempty, regexp(c, latexOnly, 'once'));

if any(looksLatex(:))
    interp = 'latex';
    wrap = ~contains(c, '$') & ~cellfun(@isempty, regexp(c, '[\\_^]', 'once'));
    c(wrap) = cellfun(@(s) ['$' s '$'], c(wrap), 'UniformOutput', false);
else
    interp = 'tex';
end

if wasChar
    labels = c{1};
else
    labels = c;
end
end
