classdef gsuaNoImageToolboxDepTest < matlab.unittest.TestCase
    %GSUANOIMAGETOOLBOXDEPTEST The shipped Functions must not depend on Image Processing Toolbox.
    %
    %   The toolbox used to pull in the Image Processing Toolbox purely for immse (plain-MSE
    %   objective in gsua_pe) and corr2 (correlation penalty in the cost/sensitivity functions).
    %   Both are now implemented in-repo (gsua_mse local to gsua_pe, gsua_corr2). This guard runs
    %   MATLAB's own dependency analysis over Functions/*.m and fails if Image Processing Toolbox
    %   reappears in the required-products list -- i.e. if someone reintroduces immse, corr2, or
    %   any other function from that product.

    methods (Test)
        function functionsDoNotRequireImageProcessingToolbox(testCase)
            testDir = fileparts(mfilename('fullpath'));
            fdir = fullfile(fileparts(testDir), 'Functions');
            d = dir(fullfile(fdir, '*.m'));
            names = fullfile({d.folder}, {d.name});
            [~, pList] = matlab.codetools.requiredFilesAndProducts(names);
            products = string({pList.Name});
            offenders = products(contains(products, 'Image Processing'));
            testCase.verifyEmpty(offenders, ...
                sprintf('Functions/ now require Image Processing Toolbox again: %s', ...
                        strjoin(offenders, ', ')));
        end
    end
end
