classdef gsuaToolboxReleaseCompatTest < matlab.unittest.TestCase
    %GSUATOOLBOXRELEASECOMPATTEST The packaged .mltbx must have no upper release bound (v1.17.1 fix).
    %
    %   v1.17.0 shipped with releaseCompatibility end="R2025b" (from param.release.end="latest",
    %   which the .prj-derived ToolboxOptions resolved to the packaging release). Installing on a
    %   newer release or on MATLAB Online then reports "not supported for this MATLAB release".
    %   buildToolbox now leaves MaximumMatlabRelease empty -- the documented way to say "no maximum
    %   restriction". MATLAB serializes that as end="latest" (its unbounded sentinel); it has no way
    %   to physically omit the attribute, so the artifact must have NO CONCRETE pinned release
    %   (end="latest", not end="R2025b") while keeping the R2017a minimum and MATLAB Online. This
    %   builds a throwaway .mltbx and inspects its metadata/systemRequirements.xml.

    properties
        ToolboxRoot
        OutFile
    end

    methods (TestClassSetup)
        function locate(testCase)
            testDir = fileparts(mfilename('fullpath'));
            testCase.ToolboxRoot = fileparts(testDir);
        end
    end

    methods (TestMethodTeardown)
        function cleanOutput(testCase)
            if ~isempty(testCase.OutFile) && isfile(testCase.OutFile)
                delete(testCase.OutFile);
            end
        end
    end

    methods (Test)
        function artifactHasNoMaxReleaseButKeepsOnline(testCase)
            % packageToolbox reads GSUA-UCI.prj relative to pwd, so build from the toolbox root.
            oldDir = cd(testCase.ToolboxRoot); %#ok<NASGU>
            restore = onCleanup(@() cd(oldDir)); %#ok<NASGU>

            testCase.OutFile = [tempname '.mltbx'];
            buildToolbox('0.0.1', testCase.OutFile);   % throwaway version; must be Major.Minor.Bug
            testCase.assertTrue(isfile(testCase.OutFile), 'buildToolbox did not produce the .mltbx.');

            tmp = tempname; mkdir(tmp);
            unpack = onCleanup(@() rmdir(tmp,'s')); %#ok<NASGU>
            unzip(testCase.OutFile, tmp);
            xmlFile = fullfile(tmp,'metadata','systemRequirements.xml');
            testCase.assertTrue(isfile(xmlFile), 'No metadata/systemRequirements.xml in the package.');
            xml = fileread(xmlFile);

            relTok = regexp(xml, '<releaseCompatibility[^>]*', 'match', 'once');
            testCase.assertNotEmpty(relTok, 'No releaseCompatibility element found.');
            % The regression was a CONCRETE pinned upper release (end="R2025b"); it must be gone.
            % Empty MaximumMatlabRelease serializes to end="latest" (unbounded) or no end at all.
            testCase.verifyEmpty(regexp(relTok, 'end\s*=\s*"R\d{4}[ab]"', 'once'), ...
                sprintf('releaseCompatibility must not pin a concrete max release, got: %s', relTok));
            % Minimum release preserved.
            testCase.verifyNotEmpty(regexp(relTok, 'start\s*=\s*"R2017a"', 'once'), ...
                sprintf('releaseCompatibility must keep start="R2017a", got: %s', relTok));
            % MATLAB Online preserved.
            platTok = regexp(xml, '<platformCompatibility[^>]*', 'match', 'once');
            testCase.assertNotEmpty(platTok, 'No platformCompatibility element found.');
            testCase.verifyNotEmpty(regexp(platTok, 'MATLABOnline\s*=\s*"true"', 'once'), ...
                sprintf('platformCompatibility must keep MATLABOnline="true", got: %s', platTok));
        end
    end
end
