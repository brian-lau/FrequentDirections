suite = matlab.unittest.TestSuite.fromFolder(pwd);
runner = matlab.unittest.TestRunner.withTextOutput;
dir = fileparts(which('FrequentDirections'));
%runner.addPlugin(matlab.unittest.plugins.CodeCoveragePlugin.forFolder(dir));
result = runner.run(suite);
table(result)