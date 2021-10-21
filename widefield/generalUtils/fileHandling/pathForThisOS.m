function pathOut = pathForThisOS(pathIn)
% pathOut = pathForThisOS(pathIn)
% returns correct path depending on OS
%

pathOut = pathIn;
pathOut(pathOut=='/') = filesep;
pathOut(pathOut=='\') = filesep;