%% Load ThorSync Episode HDF5 file: 
% - all will be exported to base workspace with unique names.
% - Attributes can be any combination of the following: 
%   [syncDataOut=LoadSyncEpisode(fnam)      -- Load the data from the h5 file at the fully qualified path 'fnam'
%                                           -- syncDataOut - Structure containing all of the data sets in the Thorsync dataset
%   Note: This version does not allow opening of only part of the sync episode. Any other arguments besides the file name are ignored.
% This is modified for use as a function from LoadSyncFunction.m (c) 2014
% --- Copyright (c) 2018, Thorlabs, Inc. All rights reserved. ---

%% Copyright (c) 2018, Thorlabs, Inc. All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the distribution
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function [syncDataOut]=LoadSyncEpisodeFunction(fnam, lineNames)
%[syncDataOut]=LoadSyncEpisodeFunction(fnam)
% loads all Thorsync data in a .h5 file specified as fnam
%
% [syncDataOut]=LoadSyncEpisodeFunction(fnam, lineNames)
% loads a subset of Thorsync data specified in lineNames
%
%
% 17/7/20 DS added 2nd input
%% Select one h5 file:

% [filename, pathname] = uigetfile('*.h5', 'Pick a Sync Episode file');
% if isequal(filename,0) || isequal(pathname,0)
   % disp('User pressed cancel')
   % return
% else
   % disp(['User selected ', fullfile(pathname, filename)])
% end

if nargin < 2
    lineNames = [];
end

[pathname,name,ext]=fileparts(fnam);
filename=[name ext];
pathname=[pathname filesep];

%% Load params from XML:
clockRate = 20000000;
sampleRate = LoadSyncXML(pathname);

%% Start loading HDF5:
pathandfilename = strcat(pathname,filename);
info = h5info(pathandfilename);

%% Parse input:
props = {'start','length','interval'};
data = {[1,1],[1 Inf],[1 1]};

% if(~isempty(varargin))
    % assert(rem(length(varargin),2)==0 && iscellstr(varargin(1:2:end)), 'Inputs failed to conform to expected string-value pair format, eg> ''start'',1');
    % %foundProps = intersect(varargin(1:2:end), props); 
    % IdxCell = cellfun(@(x) strcmpi(x,props),varargin(1:2:end),'UniformOutput',false);
    % val = double(cell2mat(varargin(2:2:end)))*sampleRate;
    % for i=1:length(val)
        % data{cell2mat(IdxCell(i))>0} = [1 val(i)];
    % end
% end

%% Read HDF5:

for j=1:length(info.Groups)
    for k = 1:length(info.Groups(j).Datasets)
        datasetPath = strcat(info.Groups(j).Name,'/',info.Groups(j).Datasets(k).Name);
        datasetName = info.Groups(j).Datasets(k).Name;
        
        if ~(strcmp(info.Groups(j).Name,'/Global')) && ~isempty(lineNames) ...
                && ~ismember(info.Groups(j).Datasets(k).Name, lineNames)
            continue;
        end
        datasetName(isspace(datasetName))='_';   
        datasetValue = h5read(pathandfilename,datasetPath,data{1},data{2},data{3})';
        %h5read is inherently slow...particularly when reading over network
        
        % load digital line in binary:
        if(strcmp(info.Groups(j).Name,'/DI'))
            datasetValue(datasetValue>0) = 1;
        end
        % create time variable out of gCtr, 
        % account for 20MHz sample rate:
        if(strcmp(info.Groups(j).Name,'/Global'))
            datasetValue = double(datasetValue)./clockRate;
            datasetName = 'time';
        end
        assignStr = UniqueName(datasetName);
		syncDataOut.(assignStr)=datasetValue;
        % assignin('base',assignStr,datasetValue);
    end
end

end


function outStr = UniqueName(str)
%% Generate unique name for variable to be exported.

cellfind = @(string)(@(cell_contents)(strcmp(string,cell_contents)));
vars = evalin('base','who');
index = 1;
unique = false;
cmpStr = str;

while (~unique)
    ret = cellfun(cellfind(cmpStr),vars);
    if(~any(ret))
        outStr = cmpStr;
        unique = true;
    else
        cmpStr = strcat(str,num2str(index,'%03d'));
        index=index+1;
    end
end

end

%% Load ThorSync XML file: 
% - input file path, output sample rate.
% --- Copyright (c) 2014, Thorlabs, Inc. All rights reserved. ---

%% Copyright (c) 2014, Thorlabs, Inc. All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the distribution
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function sampleRate = LoadSyncXML(varargin)
%% Load dll:
xmlFile = strcat(varargin{1}, 'ThorRealTimeDataSettings.xml');
assert(exist(xmlFile,'file')>0,'ThorRealTimeDataSettings.xml was not found. ');
dataStruct = xml2struct(xmlFile);

if(~isempty(dataStruct))
    BrdID = cellfun(@(x) strcmpi(x.Attributes.active,'1'),dataStruct.RealTimeDataSettings.DaqDevices.AcquireBoard);
    sampleID = cellfun(@(x) strcmpi(x.Attributes.enable,'1'),dataStruct.RealTimeDataSettings.DaqDevices.AcquireBoard{BrdID}.SampleRate);
    sampleRate = dataStruct.RealTimeDataSettings.DaqDevices.AcquireBoard{BrdID}.SampleRate{sampleID>0}.Attributes.rate;
end

end


%% Convert xml file into a MATLAB structure
% [ s ] = xml2struct( file )
%
% A file containing:
% <XMLname attrib1="Some value">
%   <Element>Some text</Element>
%   <DifferentElement attrib2="2">Some more text</Element>
%   <DifferentElement attrib3="2" attrib4="1">Even more text</DifferentElement>
% </XMLname>
%
% Will produce:
% s.XMLname.Attributes.attrib1 = "Some value";
% s.XMLname.Element.Text = "Some text";
% s.XMLname.DifferentElement{1}.Attributes.attrib2 = "2";
% s.XMLname.DifferentElement{1}.Text = "Some more text";
% s.XMLname.DifferentElement{2}.Attributes.attrib3 = "2";
% s.XMLname.DifferentElement{2}.Attributes.attrib4 = "1";
% s.XMLname.DifferentElement{2}.Text = "Even more text";
%
% Please note that the following characters are substituted
% '-' by '_dash_', ':' by '_colon_' and '.' by '_dot_'
%
% Written by W. Falkena, ASTI, TUDelft, 21-08-2010
% Attribute parsing speed increased by 40% by A. Wanner, 14-6-2011
% Added CDATA support by I. Smirnov, 20-3-2012
%
% Modified by X. Mo, University of Wisconsin, 12-5-2012

%% Copyright (c) 2010, Wouter Falkena All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the distribution
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function [ s ] = xml2struct( file )
%% Function:
    if (nargin < 1)
        clc;
        help xml2struct
        return
    end
    
    if isa(file, 'org.apache.xerces.dom.DeferredDocumentImpl') || isa(file, 'org.apache.xerces.dom.DeferredElementImpl')
        % input is a java xml object
        xDoc = file;
    else
        %check for existance
        if (exist(file,'file') == 0)
            %Perhaps the xml extension was omitted from the file name. Add the
            %extension and try again.
            if (isempty(strfind(file,'.xml')))
                file = [file '.xml'];
            end
            
            if (exist(file,'file') == 0)
                error(['The file ' file ' could not be found']);
            end
        end
        %read the xml file
        xDoc = xmlread(file);
    end
    
    %parse xDoc into a MATLAB structure
    s = parseChildNodes(xDoc);
    
end

% ----- Subfunction parseChildNodes -----
function [children,ptext,textflag] = parseChildNodes(theNode)
    % Recurse over node children.
    children = struct;
    ptext = struct; textflag = 'Text';
    if hasChildNodes(theNode)
        childNodes = getChildNodes(theNode);
        numChildNodes = getLength(childNodes);

        for count = 1:numChildNodes
            theChild = item(childNodes,count-1);
            [text,name,attr,childs,textflag] = getNodeData(theChild);
            
            if (~strcmp(name,'#text') && ~strcmp(name,'#comment') && ~strcmp(name,'#cdata_dash_section'))
                %XML allows the same elements to be defined multiple times,
                %put each in a different cell
                if (isfield(children,name))
                    if (~iscell(children.(name)))
                        %put existsing element into cell format
                        children.(name) = {children.(name)};
                    end
                    index = length(children.(name))+1;
                    %add new element
                    children.(name){index} = childs;
                    if(~isempty(fieldnames(text)))
                        children.(name){index} = text; 
                    end
                    if(~isempty(attr)) 
                        children.(name){index}.('Attributes') = attr; 
                    end
                else
                    %add previously unknown (new) element to the structure
                    children.(name) = childs;
                    if(~isempty(text) && ~isempty(fieldnames(text)))
                        children.(name) = text; 
                    end
                    if(~isempty(attr)) 
                        children.(name).('Attributes') = attr; 
                    end
                end
            else
                ptextflag = 'Text';
                if (strcmp(name, '#cdata_dash_section'))
                    ptextflag = 'CDATA';
                elseif (strcmp(name, '#comment'))
                    ptextflag = 'Comment';
                end
                
                %this is the text in an element (i.e., the parentNode) 
                if (~isempty(regexprep(text.(textflag),'[\s]*','')))
                    if (~isfield(ptext,ptextflag) || isempty(ptext.(ptextflag)))
                        ptext.(ptextflag) = text.(textflag);
                    else
                        %what to do when element data is as follows:
                        %<element>Text <!--Comment--> More text</element>
                        
                        %put the text in different cells:
                        % if (~iscell(ptext)) ptext = {ptext}; end
                        % ptext{length(ptext)+1} = text;
                        
                        %just append the text
                        ptext.(ptextflag) = [ptext.(ptextflag) text.(textflag)];
                    end
                end
            end
            
        end
    end
end

% ----- Subfunction getNodeData -----
function [text,name,attr,childs,textflag] = getNodeData(theNode)
    % Create structure of node info.
    
    %make sure name is allowed as structure name
    name = toCharArray(getNodeName(theNode))';
    name = strrep(name, '-', '_dash_');
    name = strrep(name, ':', '_colon_');
    name = strrep(name, '.', '_dot_');

    attr = parseAttributes(theNode);
    if (isempty(fieldnames(attr))) 
        attr = []; 
    end
    
    %parse child nodes
    [childs,text,textflag] = parseChildNodes(theNode);
    
    if (isempty(fieldnames(childs)) && isempty(fieldnames(text)))
        %get the data of any childless nodes
        % faster than if any(strcmp(methods(theNode), 'getData'))
        % no need to try-catch (?)
        % faster than text = char(getData(theNode));
        text.(textflag) = toCharArray(getTextContent(theNode))';
    end
    
end

% ----- Subfunction parseAttributes -----
function attributes = parseAttributes(theNode)
    % Create attributes structure.

    attributes = struct;
    if hasAttributes(theNode)
       theAttributes = getAttributes(theNode);
       numAttributes = getLength(theAttributes);

       for count = 1:numAttributes
            %attrib = item(theAttributes,count-1);
            %attr_name = regexprep(char(getName(attrib)),'[-:.]','_');
            %attributes.(attr_name) = char(getValue(attrib));

            %Suggestion of Adrian Wanner
            str = toCharArray(toString(item(theAttributes,count-1)))';
            k = strfind(str,'='); 
            attr_name = str(1:(k(1)-1));
            attr_name = strrep(attr_name, '-', '_dash_');
            attr_name = strrep(attr_name, ':', '_colon_');
            attr_name = strrep(attr_name, '.', '_dot_');
            attributes.(attr_name) = str((k(1)+2):(end-1));
       end
    end
end

