classdef Imager
    % IMAGER a class specificying the properties of an imager
    %   
    % 2013-10-02 Matteo Carandini
    % 2014-02-28 DS added "HardwareBinning" property. moved under +tools directory
    % 2015-08-27 DS added "Info" property
    
    properties
        Type                 % e.g. 'PCO'
        DataDir              % the place where the data are kept
        FileString = '';     % eg '_cam2'
        MmPerPixel = [];     % info about the sensor
        Magnification = NaN; % this is AFTER hardware binning (but before software binning) useful for scale bars
        HardwareBinning = 1;       % hard ware binning factor. cf. GetCameraInfo % DS on 28.2.14
        Info = '';           % any other info specific to each camera. eg. scanmode
    end
    
    methods
        function I = Imager( Type, DataDir, FileString, Magnification, HardwareBinning )
%             I = Imager()
%             I = Imager( Type, DataDir )
%             I = Imager( Type, DataDir, FileString )
%             I = Imager( Type, DataDir, FileString, Magnification )
%             I = Imager( Type, DataDir, FileString, Magnification, HardwareBinning )
            if nargin < 5
                HardwareBinning = 1;
            end
            if nargin < 4
                Magnification = NaN;
            end
            if nargin < 3
                FileString = '';
            end
            
            if nargin < 1
                ImagerTypes = {'PCO','PhotonFocus','Dalsa','TwoPhoton'};
                [iType,ok] = listdlg(...
                    'ListString',ImagerTypes, ...
                    'SelectionMode','Single',...
                    'PromptString','Select imager type');
                if ~ok, return; end
                Type    = ImagerTypes{iType};
                
                %% TO DO: should use uicontrol.ratiobutton
                DataDir = uigetdir('\\zserver2.ioo.ucl.ac.uk\Data\','Choose a raw data directory (e.g. GCAMP)'); 
                
                answer = inputdlg({'Extra string to append?','Magnification','Hardware binning factor'},'',1,{'_cam2','NaN', '1'});
                if ~isempty(answer)
                    FileString = answer{1}; 
                    Magnification = str2double(answer{2});
                    HardwareBinning = str2double(answer{3});
                end
            end

            switch Type
                case {'PCO' 'PCO_tif'}
                    I.MmPerPixel = 1/(10 ^3 /6.5) * HardwareBinning / Magnification; % DS on 28.2.14             
                case 'PhotonFocus'
                    I.MmPerPixel = NaN; % please enter appropriate value
                case 'Dalsa'
                    I.MmPerPixel = NaN; % please enter appropriate value                    
                case 'TwoPhoton'
                    I.MmPerPixel = NaN; % please leave NaN
                otherwise
                    error('I don''t know this type of imager');
            end
            
            I.Type = Type;
            I.DataDir = DataDir; 
            I.FileString = FileString; 
            I.Magnification = Magnification; 
            I.HardwareBinning = HardwareBinning;
        end
    end
    
end

