function info = read_polyscan_seq_pngs(seqFile)
%READ_POLYSCAN_SEQ_PNGS  Extract and resolve PNG files from PolyScan .seq
%
% info = READ_POLYSCAN_SEQ_PNGS(seqFile)
%
% Output struct fields:
%   info.pngNames        - PNG filenames in sequence order
%   info.found           - logical array (true if resolved)
%   info.fullPaths       - resolved full paths ('' if not found)
%   info.seqFolder       - folder containing the .seq file
%   info.missingNames    - filenames not found on disk
%
% Resolution strategy:
%   1) Same folder as .seq
%   2) Subfolders of .seq folder (recursive)

    arguments
        seqFile (1,:) char
    end

    % --- Resolve seq path
    seqFile = char(seqFile);
    if ~isfile(seqFile)
        error('SEQ file not found: %s', seqFile);
    end

    seqFolder = fileparts(seqFile);
    if isempty(seqFolder)
        seqFolder = pwd;
    end

    % --- Read raw binary
    fid = fopen(seqFile, 'rb');
    raw = fread(fid, inf, 'uint8=>char')';
    fclose(fid);

    % --- Extract PNG names (case-insensitive)
    expr = '(?i)[A-Za-z0-9_\-\.]+\.png';
    matches = regexp(raw, expr, 'match');

    if isempty(matches)
        warning('No PNG filenames found in %s', seqFile);
        info = struct([]);
        return
    end

    % --- Deduplicate while preserving order
    pngNames = {};
    for i = 1:numel(matches)
        if ~any(strcmpi(matches{i}, pngNames))
            pngNames{end+1} = matches{i}; %#ok<AGROW>
        end
    end

    % --- Attempt resolution
    n = numel(pngNames);
    found     = false(1, n);
    fullPaths = strings(1, n);

    % 1) Same directory as .seq
    for i = 1:n
        candidate = fullfile(seqFolder, pngNames{i});
        if isfile(candidate)
            found(i)     = true;
            fullPaths(i) = candidate;
        end
    end

    % 2) Recursive search under seq folder (only if needed)
    if ~all(found)
        allPngs = dir(fullfile(seqFolder, '**', '*.png'));
        for i = 1:n
            if found(i), continue; end
            idx = find(strcmpi(pngNames{i}, {allPngs.name}), 1);
            if ~isempty(idx)
                found(i)     = true;
                fullPaths(i) = fullfile(allPngs(idx).folder, allPngs(idx).name);
            end
        end
    end

    % --- Collect missing
    missingNames = pngNames(~found);

    % --- Package output
    info.pngNames     = pngNames;
    info.found        = found;
    info.fullPaths    = fullPaths;
    info.seqFolder    = seqFolder;
    info.missingNames = missingNames;

end
