function convertOpenEphysToDat_2(filename,varargin)

% Convert OpenEphys .continuous files to ONE .dat file
%
% INPUTS: 
%   - filename: name of the file to save. STRING.
%   Varargin: 
%   - 'datadir': path name for reading data. Default current directory.
%   - 'resdir': path name for writing result files. Default current
%   directory.
%   - 'nChan': Number of channels. Default 32.
%   - 'processor': processor number. Default 101.
% 
% EXAMPLES
% convertOpenEphysToDat_2('continuous')
% convertOpenEphysToDat_2('continuous','nChan',64,'processor',100)
%
% -------------------------------------------------------------------------
% Cecília Pardo-Bellver, 2021
% Laboratory of Network Neurophysiology
% Institute of Experimental Medicine, Hungary.
%
% Based in convertOpenEphysToRawBInary from Kilosort2.
% Modification in tis version: the naming doesn't include the "_CH" part
% -------------------------------------------------------------------------

% Default arguments
prs = inputParser;
addOptional(prs,'datadir',cd,@(s)isempty(s)|isdir(s))  % data directory
addOptional(prs,'resdir','',@(s)isempty(s)|isdir(s))   % results directory
addOptional(prs,'nChan',32,@isnumeric)   % Number of channels (default: 32 channels)
addOptional(prs,'processor',101,@isnumeric) % Processor number, default 101
parse(prs,varargin{:})
params = prs.Results;

fname       = fullfile(params.datadir, sprintf('%s.dat', filename)); 
if ~isempty(params.resdir)
    fname       = fullfile(params.resdir, sprintf('%s.dat', filename));
end
fidout      = fopen(fname, 'w');
 
clear files
for jj = 1:params.nChan
   files{jj} = dir(fullfile(params.datadir, sprintf([num2str(params.processor) '_%d.continuous'], jj)));
end

nblocks = cellfun(@(x) numel(x), files);
if numel(unique(nblocks))>1
   error('different number of blocks for different channels!') 
end
 
nBlocks     = unique(nblocks);
nSamples    = 1024;  % fixed to 1024 for now!

fid = cell(params.nChan, 1);

for kk = 1:nBlocks
    for jj = 1:params.nChan
        fid{jj}             = fopen(fullfile(params.datadir, files{jj}(kk).name));
        % discard header information
        fseek(fid{jj}, 1024, 0);
    end
     
    nsamps = 0;
    flag = 1;
    
    while 1
        samples = zeros(nSamples * 1000, params.nChan, 'int16');
        for jj = 1:params.nChan
            collectSamps = zeros(nSamples * 1000, 1, 'int16'); 
            rawData      = fread(fid{jj}, 1000 * (nSamples + 6),...
                '1030*int16', 10, 'b');
            nbatches     = ceil(numel(rawData)/(nSamples+6));
            for s = 1:nbatches
                rawSamps = rawData((s-1) * (nSamples + 6) +6+ [1:nSamples]);
                collectSamps((s-1)*nSamples + [1:nSamples]) = rawSamps;
            end
            samples(:,jj)         = collectSamps;
        end
        
        if nbatches<1000
            flag = 0;
        end
        if flag==0
            samples = samples(1:s*nSamples, :);
        end
       
        samples         = samples';
        fwrite(fidout, samples, 'int16');
        
        nsamps = nsamps + size(samples,2);
        
        if flag==0
            break;
        end
    end
    params.nSamplesBlocks(kk) = nsamps;
    
    for jj = 1:params.nChan
       fclose(fid{jj}); 
    end
    
end
    
fclose(fidout);

end