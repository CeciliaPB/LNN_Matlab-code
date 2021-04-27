function ops = convertOpenEphysToDat(filename,varargin)

% -------------------------------------------------------------------------
% Cec�lia Pardo-Bellver, 2021
% Laboratory of Network Neurophysiology
% Instinute of Experimental Medicine, Hungary.
%
% MATLAB toolboxes: - Statistics and Machine Learning Toolbox 
% -------------------------------------------------------------------------

% Default arguments
prs = inputParser;
addOptional(prs,'datadir',cd,@(s)isempty(s)|isdir(s))  % data directory
addOptional(prs,'resdir','',@(s)isempty(s)|isdir(s))   % results directory
addOptional(prs,'nChan',32,@isnumeric)   % Number of channels (default: 32 channels)
addOptional(prs,'processor',101,@isnumeric) % Processor number, default 101
addParameter(prs,'reference','common_avg',@(s)ischar(s)|isempty(s))   % switch for referencing
parse(prs,varargin{:})
ops = prs.Results;

fname       = fullfile(ops.datadir, sprintf('%s.dat', filename)); 
if ~isempty(ops.resdir)
    fname       = fullfile(ops.resdir, sprintf('%s.dat', filename));
end
fidout      = fopen(fname, 'w');
 
clear files
for jj = 1:ops.nChan
   files{jj} = dir(fullfile(ops.datadir, sprintf([num2str(ops.processor) '_CH%d.continuous'], jj)));
end

nblocks = cellfun(@(x) numel(x), files);
if numel(unique(nblocks))>1
   error('different number of blocks for different channels!') 
end
 
nBlocks     = unique(nblocks);
nSamples    = 1024;  % fixed to 1024 for now!

fid = cell(ops.nChan, 1);

tic
for kk = 1:nBlocks
    for jj = 1:ops.nChan
        fid{jj}             = fopen(fullfile(ops.datadir, files{jj}(kk).name));
        % discard header information
        fseek(fid{jj}, 1024, 0);
    end
     
    nsamps = 0;
    flag = 1;
    
    while 1
        samples = zeros(nSamples * 1000, ops.nChan, 'int16');
        for jj = 1:ops.nChan
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
    ops.nSamplesBlocks(kk) = nsamps;
    
    for jj = 1:ops.nChan
       fclose(fid{jj}); 
    end
    
end
    
fclose(fidout);

toc