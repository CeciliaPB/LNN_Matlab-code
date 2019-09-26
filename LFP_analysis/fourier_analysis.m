function [FFT, fq] = fourier_analysis(signal, freq, varargin)

% Provides the spectrum density ( i.e. frequency) of the time-domain signal.
%
% INPUTS: 
%   - signal: signal to analyse, a vector
%   - freq: frequency range to analyse, an interval: [f1 f2]
%   Varargins
%   - fr: sampling rate. Use as 'fr',value. Default 30000.
%   - n:  downsampling factor. Use as 'n',value. Default 3. 
%   - mother: the mother wavelet function. Use as 'mother','name'. Default
%             Matlab's 'fft' with Butterworh filter 'Butter'. It can be
%             changed to 'Multi-taper'.
%
% OUTPUTS: 
%   - FFT: results of the fast fourier transform.
%   - fq: the frequency vector.
%   - The plot
%   If no output vars are provided, then only the plot.
% 
% Examples: 
%   [FFT, fq] = fourier_analysis(signal,[0.5 300]);
%   fourier_analysis(signal,[0.5 300],'mother','Multi-taper');
%   fourier_analysis(signal,[0.5 300]);
% 
% -------------------------------------------------------------------------
% Cecília Pardo-Bellver, 2019
% Laboratory of Network Neurophysiology
% Institute of Experimantal Medicine, Hungary.
%
% Multi-taper script was provided by Chronux see website
% http://chronux.org/ for more information.
%
% MATLAB toolboxes: Signal Processing Toolbox.
% -------------------------------------------------------------------------

% Optional input variables ------------------------------------------------
% Params, defalut values:
fr = 30000; 
n = 3;
mother = 'Butter';

if nargin
  for i=1:2:size(varargin,2)
    switch varargin{i}
      case 'fr'
        fr = varargin{i+1};
      case 'n' 
        n = varargin{i+1};
      case 'mother'
        mother = varargin{i+1};
      otherwise
        errordlg('Unknown varargin');
    end
  end
else
end
%--------------------------------------------------------------------------

switch mother
    case 'Butter'    
        signal2 = downsample(signal,n);
        
        wn = freq/(fr/2);   %bandpass
        [b,a] = butter(2,wn); % 2nd order  
        ft = filter(b,a,signal2);
        N = length(ft); % Lower values truncate the signal vector.
        f = linspace(0,fr/n,N);
        F = fft(ft,N);
        figure; 
        plot(f,smooth(abs(F),0.0001,'lowess')); % 0.0001
        xlim([min(freq) max(freq)]);

    case 'Multi-taper' % Multi-taper fourier transform - continuous data 
        signal2 = downsample(signal,n);
        
        DEF = input('Params are required. Use default params? Y/N [Y]: ','s');
        if isempty(DEF)
            DEF = 'Y';
        end
        switch DEF
            case {'Y','y'} 
                params.tapers = [4 12];                  
            case {'N','n'}
                params.tapers = input('Enter "tapers":'); % A 1x2 double
                if n == 3
                    n = input('Enter downsampling factor: [Default 3]'); % A 1x1 double
                    if isempty(n)
                        n = 3;
                    end
                elseif n ~= 3
                end
        end
        
        params.Fs = fr/n;
        params.fpass = freq;

        [F,f] = mtspectrumc(signal2,params);
        figure;
        plot(f,smooth(F,0.001,'lowess')); % Adjust smooth level (0.02)
        xlim([min(freq) max(freq)]);
end

if nargout > 0
    FFT = F;
    fq = f;
else
end

end