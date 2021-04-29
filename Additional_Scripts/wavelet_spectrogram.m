function [wavelet_signal]=wavelet_spectrogram (signal,params,band,varargin)
% ------------------------------------------------------------------------------------------
% Wavelet spectrogram
% Usage: function [wsignal]=wavelet_bandratio (signal,params,band,varargin)
%        Input:    signal, raw signal
%                  params, previous wavelet parameters 
%                  band, frequency band as [a b]
%        Options: 'nonplot', hidding plot option, plot by default
%                 'zscore',  zscore normalitzation, i.e., mean 0 and
%                            standard deviation 1
%                 'zscore2', zscore normalization scaled to the interval 
%                           [0,max(power)]
%                 'coi',     con of influence
%                 'sig',     countour plot of the 95% of significance
%                 
%        Output:  wsignal,   wavelet data
%                            power matrix, 
%                            powernorm, zscore normalized matrix
%                            other data, see wsignal output in workspace
%                 
%        Example: [wHPC]=wavelet_spectrogram (HPC,[3 6]);
%        Example: [wHPC]=wavelet_spectrogram (HPC,[3 6],'nonplot','zscore');
%        Example: [wHPC]=wavelet_spectrogram (HPC,[3 6],'zscore','coi','sig');
% Vicent Teruel, Neuronal Circuits Lab, september 2012
% Modified in may 2013
% ------------------------------------------------------------------------------------------
 
warning('off');  
grafico=true;
%scale_div=1;
z_score=false;
z_score2=false;
z_score3=false;
coi=false;
sig=false;
average=false;
scale_normalized=true;
 
 
    if band(1)==0 || band(2) == 0
        fprintf('\n\tFrequencies choosed should be different to 0 \n');
        fprintf('\t Frequency changed to 0.5 \n');
        if band(1)==0 band(1)=0.5;, end;
        if band(2)==0 band(2)=0.5;, end;
        if band(2)<band(1)
            temp =band(2);
            band(2)=band(1);
            band(1)=temp;
        end
         
    end;
 
%----------input variables -----------
Args=struct('nonplot',0, ...
       'zscore',0,...
       'zscore2',0,...
       'zscore3',0,...
       'coi',0,...
       'sig',0,...
       'mean',0,...
       'nonscaled',0); 
 
Args=parseArgs(varargin,Args, ...
          {'nonplot';'zscore';'zscore2';'zscore3';'coi';'sig';'mean';...
          'nonscaled'});
if Args.nonplot==1 grafico=false; end
if Args.zscore==1 z_score=true; end
if Args.zscore2==1 z_score2=true; end
if Args.zscore3==1 z_score3=true; end
if Args.zscore3==1 z_score3=true; end
if Args.nonscaled==1 scale_normalized=false; end
 
if Args.coi==1 coi=true; end
if Args.sig==1 sig=true; end
%---------------------------------------
   colormap jet
   dt=params.dt;
   n=size(signal.data,1);
    
   sigma2=var(signal.data);
   %signal.data = (signal.data - mean(signal.data))/sqrt(sigma2) ;
 
   wavelet_signal.signal=signal.data;
   wavelet_signal.fs = signal.fs;
   wavelet_signal.params=params;
   wavelet_signal.variance = sigma2;
    
   [wavelet_signal.wave,wavelet_signal.period,wavelet_signal.scale,wavelet_signal.coi] = ...
       wavelet(signal.data,params.dt,params.pad,params.dj,params.s0,params.j1,params.mother);
    
   wavelet_signal.power     = abs(wavelet_signal.wave).^2 ;
   wavelet_signal.amplitud  = abs(wavelet_signal.wave);
   wavelet_signal.phase     = angle(wavelet_signal.wave);
    
   wavelet_signal.t = signal.t;
       
   yscale=log2(wavelet_signal.period);
   Yticks = 2.^(fix(log2(min(wavelet_signal.period))):fix(log2(max(wavelet_signal.period))));
     
   power2 = zeros(size(wavelet_signal.power));
   scales = (wavelet_signal.scale')*(floor(ones(1,n)));  % expand scale --> (J+1)x(N) array
   if scale_normalized
       %power2 = wavelet_signal.power./scales;
       power2 = bsxfun(@rdivide, wavelet_signal.power,scales);
   else
       power2 = wavelet_signal.power;
   end
   wavelet_signal.power = power2;%
         
   if z_score || z_score2
       power2 = zscore(power2);
       wavelet_signal.powernorm = power2;
   end
    
   if z_score2
      power2=(min(min(power2))-power2)./min(min(power2)); 
   end
    
   if z_score3
      for i=1:size(power2,1)
          power2(i,:)=zscore(power2(i,:));
      end
   end
    
   if grafico 
       h=imagesc(signal.t,yscale,power2);
       wavelet_signal.h = h;
   end
    
 
if coi
hold on;
        tt=[signal.t([1 1])-dt*.5;signal.t;signal.t([end end])+dt*.5];
        hcoi=fill(tt,log2([wavelet_signal.period([end 1]) wavelet_signal.coi wavelet_signal.period([1 end])]),'w');
        set(hcoi,'alphadatamapping','direct','facealpha',.5)
hold off;
end
 
if sig
        signif = wave_signif(1,dt,wavelet_signal.scale,-1,-1,0.95);   
        sig95 = (signif')*(ones(1,n));  % expand signif --> (J+1)x(N) array
        sig95 = bsxfun(@rdivide, wavelet_signal.power,sigma2*sig95);
        %sig95 = wavelet_signal.power ./ (sigma2*sig95);
        wavelet_signal.sig95=sig95;
        if grafico==1
        hold on;
            [c,h] = contour(signal.t,yscale,sig95,[1 1],'k'); %#ok
            set(h,'linewidth',1)
        hold off
        end
end
 
if grafico
set(gca,'YLim',([-log2(band(2)) -log2(band(1))]), ...
    'YDir','reverse', ...
    'YTick',log2(Yticks(:)), ...
    'YTickLabel',num2str(1./Yticks'), ...
    'layer','top');
 
colorbar;
     
    title('Wavelet spectrogram');
    ylabel('Frequency (Hz)');
    xlabel('t (s)');
end   
 
warning('on');
 
     
end