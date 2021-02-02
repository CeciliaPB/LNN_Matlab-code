function wfData = wf_DataExtraction_spike2

    % Find the mat files, check if they have calculated wf
    [name, ~] = uigetfile('*.mat*', 'Select a Mat File', 'MultiSelect', 'off');
        if (name == 0)
          return;
        end
    
    load(name,'wf');
    if exist('wf','var') == 0
        wf = waveform_analysis_spike2('features','plot');
    end
    
    % load the variables inside wf
    names = fieldnames(wf);
    for ii=1:length(names)
        eval([names{ii} '=wf.' names{ii} ]);
    end
    
    if pk1(1,2)>pk2(1,2) %#ok<NODEF> % In case there are no two valleys
        b     = a;      %#ok<NODEF>
        a     = 0;
        pk3   = pk1;    %#ok<NASGU>
        pk1   = [0,0];  %#ok<NASGU>
        d     = abs(c); %#ok<NODEF>
        c     = NaN;
        f     = e;      %#ok<NODEF>
        e     = NaN;
        first = NaN;
    elseif isempty(nonzeros(pk1)) == 1
        c     = NaN;
        e     = NaN;
        first = NaN;
    elseif isempty(nonzeros(pk3)) == 1 %#ok<NODEF> % In case there are no two valleys
        b     = a;      %#ok<NODEF>
        a     = 0;
        pk3   = pk2;    %#ok<NASGU>
        pk2   = pk1;    %#ok<NASGU>
        pk1   = [0,0];  %#ok<NASGU>
        d     = c;      %#ok<NODEF>
        c     = NaN;
        f     = e;      %#ok<NODEF>
        e     = NaN;
        first = NaN;
    end    
    
    % Calculations
    sym  = abs(a)/(abs(a)+abs(b));
    amp1 = abs(e); % pk2-pk1
    amp2 = abs(f); % pk3-pk2
    if abs(c) > 1
        pkdist1 = (15 * 0.5 )/c; % in ms, distance pk1 to pk2
        pkdist2 = (15 * 0.5 )/d; % in ms, distance pk2 to pk3
        halfwidth = (15 * 0.5 )/g; % in ms, halfwith pk2
        pkwidth1  = abs((15 * 0.5 )/first); % in ms, width pk1
        pkwidth2  = abs((15 * 0.5 )/last); % in ms, width pk3
    else
        pkdist1 = c*1000; % in ms, distance pk1 to pk2
        pkdist2 = d*1000; % in ms, distance pk2 to pk3
        halfwidth = g*1000; % in ms, halfwith pk2
        pkwidth1  = abs(first*1000); % in ms, width pk1
        pkwidth2  = abs(last*1000); % in ms, width pk3
    end
    
    wfData = table;
    wfData.sym       = sym;
    wfData.amp1      = amp1;
    wfData.amp2      = amp2;
    wfData.pkdist1   = pkdist1;
    wfData.pkdist2   = pkdist2;
    wfData.halfwidth = halfwidth;
    wfData.pkwidth1  = pkwidth1;
    wfData.pkwidth2  = pkwidth2;
    
    disp(['Save as ' name '?']);
    saveID = input('[Y/ N]:','s');
    switch saveID
        case {'Y', 'y',''}
            save(name(1:end-4), 'wfData', '-append');
        case {'N', 'n'}
            name2 = input('Choose file name:','s');
            save(name2, 'wfData', '-append');
    end
    
    clearvars -except wfData name
end
