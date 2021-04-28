
function load_open_ephys_binary_noADC(folder, channel_num)

% -------------------------------------------------------------------------
% Dániel Magyar, 2021
% Laboratory of Network Neurophysiology
% Institute of Experimental Medicine, Hungary.
% -------------------------------------------------------------------------

info = load_open_ephys_binary([folder '\structure.oebin'],'continuous',1, 'mmap');

first_channel = info.Data.Data(1).mapped(1,1:end);
incoming_data = zeros(channel_num,length(first_channel),'int16');
incoming_data(1,:) = int16(first_channel);

for ii = 2:channel_num
    [next_channel] = info.Data.Data(1).mapped(ii,1:end);
    incoming_data(ii,:) = int16(next_channel);
end

fid = fopen('continuous_noADC.dat','w+');
fwrite(fid,incoming_data(:,:),'int16');
fclose(fid);
disp('Done!')
end