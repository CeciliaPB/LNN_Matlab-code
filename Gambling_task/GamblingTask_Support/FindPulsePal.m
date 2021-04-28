function [PulsePort, RawString, Status, PortLocations, nPorts, Clip, x] = FindPulsePal
%FINDPULSEPAL  PulsePal automatic recognition
% This code finds the current COMport of the PulsePal via searching
% it in the raw string which contains all the available instruments (name,
% port, vendorID, productID etc..)
%
% To find ALL COMport connected devices use: 
% system(['wmic path Win32_SerialPort'])
% This will give a list. Check the one where the PulsePal is connected to
% and change the "PNPDeviceID LIKE ''%7&1755%''" to the one actually
% corresponding to the PulsePal: "PNPDeviceID LIKE ''%HEREYOURid%''".
% -------------------------------------------------------------------------
% Modified by Cecília Pardo-Bellver, 2021
% Laboratory of Network Neurophysiology
% Institute of Experimental Medicine, Hungary.
% -------------------------------------------------------------------------

% delete(instrfindall); %deleting former

[Status, RawString] = system('wmic path Win32_SerialPort Where "PNPDeviceID LIKE ''%7&1755%''" Get DeviceID'); %Search for the available ports at the moment
PortLocations = strfind(RawString, 'COM'); %searching the COM word from the former string
PulsePort = cell(1,100);
nPorts = length(PortLocations);
for x = 1:nPorts
    Clip = RawString(PortLocations(x):PortLocations(x)+6);
    PulsePort{x} = Clip(1:find(Clip == 32,1, 'first')-1);
end
PulsePort = PulsePort(1:nPorts);