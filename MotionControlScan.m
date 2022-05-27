try % In case the oscilloscope is being used, stop it
    % Stop the device
    [status.stop] = invoke(ps5000aDeviceObj, 'ps5000aStop');
    % Disconnect device
    % Disconnect device object from hardware.
    disconnect(ps5000aDeviceObj);
    delete(ps5000aDeviceObj);
end

clc
clear all
close all

%% user-defined parameters
filename          ='carbonwishboneflatsetup';  %descriptive name, 4 random digits will be automatically appended to this
pointsBefore      = 64;  %points recorded before trigger
pointsAfter       = 4096;  %points recorded after trigger
triggerVoltage    = -3500;  %trigger voltage (mV)
startY            = -30;  %botton of scan (mm)
endY              = 30;  %top of scan (mm)
startX            = -40;  %left of scan (mm)
endX              = 40;  %right of scan (mm)
scanResolution    = 1;  %scan resoluiton (mm)
startData         = 500;  %begining of waveform data to be saved (also affects preview)
endData           = 2500;  %end of waveform data to be saved (also affects preview)
minPeakProminence = 50;  %only affects preview visualization, can be changed in post-process
rmsWindow         = 30;  %tune to acheive one envelope arc per pulse

%% connect to scanner and configure

s = serialport("COM3",115200);%COM 5 and 250000 for connecting to lulzbot
pause(5)% allow time for the PC to connect to the stages
%% home scanner and move to start position
reply="";%initialize
while (s.NumBytesAvailable>0)
    reply=readline(s);%returns the recieved text
end
%writeline(s,"$H"); %uncomment this line to home the printer. This is usually uneccesary
userinput = input('is the scanner at the correct start position?(yes/no)','s')
if ~strcmp(userinput,'yes')
    error('answer was not yes (case sensitive)')
end

writeline(s,"$120=200");% set acceleration to a high value
pause(1) %give a moment for the previous line to take effect
writeline(s,"$121=200");% set acceleration to a high value
writeline(s,"G92 X0 Y0 Z0"); % resets origin coordinates to current position
%% PicoScope 5000 Series (A API) Instrument Driver Oscilloscope Rapid Block Data Capture Example
%%%%%%%%%%%%%%%%%%%%%%%%
% The following code is based on the example code from PicoScopeSDK
%%%%%%%%%%%%%%%%%%%%%%%%
% This is an example of an instrument control session using a device
% object. The instrument control session comprises all the steps you are
% likely to take when communicating with your instrument.
%
% These steps are:
%
% # Create a device object
% # Connect to the instrument
% # Configure properties
% # Invoke functions
% # Disconnect from the instrument
%
% To run the instrument control session, type the name of the file,
% PS5000A_ID_Rapid_Block_Example, at the MATLAB command prompt.
%
% The file, PS5000A_ID_RAPID_BLOCK_EXAMPLE.M must be on your MATLAB PATH.
% For additional information on setting your MATLAB PATH, type 'help
% addpath' at the MATLAB command prompt.
%
% *Example:*
%    PS5000A_ID_Rapid_Block_Example;
%
% *Description:*
%     Demonstrates how to call functions in order to capture a series of
%     waveforms using rapid block mode on a PicoScope 5000 Series
%     Oscilloscope using the underlying 'A' API library functions.
%
% *See also:* <matlab:doc('icdevice') |icdevice|> | <matlab:doc('instrument/invoke') |invoke|>
%
% *Copyright:* © 2013-2018 Pico Technology Ltd. See LICENSE file for terms.

%% Suggested input test signals
% This example was published using the following test signal:
%
% * Channel A: 4 Vpp Swept sine wave (Start: 10 kHz, Stop: 100 kHz, Sweep type: Up, Increment Time: 1 ms, Increment type: Linear, Mode: Continous)
% * Channel B: 2 Vpp Swept square wave (Start: 10 kHz, Stop: 50 kHz, Sweep type: Up, Increment: 5 kHz, Increment Time: 1 ms)

%% Load configuration information

PS5000aConfig;

%% Device connection

% Check if an Instrument session using the device object |ps5000aDeviceObj|
% is still open, and if so, disconnect if the User chooses 'Yes' when prompted.
if (exist('ps5000aDeviceObj', 'var') && ps5000aDeviceObj.isvalid && strcmp(ps5000aDeviceObj.status, 'open'))

    openDevice = questionDialog(['Device object ps5000aDeviceObj has an open connection. ' ...
        'Do you wish to close the connection and continue?'], ...
        'Device Object Connection Open');

    if (openDevice == PicoConstants.TRUE)

        % Close connection to device.
        disconnect(ps5000aDeviceObj);
        delete(ps5000aDeviceObj);

    else

        % Exit script if User selects 'No'.
        return;

    end

end

% Create a device object.
ps5000aDeviceObj = icdevice('picotech_ps5000a_generic', '');

% Connect device object to hardware.
connect(ps5000aDeviceObj);

%% Set channels
% Default driver settings applied to channels are listed below - use the
% Instrument Driver's |ps5000aSetChannel()| function to turn channels on or
% off and set voltage ranges, coupling, as well as analog offset.

% In this example, data is collected on channels A and B. If it is a
% 4-channel model, channels C and D will be switched off if the power
% supply is connected.

% Channels       : 0 - 1 (ps5000aEnuminfo.enPS5000AChannel.PS5000A_CHANNEL_A & PS5000A_CHANNEL_B)
% Enabled        : 1 (PicoConstants.TRUE)
% Type           : 1 (ps5000aEnuminfo.enPS5000ACoupling.PS5000A_DC)
% Range          : 8 (ps5000aEnuminfo.enPS5000ARange.PS5000A_5V)
% Analog Offset  : 0.0 V

% Channels       : 2 - 3 (ps5000aEnuminfo.enPS5000AChannel.PS5000A_CHANNEL_C & PS5000A_CHANNEL_D)
% Enabled        : 0 (PicoConstants.FALSE)
% Type           : 1 (ps5000aEnuminfo.enPS5000ACoupling.PS5000A_DC)
% Range          : 8 (ps5000aEnuminfo.enPS5000ARange.PS5000A_5V)
% Analog Offset  : 0.0 V

% Find current power source
[status.currentPowerSource] = invoke(ps5000aDeviceObj, 'ps5000aCurrentPowerSource');

if (ps5000aDeviceObj.channelCount == PicoConstants.QUAD_SCOPE && status.currentPowerSource == PicoStatus.PICO_POWER_SUPPLY_CONNECTED)

    [status.setChB] = invoke(ps5000aDeviceObj, 'ps5000aSetChannel', 1, 0, 1, 8, 0.0);
    [status.setChC] = invoke(ps5000aDeviceObj, 'ps5000aSetChannel', 2, 0, 1, 8, 0.0);
    [status.setChD] = invoke(ps5000aDeviceObj, 'ps5000aSetChannel', 3, 0, 1, 8, 0.0);

end

%% Set device resolution

% resolution : 12bits

[status.setDeviceResolution, resolution] = invoke(ps5000aDeviceObj, 'ps5000aSetDeviceResolution', 15);
%% Verify timebase index and maximum number of samples
% Use the |ps5000aGetTimebase2()| function to query the driver as to the
% suitability of using a particular timebase index and the maximum number
% of samples available in the segment selected, then set the |timebase|
% property if required.
%
% To use the fastest sampling interval possible, enable one analog
% channel and turn off all other channels.
%
% Use a while loop to query the function until the status indicates that a
% valid timebase index has been selected. In this example, the timebase
% index of 4 is valid.

% Initial call to ps5000aGetTimebase2() with parameters:
%
% timebase      : 4
% segment index : 0

status.getTimebase2 = PicoStatus.PICO_INVALID_TIMEBASE;
timebaseIndex = 4;

while (status.getTimebase2 == PicoStatus.PICO_INVALID_TIMEBASE)

    [status.getTimebase2, timeIntervalNanoseconds, maxSamples] = invoke(ps5000aDeviceObj, ...
        'ps5000aGetTimebase2', timebaseIndex, 0);

    if (status.getTimebase2 == PicoStatus.PICO_OK)

        break;

    else

        timebaseIndex = timebaseIndex + 1;

    end

end

fprintf('Timebase index: %d, sampling interval: %d ns\n', timebaseIndex, timeIntervalNanoseconds);

% Configure the device object's |timebase| property value.
set(ps5000aDeviceObj, 'timebase', timebaseIndex);

%% Set simple trigger
% Set a trigger on channel A, with an auto timeout - the default value for
% delay is used. The device will wait for a rising edge through
% the specified threshold unless the timeout occurs first.

% Trigger properties and functions are located in the Instrument
% Driver's Trigger group.

triggerGroupObj = get(ps5000aDeviceObj, 'Trigger');
triggerGroupObj = triggerGroupObj(1);

% Set the |autoTriggerMs| property in order to automatically trigger the
% oscilloscope after 1 second if a trigger event has not occurred. Set to 0
% to wait indefinitely for a trigger event.

%set(triggerGroupObj, 'autoTriggerMs', 1);

% Channel     : 0 (ps5000aEnuminfo.enPS5000AChannel.PS5000A_CHANNEL_A)
% Threshold   : 500 mV
% Direction   : 2 (ps5000aEnuminfo.enPS5000AThresholdDirection.PS5000A_RISING)
% Delay       : 0
% Auto trigger: 0 (wait indefinitely)

[status.setSimpleTrigger] = invoke(triggerGroupObj, 'setSimpleTrigger', 0, -3500, 2, 0, 1);

%% Set block parameters and capture data
% Capture a block of data and retrieve data values for channels A and B.

% Block data acquisition properties and functions are located in the
% Instrument Driver's Block group.

blockGroupObj = get(ps5000aDeviceObj, 'Block');
blockGroupObj = blockGroupObj(1);

% Set pre-trigger and post-trigger samples as required - the total of this
% should not exceed the value of |maxSamples| returned from the call to
% |ps5000aGetTimebase2()|. The number of pre-trigger samples is set in this
% example but default of 10000 post-trigger samples is used.

% Set pre-trigger samples.
set(ps5000aDeviceObj, 'numPreTriggerSamples', pointsBefore);
% Set post-trigger smaples
set(ps5000aDeviceObj, 'numPostTriggerSamples', pointsAfter);



%% Loop motion and scan
% run commands
ready=false;
GcodeCheck=false;
%reply="";%initialize
q=1;
k=1;%initialize count
figure1=figure(1);
for x=startX:scanResolution:endX
    writeline(s,append("G0 X",num2str(x)," Y",num2str(startY-2)));% addresses backlash
    for y=startY:scanResolution:endY
        stay=true;
        while (stay)
            while (s.NumBytesAvailable>0)
                reply=readline(s); %returns the recieved text
                if contains(reply,'error') %means that the scanner is in motion
                    ready=false;
                    writeline(s,'$$') % this Gcode returns and error while the scanner is in motion
                    GcodeCheck=true;
                elseif contains(reply,'ok') && GcodeCheck %means that the scanner is idle
                    ready=true;
                    GcodeCheck=false;
                elseif contains(reply,'ok') && ~GcodeCheck %means that the scanner began a motion
                    ready=false;
                    writeline(s,'$$')
                    GcodeCheck=true;
                else
                    ready=false;
                end
            end
            %send command to scanner
            if (ready)
                %take a sample
                % This example uses the |runBlock()| function in order to collect a block of
                % data - if other code needs to be executed while waiting for the device to
                % indicate that it is ready, use the |ps5000aRunBlock()| function and poll
                % the |ps5000aIsReady()| function.

                % Capture a block of data:
                %
                % segment index: 0 (The buffer memory is not segmented in this example)

                [status.runBlock] = invoke(blockGroupObj, 'runBlock', 0);

                writeline(s,append("G0 X",num2str(x)," Y",num2str(y)));
                %writeline(s,'G4 1000')
                %each command either changes x or y
                %if y is odd, x increases and vice versa
                %if x is outside of bounds, change y and reverse direction
                ready=false;
                

                % Retrieve data values:

                startIndex              = 0;
                segmentIndex            = 0;
                downsamplingRatio       = 1;
                downsamplingRatioMode   = ps5000aEnuminfo.enPS5000ARatioMode.PS5000A_RATIO_MODE_NONE;

                % Provide additional output arguments for other channels e.g. chC for
                % channel C if using a 4-channel PicoScope.
                [numSamples, overflow, chA] = invoke(blockGroupObj, 'getBlockData', startIndex, segmentIndex, ...
                    downsamplingRatio, downsamplingRatioMode);

                %save data to variables
                data=chA(startData:min(endData,length(chA)));
                for j=1:length(data)
                    %data(:,k) = chA;
                    X(k)=x;
                    Y(k)=y;
                    Z(k)=j;
                    C(k)=data(j);
                    k=k+1;
                end
                subplot(2,1,1)
                findpeaks(envelope(data,rmsWindow,'rms'),'MinPeakProminence',minPeakProminence);
                [pks,locs]=findpeaks(envelope(data,rmsWindow,'rms'),'MinPeakProminence',minPeakProminence);
                for l=1:length(pks)
                    Xpeaks(q)=x;
                    Ypeaks(q)=y;
                    Cpeaks(q)=pks(l);
                    Zpeaks(q)=locs(l);
                    q=q+1;
                end
                title('A-scan Peaks')
                %plot(chA)
                drawnow
                stay=false;
            end
        end
    end
    subplot(2,1,2)
    %colormap(gray)
    scatter3(Xpeaks,Ypeaks,-Zpeaks,[],Cpeaks,'.')%this gets slower over time
    title('Preview')
    xlabel('mm')
    ylabel('mm')
    zlabel(append(num2str(timeIntervalNanoseconds),' Nanoseconds'))
    drawnow
end


%% Process data
% In this example the data values returned from the device are displayed in
% plots in a Figure.

figure2 = figure('Name','PicoScope 5000 Series (A API) Example - Block Mode Capture', ...
    'NumberTitle','off');

% Calculate time (nanoseconds) and convert to milliseconds.
% Use |timeIntervalNanoseconds| output from the |ps5000aGetTimebase2()|
% function or calculate it using the main Programmer's Guide.
% Take into account the downsampling ratio used.

timeNs = double(timeIntervalNanoseconds) * downsamplingRatio * double(0:numSamples - 1);
timeMs = timeNs / 1e6;

% Channel A
axisHandleChA = subplot(1,1,1);
plot(timeMs, chA, 'b');
title(axisHandleChA, 'Channel A');
xlabel(axisHandleChA, 'Time (ms)');
ylabel(axisHandleChA, 'Voltage (mV)');
grid(axisHandleChA, 'on');

%% Reset position

writeline(s,"G0 X0 Y0")

%% Save data
save(append(filename,num2str(round(1000+8999*rand)),'.mat'))

%% Stop the device

[status.stop] = invoke(ps5000aDeviceObj, 'ps5000aStop');

%% Disconnect device
% Disconnect device object from hardware.

disconnect(ps5000aDeviceObj);
delete(ps5000aDeviceObj);