
% testscript to communicate with the arduino uno
portname = '/dev/cu.usbmodem1411';
Arduino_Serial = serial(portname, 'BaudRate', 9600, 'DataBits', 8, 'StopBits', 1, 'Timeout', 1, 'DataTerminalReady', 'off');
set(Arduino_Serial, 'OutputBufferSize', 8000);
set(Arduino_Serial, 'InputBufferSize', 50000);
fopen(Arduino_Serial);
        
%fwrite(Arduino_Serial, char(50));

%fclose(instrfind)