function syncNIDAQBoards(master, slave)

% This function will tie the board clocks for many boards together.
% The master is an object from the main board
% The slave is an array of objects, one from each of the slave boards.
% This uses the internal RTSI Bus, so they need to be wired together inside the computer.
% Use master as fast board, and slow boards as slaves.



if nargin ~= 2
	error('syncNIDAQBoards: must supply a master and slave objects to sync boards.');
end

InfoMaster=daqhwinfo(master);
MasterID=InfoMaster.ID;
daqmex(master,'call', 'select_signal', 32100, 12170,15900); % output clock from master on RTSI Pin 1

for i = 1:length(slave)
	InfoSlave=daqhwinfo(slave(i));
	SlaveID=InfoSlave.ID;
	if MasterID == SlaveID	% Check to make sure board IDs are different for master and slave
		error('syncNIDAQBoards: Master and slave cannot be from the same board. Check daqhwinfo(obj). ');
	else
		daqmex(slave(i),'call', 'select_signal', 12170, 32100,15900);% slave slow board clock on RTSI Pin 1
		disp(['Tied Clocks for NIDAQ boards # ' num2str(MasterID) ' and # ' num2str(SlaveID) ' together.']); 
	end
end

