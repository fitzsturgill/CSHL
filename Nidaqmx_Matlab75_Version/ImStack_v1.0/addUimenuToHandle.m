function addUimenuToHandle(handle, uimenuhandle)
% This function adds the uimenuhandle UI COntextmenu (accessed by
% right clicking on the object) to the handle if possible.
try
    set(handle,'UIContextMenu',uimenuhandle);
catch
    lasterr;
    disp('addUimenuToHandle: Handle does not accept a UIContext Menu.');
end