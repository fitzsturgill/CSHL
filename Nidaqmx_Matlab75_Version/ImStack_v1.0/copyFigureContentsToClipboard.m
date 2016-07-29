function copyFigureContentsToClipboard(figurehandle)
global gh state

eval(['print -dmeta -' figurehandle ' -noui']);