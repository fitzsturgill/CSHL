function startpoint = getStartPointOfLine(axis, line)

startpoint = startEndMidOfDrawnLine(axis, line);
startpoint = startpoint(1,1:2);
