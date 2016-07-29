function endpoint = getEndPointOfLine(axis, line)

endpoint = startEndMidOfDrawnLine(axis, line);
endpoint = endpoint(2,1:2);
