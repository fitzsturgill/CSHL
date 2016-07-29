function count = fwriteu8(fname, A)
%FWRITEU8 Append bytes to a file.
%   COUNT = FWRITEU8(FNAME, A) appends the bytes in the uint8
%   array A to the file specified by the string FNAME, returning
%   the count of bytes successfully written.

%   Copyright (c) 1984-98 by The MathWorks, Inc.
%   $Revision: 1.3 $  $Date: 1997/11/21 23:36:00 $

error('Missing MEX-file FREADU8');
