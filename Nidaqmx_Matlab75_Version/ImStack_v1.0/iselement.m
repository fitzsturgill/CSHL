function out = iselement(inputArray, value)
% 1 if value is contained in the array inputArray, 0 otherwise.

c = (inputArray == value);
out = any(c);
