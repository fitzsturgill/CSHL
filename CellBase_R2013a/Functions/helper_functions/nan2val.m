function Y = nan2val(X, val)
%NAN2ZERO   Replace NaNs with zeros.
%   Y = NAN2ZERO(X) replaces NaNs in X with zeros.
%
%   See also ZERO2NAN.

Y = X;
Y(isnan(X)) = val;