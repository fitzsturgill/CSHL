function X_rotated = random_rotation(X)
% Performs a random rotation of X (samples x dimensions)
%
%
% Alex Vaughan, 2015

% Must have dims > 2
assert(size(X,2)>2)

% Create random orthogonal matrix
[Q,R] = qr(randn(size(X,2)));
Q = Q*diag(sign(diag(R)));

% If det(Q) == -1, swap first two columns
% to change sign of determinant. We occasionally
% see floating point errors here, so compute to ~20 bits of precision
if round(det(Q)*1E6)/1E6 == -1,
    Q = Q(:,[2 1 3:end]);
end
assert(round(det(Q)*1E6)/1E6 == 1)

% Multiply input matrix X by your orthogonal matrix to get a rotation.
X_rotated = X * Q;
