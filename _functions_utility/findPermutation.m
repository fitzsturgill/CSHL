function B = findPermutation(A,x,s)

% Lazily finds permutations of the elements in A up to arbitrary size.
% s is a short-circuit to use the results of perms once the permutation array is small enough
% - if s is an integer, than findPermutation will call perms(A) at the point when numel(A) <= s during recursion
% - if s is an array of size = n! x n, then that array is used as if it  were the output of perms(A) at the same point
% Implemented from here: http://stackoverflow.com/questions/352203/generating-permutations-lazily
% Note that this is MUCH slower than 

[~,maxsize] = computer;
N = numel(A);
%x = (factorial(N)-x+1;
%x = x-1;

if nargin < 3,
    s = 0; %This is efficient on my computer.
end

% Short-circuit because perms is much faster
if numel(s) == 1
    if N <= s,
        %fprintf('Short circuiting at N = %.0f\n',N)
        p = perms(A);
        B = p(x,:);
        return
    end
else
    assert( size(s,1) == factorial(size(s,2)), 'Error: s is poorly shaped.  Should be shaped as result of perms (ie., n! x n)')
    B = s(x,:);
    return
end



if N == 1,
    B = A;
else
    ind = floor(mod((x-1)/factorial(N-1),N))+1;               % set first member of output array  to A (( x / (N-1)! ) mod N)
    ind = N - floor(mod((x-1)/factorial(N-1),N));               % set first member of output array  to A (( x / (N-1)! ) mod N)
    %ind = N - floor(mod((x-1)/factorial(N-1),N));
    %ind = N - floor(mod(x/factorial(N-1),N));   
    
    % set first member of output array  to A (( x / (N-1)! ) mod N)
    B = [ A(ind) ...
          findPermutation( A([1:ind-1 ind+1:end]),...       % pass remainder of input array
                           mod((x-1),factorial(N-1))+1,...        % ask for permutation ( x mod (N-1)! ) of the remaining elements of A
                           s),...                           % pass the short-circuit variable
        ];        
end




% To find permutation x of array A, where A has N elements:
% 0. if A has one element, return it
% 1. set p to ( x / (N-1)! ) mod N
% 2. the desired permutation will be A[p] followed by
%    permutation ( x mod (N-1)! )
%    of the elements remaining in A after position p is removed

% perm 13 of ABCD: {p = (13 / 3!) mod 4 = (13 / 6) mod 4 = 2; ABCD[2] = C}
% C followed by perm 1 of ABD {because 13 mod 3! = 13 mod 6 = 1}
%   perm 1 of ABD: {p = (1 / 2!) mod 3 = (1 / 2) mod 2 = 0; ABD[0] = A}
%   A followed by perm 1 of BD {because 1 mod 2! = 1 mod 2 = 1}
%     perm 1 of BD: {p = (1 / 1!) mod 2 = (1 / 1) mod 2 = 1; BD[1] = D}
%     D followed by perm 0 of B {because 1 mod 1! = 1 mod 1 = 0}
%       B (because there's only one element)
%     DB
%   ADB
% CADB



end


%% Testing code


% pNum = 8
% tic
% permsResult = perms(1:pNum);
% toc
% findPermutationResult = zeros(size(permsResult));
% tic
% for i = 1:factorial(pNum),
%     findPermutationResult(i,:) = findPermutation(1:pNum,i);
% end
% toc
% 
% figure;  clf;  hold on;
% for s = 1:9,
%     tic;
%     for i = 1:1000,
%         findPermutation(1:10,i,s);
%     end;
%     t = toc;
%     plot(s,t,'o')
% end