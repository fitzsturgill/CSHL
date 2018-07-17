%%
figure(2);  clf;  hold on;
pNum = 10
for s = 1:pNum,
    
    tic
    sPerm = perms(1:min(10,s));
    toc
    
    tic;
    for i = 1:1000,
        findPermutation(1:pNum,mod(i,s)+1,sPerm);
    end;
    t = toc
    plot(s,t,'o')
    
end

%%
figure(1); clf;
pNum = 4
tic
permsResult = perms(1:pNum);
toc
findPermutationResult = zeros(size(permsResult));
tic
for i = 1:factorial(pNum),
    findPermutationResult(i,:) = findPermutation(1:pNum,i);
    %findPermutation(1:pNum,i)
    %pause
end
toc

subplot(1,3,1); imagesc(permsResult);
subplot(1,3,2); imagesc(findPermutationResult);
subplot(1,3,3); imagesc(permsResult- findPermutationResult);