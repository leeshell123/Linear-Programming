function [acp,bcp] = addCuttingPlane(A1,b1,BAS,NBAS,nonIntIndex);
%ADDCUTTINGPLANE Summary of this function goes here
%   Detailed explanation goes here

%% Get the size and select the submatrices for Ab, Ai
[m,n] = size(A1);
A  = [A1 eye(m)];
b = b1;
%A  = round([A1 eye(m)]*1e+6);
%b = round(b1*1e+6);
Ab = A(:,BAS);
Ai = A(:,NBAS);

%% Compute the cutting plane in terms of nonbasic variables
bhat = Ab \ b; % reconstruct the dictionary
Ahat = -Ab \ Ai;

index1 = find((bhat - round(bhat)).^2 > 1e-6, m); % find the index of the first non integer variable
index2 = find(sum(((Ahat-round(Ahat)).^2)')> 1e-6, m);
index = intersect(index1, index2');
if (isempty(index))
    acp = [];
    bcp = [];
    return;
end
temp = index;
index = [];
index = temp(1);

acp = Ahat(index,:);
acp = (acp + floor(-acp)) .* ((acp - floor(acp)).^2>1e-6);

bcp = bhat(index,1);
bcp = (-bcp + floor(bcp)) .* ((bcp - floor(bcp)).^2>1e-6);

%% Express the cutting plabe in terms of decision variables
temp = zeros(1,m+n);
temp(NBAS) = acp;
for i=1:m
    temp([1:n]) = temp([1:n]) - temp(i+n) * A1(i,[1:n]);
    bcp = bcp - temp(i+n) * b1(i);
end
acp = [];
acp = temp([1:n]);

end

