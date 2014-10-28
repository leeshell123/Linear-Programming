function [X,FVAL,EXITFLAG,BAS,NBAS] = mylinprog(f,A,b,Aeq,beq,LB,UB)
% Inputs:
% P: a m x n vector is the problem matrix
% q: a m x 1 vector is the RHS constant. Note that we will guarantee q(i) >= 0 to avoid initialization phase.
% r: a n x 1 vector of objective coefficients
%
% Outputs:
% optVal: the optimal objective value
% optSolution: a n x 1 vector of values for the optimal solution for original problem variables
% 

%% Sort the problem in standard form
[mA,nA] = size(A);
[mLB,nLB] = size(LB);
[mUB,nUB] = size(UB);
[mAeq,nAeq] = size(Aeq);

if(mLB==nA)
    A = [A; -eye(mLB)];
    b = [b; -LB];
end

if(mUB==nA)
    A = [A; eye(mUB)];
    b = [b; UB];
end

if(mAeq==nA)
    A = [A; -Aeq; Aeq]
    b = [b; zeros(nA,1); zeros(nA,1)];
end

[m,n] = size(A);

%% Filter out malformed inputs


%% Create augmented matrices with slack variables added.
A  =[ A eye(m)];
c = [ -f; zeros(m,1)]; % inverse the 'f' before calling pivot.m

%% start off the basis and non basic indices
bas1 = (n+1:m+n)'; %% basis is {n+1..n+m}
nonbas1 = (1:n)'; %% nonbasic are {1..n}
isFinal = 0; %% Flags for isFinal
isUnbounded = 0; %% isUnbounded

%% Initialization phase

%% Repeatedly pivot to get the solution
while (isFinal == 0 && isUnbounded == 0)
    % Call pivot function to update the dictionary
    [bas2,nonbas2,objVal,isFinal,isUnbounded] = pivot(A,b,c,bas1,nonbas1);
    % Prepare for one more iteration and test.
    bas1 = bas2;
    nonbas1 = nonbas2;
end

%% Compute the optimal solution and set the exit flag
if (isFinal)
    temp = zeros(m+n, 1);
    temp(bas2) = A(:, bas2) \ b;
    X = temp([1:n],1);
    FVAL = f' * X;
    EXITFLAG = 1;
    BAS = bas2;
    NBAS = nonbas2;
else
    X = [];
    FVAL = inf;
    EXITFLAG = -3;
    BAS = bas2;
    NBAS = nonbas2;
end;

end
