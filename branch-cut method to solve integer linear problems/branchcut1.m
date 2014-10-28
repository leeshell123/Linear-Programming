function x_best = branchcut1(c,A,b,lb,ub)
% branch-cut method 1 add 5 cutting planes to the problem, if still no
% integer solution, then use branch-and-bound to solve the problem with
% cuttingplanes.

if nargin<5
    error('BCPROG: Too few input arguments')
end


%delete all-zero rows
b = b(any(A,2),:);
A = A(any(A,2),:);
[X,FVAL,EXITFLAG,P,q] = intlinprog(c,A,b,[],[],lb,ub,'CP',5);

if(EXITFLAG==1)
    %disp('Integer solution founded by cutting plane method');
    %disp('fval=%f',f);
    x_best=X;        
elseif(EXITFLAG==-8)   
    x_best=bbbprog(c,P,q,lb,ub);    
else
    disp('LP relexation unbounded or infeasible');
end
end


