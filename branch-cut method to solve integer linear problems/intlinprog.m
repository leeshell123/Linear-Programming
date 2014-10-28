function [X,FVAL,EXITFLAG,P,q] = intlinprog(f,A,b,Aeq,beq,LB,UB,ALG,DEP)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%% Filter out the malformed inputs

%% Initialization
P = A;
q = b;

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
    A = [A; -Aeq; Aeq];
    b = [b; zeros(nA,1); zeros(nA,1)];
end

[m,n] = size(A);
A  = round(A*1e+6);
b = round(b*1e+6);

%% Solve the LP relaxation of the given ILP and export the final dictionary
[X,FVAL,EXITFLAG,BAS,NBAS] = mylinprog(f,A,b,Aeq,beq,[],[]);
if (EXITFLAG==0)
    return;
end
if (EXITFLAG==-2)
    return;
end
if (EXITFLAG==-3)
    return;
end
if (EXITFLAG==1)
    % Is it an integer solution, if yes then stop.
    nonIntIndex = [];
    nonIntIndex = find((X - round(X)).^2 > 1e-6, 1); % find a non integer variable
    if (isempty(nonIntIndex))
        return;
    end
end

%% Repeatedly adding cutting plane for fixed many times if possible
% Expand matrix A and fomulate the ILP into standard form
if (ALG=='CP')
    EXITFLAG = -8;
    A1 = A;
    b1 = b;
    for i = 1:DEP
       %% Call the core function _addCuttingPlane and get new dictionary
        [acp,bcp] = addCuttingPlane(A1,b1,BAS,NBAS);
        A2 = [A1; acp];
        b2 = [b1; bcp];
        P = [P; acp/1e+6];
        q = [q; bcp/1e+6];
       %% Solve the dual LP relaxation of the new ILP
        [X,FVAL,FLAG,BAS,NBAS] = mylinprog(f,A2,b2,[],[],[],[]);
        %[X,FVAL,EXITFLAG] = linprog(f,A2,b2,[],[],zeros(16,0),[]);
        %[X,FVAL,EXITFLAG] = linprog(b2,-A2',f,[],[],zeros(7,0),[]);
       %% Analyze the LP relaxation
        if (FLAG==0)
            return;
        end
        if (FLAG==-2)
            EXITFLAG = -2;
            return;
        end
        if (FLAG==-3)
            EXITFLAG = -3;
            return;
        end
        if (FLAG==1)
            % Is it an integer solution, if yes then stop.
            nonIntIndex = [];
            nonIntIndex = find((X - round(X)).^2 > 1e-6, 1); % find a non integer variable
            if (isempty(nonIntIndex))
                EXITFLAG = 1;
                return;
            end
        end
       %% Add an additional cutting plane
        A1 = A2;
        b1 = b2;
    end
end

end

