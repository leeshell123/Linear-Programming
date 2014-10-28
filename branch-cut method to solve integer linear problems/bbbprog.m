function x_best = bbbprog(c,A,b,lb,ub)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


stopTol = 1e-8; %Stopping tolerance of the gap (f_integer-f_lp)/(f_integer+f_lp)
maxNodes = 1e4; %Maximum number of nodes in the branch and bound tree to visit
intTol = 1e-6; %Integer tolerance

if nargin<5
    error('MIPROG: Too few input arguments')
end

[m,n]=size(A)
yidx=true(n,1);

%delete all-zero rows
b = b(any(A,2),:);
A = A(any(A,2),:);

%Assume no initial best integer solution
%Add your own heuristic here to find a good incumbent solution, store it in
%intBest,y_best,x_best
intBest = inf;
y_best = [];
x_best = [];
numIntSol = double(~isempty(y_best));

f = inf(maxNodes,1);
f(1) = 0;
prtOpt = inf;

%The nodeOrder in which the problems shall be solved
nodeOrder = [1];
%The indices of the problems that have been visited
visited = nan(maxNodes,1);

S=lb;
D=ub;
i=0;

%% Branch && bound
while i==0 || isinf(intBest) || (~isempty(nodeOrder) &&  ((intBest-min(prtOpt(nodeOrder)))/abs(intBest+min(prtOpt(nodeOrder))) > stopTol) &&  i<maxNodes)
    %Is the parent node less expensive than the current best
    if i==0 || prtOpt(nodeOrder(1))<intBest
        %Solve the LP-relaxation problem
        i=i+1;
        %lower bound
        s = S(:,nodeOrder(1));
        %upper bound
        d = D(:,nodeOrder(1));
        [x,this_f,flag] = linprog(c,A,b,[],[],s,d);

        %Visit this node
        visited(i) = nodeOrder(1);
        nodeOrder(1) = [];
        if (flag~=1)
            %infeasible, dont branch
            if i==1
                error('MIPROG: Infeasible or unbounded initial LP problem.')
            end
            f(i) = inf;
        elseif flag==1
            %convergence
            f(i) = this_f;  
            if this_f<intBest
                y = x(yidx);
                %fractional part of the integer variables -> diff
                diff = abs(round(y)-y);
                if all(diff<intTol)
                    %all fractions less than the integer tolerance
                    %we have integer solution
                    numIntSol = numIntSol+1;
                    intBest = this_f;
                    y_best = round(x(yidx));
                    x_best = round(x);
                else
                    [~,branch_idx] = max(diff,[],1);
                    
                    %Branch into two subproblems
                    s1 = s;
                    s2 = s;
                    d1 = d;
                    d2 = d;
                    if(s1(branch_idx)<ceil(y(branch_idx)))
                        s1(branch_idx)=ceil(y(branch_idx));
                    end
                    
                    if(d2(branch_idx)>floor(y(branch_idx)))
                        d2(branch_idx)=floor(y(branch_idx));
                    end

                    nsold = size(S,2);
                    
                    % add subproblems to the problem tree
                    S = [S s1 s2];
                    D = [D d1 d2];
                    prtOpt = [prtOpt f(i) f(i)];
                    
                    nsnew = nsold+2;
                    %branch on the best lp solution
                    nodeOrder = [nsold+1:nsnew nodeOrder];
                    [dum,pidx] = sort(prtOpt(nodeOrder));
                    nodeOrder=nodeOrder(pidx);
                     
                end
            end
        else
            error('BCPROG: Problem neither infeasible nor solved')
        end
    else %parent node is more expensive than current f-best -> don't evaluate this node
        nodeOrder(1) = [];
    end
end


disp(['Iteration ', num2str(sum(~isnan(visited))), '. Optimization ended.']);
if isempty(nodeOrder) || intBest<min(prtOpt(nodeOrder))
    disp('Found optimal solution!')
elseif ~isinf(intBest)
    disp('Ended optimization. ');
else
    disp(['Did not find any integer solutions']);
end
disp(['Objective function value: ',num2str(intBest)]);


end

