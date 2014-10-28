clear;
clc;
load test2.mat
lb=zeros(5,1);
ub=20* ones(5,1);
x_best1 = bbbprog(c,A,b,lb,ub)
x_best2 = branchcut1(c,A,b,lb,ub)
x_best3 = branchcut2(c,A,b,lb,ub)
x_best4 = branchcut(c,A,b,lb,ub)