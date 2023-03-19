function [x,lambda]=EqualityQPSolver(H,g,A,b,solver)
%EqualityQPSolver  Equality constrained convex QP

%Syntax: [x,lambda]=EqualityQPSolver(H,g,A,b,solver)
%solver is a flag used to switch between the different factorizations
%solver : 'LUdense'----->LUfactorization (dense)
%         'LUsparse'---->LUfactorization (sparse)
%         'LDLdense'---->LDL-factorization (dense)
%         'LDLsparse'--->LDL-factorization (sparse)
%         'RangeSpace'-->Range-Space factorization
%         'NullSpace'--->Null-Space factorization
switch solver
    case 'LUdense'
        [x,lambda]=EqualityQPSolverLUdense(H,g,A,b);
    case 'LUsparse'
        [x,lambda]=EqualityQPSolverLUsparse(H,g,A,b);
    case 'LDLdense'
        [x,lambda]=EqualityQPSolverLDLdense(H,g,A,b);
    case 'LDLsparse'
        [x,lambda]=EqualityQPSolverLDLsparse(H,g,A,b);
    case 'RangeSpace'
        [x,lambda]=EqualityQPSolverRangeSpace(H,g,A,b);
    case 'NullSpace'
        [x,lambda]=EqualityQPSolverNullSpace(H,g,A,b);
    otherwise 
        x=[];
        lambda=[];
end
