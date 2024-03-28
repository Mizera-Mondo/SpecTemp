function S = specTemp(X, ep)
%SPECTEMP Given graph spectrum templates, estimate the best matched
%topology under sparest criteria.
% Reference: Network Topology Inference from Spectral Templates

%% IMPORTANT! ADJACENCY WILL CAUSE PROBLEM SINCE NEGATIVE EIGENVALUE. SORT WILL DESTRUCT THE STRUCT.
[N, L] = size(X);
Cx = 1/L*(X*X');
[vec, val] = eig(Cx);

Z = zeros(N);
E = eye(N);
conND = ones(N) - E;
N0 = ones(1, N);
N1 = zeros(N, 1);
N1(1) = 1;
cvx_begin
    variable S(N, N) symmetric
    variable Lbd(N, N) diagonal
    minimize norm(S, 1) + ep*norm(S - vec*Lbd*vec')
    subject to
    diag(S) == 0
    S - diag(diag(S)) >= 0
    N0*S*N1 == 1
cvx_end
end
