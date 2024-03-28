N = 30;
M = 90;
L = 100000;
A = rand_ugraph(N, M, 0.2, 0.1);
% L = diag(sum(A)) - A;
H = 0.2*A^3 + 0.5*A^2 + A + 0.1*eye(N);
stmls = randn(N, L);
obsv = H*stmls;
Aest = specTemp(obsv, 1000);
close all;
imagesc(A); title('Ground Truth');
figure;
imagesc(Aest); title('Estimated');