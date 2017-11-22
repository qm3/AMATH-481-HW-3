clc; clear all; close all;


nx = 4; ny = 8;
deltax = .25; deltay = 1.5;
%dx operator is B
N = nx*ny;
adds = (1/(2*deltax)) *ones(N, 1); subtr = -1 * adds;
B = full(spdiags([adds adds subtr subtr], [8 -24 24 -8], N, N));

adds = ones(N,1);
for n = 1:4
    adds(8*n) = 0;
end
adds = (1/(2*deltay)) * adds;
subtr = -1 * adds;
adds = flipud(adds);
C = full(spdiags([adds subtr], [1 -1], N, N));
% still need to fill in the odd terms every 8/9th row.