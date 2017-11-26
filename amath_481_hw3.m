clc; clear all; close all;


nx = 4; ny = 8;
deltax = 2.5; deltay = 1.5;
%dx operator is B
N = nx*ny;
adds = (1/(2*deltax)) *ones(N, 1); subtr = -1 * adds;
B = full(spdiags([adds adds subtr subtr], [8 -24 24 -8], N, N));

adds = (1/(2*deltay)) *ones(N, 1); subtr = -1 * adds;
A = eye(4);
C_setup = full(spdiags([adds subtr], [1 -1], 8, 8));
C_setup(1, ny) = -1/(2*deltay); C_setup(end, end - (ny-1)) = 1/(2*deltay);

C = kron(A, C_setup);

%% HW 3 part 2
clc; clear all; close all;

Nx = 64; Ny = 64;
N = Nx * Ny;
tspan = 0:0.5:20;
xspan = linspace(-5, 5, Nx + 1); 
deltax = xspan(2) - xspan(1);
yspan = linspace(0, 12, Ny + 1);
deltay = yspan(2) - yspan(1);
v = 1e-3;

adds = ones(Ny, 1)* 1/(2*deltay); subs = -1 * adds;

C_setup = full(spdiags([adds subs], [1 -1], Ny, Ny));
C_setup(1, Ny) = -1/(2*deltay); C_setup(end, end - (Ny-1)) = 1/(2*deltay);
I = eye(Ny);
C = sparse(kron(I, C_setup));



adds = ones(N, 1)* 1/(2*deltax); subs = -1 * adds;
B = spdiags([adds adds subs subs], [Nx (-N + Nx) (N - Nx) -Nx], N, N);
A = two_d_lap_gen(10, 12, Nx, Ny);

x = xspan(1:Nx); y = yspan(1:Ny);
[X, Y] = meshgrid(x, y);
w_norm = .01;

w0 = exp((-X.^2 - (Y - 6).^2)./25) - w_norm;

A(1,:) = 0; A(1,1) = 1;

%tol = 1e-6 ; options = odeset('RelTol',tol,'AbsTol',tol);
[t, wsol] = ode45('rhs', tspan, w0, [], N, v, A, B, C);







