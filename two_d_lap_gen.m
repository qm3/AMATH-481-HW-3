function matrixA = two_d_lap_gen (Lx, Ly, Nx, Ny)

Nx1 =  Nx + 1 ;  Ny1 =  Ny + 1 ;
N = Nx * Ny ;

dx = Lx / Nx ;  idx2 = 1 / dx^2 ;  i2dx = 0.5 / dx ;
dy = Ly / Ny ;  idy2 = 1 / dy^2 ;  i2dy = 0.5 / dy ;

e0 = (-2 * idx2 - 2 * idy2) * ones(N,1) ;
%  \partial_xx
e3  = idx2 * ones(N,1) ;  ie3  = Ny ;
em3 = idx2 * ones(N,1) ;  iem3 = -ie3 ;
e4  = idx2 * ones(N,1) ;  ie4  = (Nx-1) * Ny ;  
em4 = idx2 * ones(N,1) ;  iem4 = -ie4 ;

%  \partial_yy
e1  = idy2 * ones(N,1)  ; ie1  = 1 ;
em1 = idy2 * ones(N,1)  ; iem1 = -ie1 ;
e2  = zeros(N,1) ;        ie2  = Ny - 1 ;
em2 = zeros(N,1) ;        iem2 = -ie2 ;

for ix = 1:Nx
    e1 ((ix-1) *  Ny + 1) = 0 ;
    em1( ix * Ny)         = 0 ;
    e2 ( ix * Ny)         = 1 * idy2 ;
    em2((ix-1) *  Ny + 1) = 1 * idy2 ;
end

vec_ediag  = [em4  em3  em2  em1  e0 e1   e2  e3  e4 ] ;
vec_iediag = [iem4 iem3 iem2 iem1 0  ie1  ie2 ie3 ie4] ;
matrixA = spdiags(vec_ediag, vec_iediag, N, N) ;

end