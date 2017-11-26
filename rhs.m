function dw = rhs(t, w0, dummy, N, v, A, L, U, B, C)

y = L\w0;
psi = U\y;

dw = v * A * w0 + (C * psi) .* (B * w0) - (B * psi) .* (C * w0);
end
