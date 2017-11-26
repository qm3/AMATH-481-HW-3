function dw = rhs(t, w0, dummy, N, v, A, B, C)

w0_column = reshape(w0, N, 1);

[L, U] = lu(A);
y = L\w0_column;
psi = U\y;

dw = v * A * w0_column + (C * psi) .* (B * w0_column) - (B * psi) .* (C * w0_column);
end
