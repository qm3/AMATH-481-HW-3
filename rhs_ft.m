function dw = rhs_ft(t, w0_col, dummy, n, KXY2D, v, A, B, C)

w0_ft = fft2(reshape(w0_col, n, n));

psi = reshape(real(ifft2(-w0_ft./KXY2D)), n^2, 1);

dw = v * A * w0_col + (C * psi) .* (B * w0_col) - (B * psi) .* (C * w0_col);

end

