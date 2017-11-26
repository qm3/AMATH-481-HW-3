function dw = rhs_ft(t, w0_ft_col, dummy, n, KXY_col,KXY2D, KX, KY, v)

w0_ft = reshape(w0_ft_col, n, n);

psi = -w0_ft./KXY2D;

psix = ifft2(i*KX.*psi);
psiy = ifft2(i*KY.*psi);

wx = ifft2(i*KX.*w0_ft);
wy = ifft2(i*KY.*w0_ft);

poisson = reshape(psix .* wy - psiy .* wx, n^2, 1);

lapl_w = reshape(ifft2(v*(- KXY2D .* w0_ft)), n^2, 1);

dw = lapl_w - poisson;

end

