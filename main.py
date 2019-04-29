from generator import *

xlim = (0.0,100.0)
ylim = (0.0,20.0)
zlim = (0.0,10.0)
n = 10
mean_alpha = 36.8 * pi / 180.0
std_alpha = 14.4 * pi / 180.0
Lf = 13.0


xs, ys, zs = generate_fiber_center(n, xlim, ylim, zlim)

alphas, betas = generate_angles(n, beta_lim = (0.0,2.0*pi), alpha_params = (mean_alpha, std_alpha))

((xi, yi, zi), (xf, yf, zf)) = compute_fiber_edges(xs, ys, zs, alphas, betas, Lf)

print(xs, ys, zs)
print(alphas, betas)

print(xi, yi, zi)
print(xf, yf, zf)

