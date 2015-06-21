# Plot radial component of Laplacian-of-Gaussian


Sigma = 0.5


LaplacianOfGaussian(R, Sigma) = \
    ((R**2 / (2.0 * Sigma**2)) - 1.0) * \
    exp(-(R**2 / (2.0 * Sigma**2))) / \
    (pi * Sigma**4)


set key bottom right
set grid
set samples 1023
set xlabel "R"
#set xtics 0.2
set ylabel "k(R)"
#set ytics 0.2


load DATA_DIR . "/colors.gp"


plot [R = 0:2] LaplacianOfGaussian(R, Sigma)
