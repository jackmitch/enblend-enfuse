# Plot example of exposure-cutoff function


Lower = 0.05
Upper = 0.97
Mu = 0.5
Sigma = 0.2


_Epsilon = 1.0 / 1024.0
Step(X) = X < 0 ? 0 : (X > _Epsilon ? 1 : 1/0)

Gaussian(Y, Mu, Sigma) = exp(-0.5 * ((Y - Mu) / Sigma)**2)

ExposureCutoffProper(Y, LowerCutoff, UpperCutoff) = \
    Y <= LowerCutoff ? 0 : (Y >= UpperCutoff ? 0 : Gaussian(Y, Mu, Sigma))

ExposureCutoff(Y, LowerCutoff, UpperCutoff) = \
    Step(Y - LowerCutoff) * \
    Step(UpperCutoff - Y) * \
    ExposureCutoffProper(Y, LowerCutoff, UpperCutoff)


set grid
set key right bottom
set samples 1023
set xlabel "Y"
set xtics 0.2
set yrange [-0.1:1.1]
set ytics 0.2


load DATA_DIR . "/colors.gp"


plot [Y = 0:1] ExposureCutoff(Y, Lower, Upper)
