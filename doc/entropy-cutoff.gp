# Plot example of entropy-cutoff function


Lower = 0.05
Upper = 0.9
Offset = 5e-3

_Epsilon = 1.0 / 1024.0
Step(X) = X < 0 ? 0 : (X > _Epsilon ? 1 : 1/0)

EntropyCutoffProper(Y, LowerCutoff, UpperCutoff) = \
    Y <= LowerCutoff ? 0.0 : (Y >= UpperCutoff ? 1.0 : Y)

EntropyCutoff(Y, LowerCutoff, UpperCutoff) = \
    Step(Y - LowerCutoff) * Step(UpperCutoff - Y) * \
    EntropyCutoffProper(Y, LowerCutoff, UpperCutoff) + \
    Step(Y - UpperCutoff)


set grid
set key right bottom
set samples 1023
set xlabel "Y"
set xtics 0.2
set ylabel "EC"
set yrange [-0.1:1.1]
set ytics 0.2


load DATA_DIR . "/colors.gp"


plot [Y = 0:1] Y + Offset, EntropyCutoff(Y, Lower, Upper) - Offset
