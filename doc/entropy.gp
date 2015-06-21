# Plot entropy with respect to base 2


_Base = 2.0
lb(x) = log(x) / log(_Base)


Entropy(x) = -x * lb(x) - (1.0 - x) * lb(1.0 - x)


set samples 1023
set xlabel "$p$"
set xtics 0.2
set ylabel "$H_2(p)$"
set ytics 0.2


load DATA_DIR . "/colors.gp"


plot [p = 0:1] Entropy(p)
