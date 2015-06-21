# Show the "Log" transform we use to convert floating-point data
# before copying it into a pyramid.


LogTransform(x) = x >= 0.0 ? 1.0 + log(1.0 + x) : 1.0 / (1.0 - x);


set grid
set samples 1023
set xlabel "L"
set ylabel "Log(L)"

unset key


load DATA_DIR . "/colors.gp"


plot [-20 : 100] LogTransform(x)
