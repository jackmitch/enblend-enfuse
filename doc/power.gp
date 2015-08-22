default_optimum = 0.5
default_width = 0.2
fwhm_gaussian = 2.3548200450309493820231386529193992755
fwhm_power = 1.6817928305074290860622509524664297901

power(y, y_opt, width) = \
        abs((y - y_opt) / (width * fwhm_gaussian / fwhm_power)) <= 1.0 ? \
        1.0 - (abs(y - y_opt) / (width * fwhm_gaussian / fwhm_power))**4 : \
        0.0


set key bmargin center horizontal
set samples 1023
set xlabel "$Y$"
set xtics 0.2
set ylabel "$w$"
set ytics 0.2


load DATA_DIR . "/colors.gp"


plot [y = 0:1] \
     power(y, default_optimum, 0.1) title "width = 0.1", \
     power(y, default_optimum, 0.2) title "default width = 0.2", \
     power(y, default_optimum, 0.4) title "width = 0.4"
