default_optimum = 0.5
default_width = 0.2
fwhm_gaussian = 2.3548200450309493820231386529193992755
fwhm_lorentzian = 2.8284271247461900976033774484193961571

lorentzian(y, y_opt, width) = 1.0 / (1.0 + ((y - y_opt) / (width * fwhm_gaussian / fwhm_lorentzian))**2 / 2.0)


set key bmargin center horizontal
set samples 1023
set xlabel "$Y$"
set xtics 0.2
set ylabel "$w$"
set ytics 0.2


load DATA_DIR . "/colors.gp"


plot [y = 0:1] \
     lorentzian(y, default_optimum, 0.1) title "width = 0.1", \
     lorentzian(y, default_optimum, 0.2) title "default width = 0.2", \
     lorentzian(y, default_optimum, 0.4) title "width = 0.4"
