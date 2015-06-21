default_optimum = 0.5
default_width = 0.2
fwhm_gaussian = 2.3548200450309493820231386529193992755
fwhm_halfsinusodial = 2.0943951023931954923084289221863352561


halfsine(y, y_opt, width) = \
        abs((y - y_opt) / (width * fwhm_gaussian / fwhm_halfsinusodial)) <= pi / 2.0 ? \
        cos((y - y_opt) / (width * fwhm_gaussian / fwhm_halfsinusodial)) : \
        0


set samples 1023
set xlabel "Y"
set xtics 0.2
set ytics 0.2


load DATA_DIR . "/colors.gp"


plot [y = 0:1] \
     halfsine(y, default_optimum, 0.1) title "width = 0.1", \
     halfsine(y, default_optimum, 0.2) title "default width = 0.2", \
     halfsine(y, default_optimum, 0.4) title "width = 0.4"
