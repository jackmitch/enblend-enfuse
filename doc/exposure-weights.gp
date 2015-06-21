default_optimum = 0.5
default_width = 0.2

fwhm_gaussian = 2.3548200450309493820231386529193992755
fwhm_lorentzian = 2.8284271247461900976033774484193961571
fwhm_halfsinusodial = 2.0943951023931954923084289221863352561
fwhm_fullsinusodial = pi
fwhm_power = 1.6817928305074290860622509524664297901


gaussian(y, y_opt, width) = exp(-0.5 * ((y - y_opt) / width)**2)
lorentzian(y, y_opt, width) = 1.0 / (1.0 + ((y - y_opt) / (width * fwhm_gaussian / fwhm_lorentzian))**2 / 2.0)
halfsine(y, y_opt, width) = \
        abs((y - y_opt) / (width * fwhm_gaussian / fwhm_halfsinusodial)) <= pi / 2.0 ? \
        cos((y - y_opt) / (width * fwhm_gaussian / fwhm_halfsinusodial)) : \
        0
fullsine(y, y_opt, width) = \
        abs((y - y_opt) / (width * fwhm_gaussian / fwhm_fullsinusodial)) <= pi ? \
        (1.0 + cos((y - y_opt) / (width * fwhm_gaussian / fwhm_fullsinusodial))) / 2.0 : \
        0
power(y, y_opt, width) = \
        abs((y - y_opt) / (width * fwhm_gaussian / fwhm_power)) <= 1.0 ? \
        1.0 - (abs(y - y_opt) / (width * fwhm_gaussian / fwhm_power))**4 : \
        0.0


set grid ytics
set samples 1023
set xlabel "Y"
set ylabel "w"
set xtics 0.2
set ytics 0.5


load DATA_DIR . "/colors.gp"


plot [y = 0:1] \
     gaussian(y, default_optimum, default_width) title "Gaussian", \
     lorentzian(y, default_optimum, default_width) title "Lorentzian", \
     halfsine(y, default_optimum, default_width) title "Half-Sine", \
     fullsine(y, default_optimum, default_width) title "Full-Sine", \
     power(y, default_optimum, default_width) title "Bi-Square"
