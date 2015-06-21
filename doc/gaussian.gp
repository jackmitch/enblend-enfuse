# Plot Gauss curve with default parameters of Enfuse


default_optimum = 0.5
default_width = 0.2

gaussian(y, y_opt, width) = exp(-0.5 * ((y - y_opt) / width)**2)


set samples 1023
set xlabel "Y"
set xtics 0.2
set ytics 0.2


load DATA_DIR . "/colors.gp"


plot [y = 0:1] \
     gaussian(y, default_optimum, 0.1) title "width = 0.1", \
     gaussian(y, default_optimum, 0.2) title "default width = 0.2", \
     gaussian(y, default_optimum, 0.4) title "width = 0.4"
