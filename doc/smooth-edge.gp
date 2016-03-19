# Plot an example of a smooth edge


set contour
set grid
set hidden3d
set style data line
set view 60, 30
set xtics 0, 1, 4
set ytics 0, 1, 4

unset key


splot DATA_DIR . "/smooth-edge.data" matrix title "smooth edge"
