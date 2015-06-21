# Plot an example of a particularly sharp edge


set contour
set grid
set hidden3d
set style data line
set view 60, 30
set xtics 0, 1, 4
set ytics 0, 1, 4


splot DATA_DIR . "/sharp-edge.data" matrix title "sharp\\_edge"
