# Plot an example of a particularly sharp edge

load "config.gp"
load "config-edge.gp"

set terminal unknown
splot "sharp-edge.data" matrix title "sharp\\_edge"

set output "sharp-edge.txt"
set terminal dumb Dumb_Width Dumb_Height enhanced
replot

set output "sharp-edge.png"
set terminal png font FreeSans 10 size Png_Width, Png_Height enhanced
replot

set output "sharp-edge.eps"
set terminal postscript eps enhanced
replot

# Newer Gnuplots have a "pdf" terminal.
set output "| ps2pdf -dEPSCrop - sharp-edge.pdf"
set terminal postscript eps enhanced
replot
