# Plot an example of a particularly sharp edge

load "config.gp"
load "config-edge.gp"

set terminal unknown
splot "smooth-edge.data" matrix title "smooth\\_edge"

set output "smooth-edge.txt"
set terminal dumb Dumb_Width Dumb_Height enhanced
replot

set output "smooth-edge.png"
set terminal png font FreeSans 10 size Png_Width, Png_Height enhanced
replot

set output "smooth-edge.eps"
set terminal postscript eps enhanced
replot

# Newer Gnuplots have a "pdf" terminal.
set output "| ps2pdf -dEPSCrop - smooth-edge.pdf"
set terminal postscript eps enhanced
replot
