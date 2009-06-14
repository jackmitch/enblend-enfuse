# Plot entropy with respect to base 2

load "config.gp"

Base = 2.0
lb(x) = log(x) / log(Base)
Entropy(x) = -x * lb(x) - (1.0 - x) * lb(1.0 - x)

set xlabel "p"
set xtics 0.2
set ylabel "H_2(p)"
set ytics 0.2

set terminal unknown
plot [p = 0:1] Entropy(p)

set output "entropy.txt"
set terminal dumb Dumb_Width Dumb_Height enhanced
replot

set output "entropy.png"
set terminal png font FreeSans 10 size Png_Width, Png_Height enhanced
replot

set output "entropy.eps"
set terminal postscript eps enhanced
replot

# Newer Gnuplots have a "pdf" terminal.
set output "| ps2pdf -dEPSCrop - entropy.pdf"
set terminal postscript eps enhanced
replot
