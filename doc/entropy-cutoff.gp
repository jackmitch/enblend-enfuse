# Plot example of entropy-cutoff function

load "config.gp"

Lower = 0.05
Upper = 0.9
EntropyCutoff(Y, LowerCutoff, UpperCutoff) = \
    Y <= LowerCutoff ? 0.0 : (Y >= UpperCutoff ? 1.0 : Y)

set grid
set key right bottom
set sample 1024
set xlabel "Y"
set xtics 0.2
set ylabel "EC"
set yrange [-0.1:1.1]
set ytics 0.2

set terminal unknown
plot [Y = 0:1] Y, EntropyCutoff(Y, Lower, Upper)

set output "entropy-cutoff.txt"
set terminal dumb Dumb_Width Dumb_Height enhanced
replot

set output "entropy-cutoff.png"
set terminal png font FreeSans 10 size Png_Width, Png_Height enhanced
replot

set output "entropy-cutoff.eps"
set terminal postscript eps enhanced
replot

# Newer Gnuplots have a "pdf" terminal.
set output "| ps2pdf -dEPSCrop - entropy-cutoff.pdf"
set terminal postscript eps enhanced
replot
