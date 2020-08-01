#!/usr/bin/gnuplot -persist
reset;
se te pngcai fontscale 1.1 size 640,480 font "Helvetica, 15" dashed;
se output "ice.png";
set xrange[2899:];
p 'ice.txt' u 1:5 w l lw 1 ti '14.11709';
q;
