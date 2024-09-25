reset
ylabel_offset_x = -15
ylabel_offset_y = 0
xlabel_offset_x = 0
xlabel_offset_y = -5
xtics_offset_x = 0
xtics_offset_y = -2
ytics_offset_x = -2
ytics_offset_y = 0
key_offset_x = -2
key_offset_y = 0
line_width = 15
p0 = "Bold Bold, 32"
p1 = "Bold Bold, 40"
set term eps size 14,10
set out "runTimes.eps"
set size ratio 1/1.4
set margins screen 0.15, screen 0.98, screen 0.95, screen 0.15
set boxwidth 0.125 absolute
set border linewidth 8
set key font p0
set ytics offset ytics_offset_x, ytics_offset_y
set xtics offset xtics_offset_x, xtics_offset_y
set xtics 1
set format x "2^{%.0f}"
off_value = 0#0.125/1.5

titleformat="%25.0f"
set tics font p0
set key title "L_1, n=38" font p0 offset key_offset_x, key_offset_y
set xlabel "Number of threads in block" font p1 offset xlabel_offset_x, xlabel_offset_y
set ylabel "Time (seconds)" font p1 offset ylabel_offset_x, ylabel_offset_y
plot [log(6)/log(2):log(1536)/log(2)] './grid_8192_order_1.txt' u (log($2)/log(2)-off_value):5:4:8:7 w candlesticks lw line_width lc 1 t sprintf(titleformat, 8192), \
'./grid_16384_order_1.txt' u (log($2)/log(2)+off_value):5:4:8:7 w candlesticks lw line_width lc 3 t sprintf(titleformat, 16384)

set out
reset
set term qt
