#!/usr/bin/gnuplot -persist

#hinit
do for [ff in system("ls ./output_steady/hinit*")] {
   se size ratio -1;
   end=strlen(ff);
   se te pngcai font "Helvetica, 15";
   set output "./output_change_ela/hinit".ff[end-7:].".png"
   se xla "Easting (100 m)";
   se yla "Northing (100 m)";
   se key above;
   p ff u 1:2:($3/1000.) w p pt 5 pale ti ff[end-7:]
}

#hfinal
do for [ff in system("ls ./output_change_ela/hfinal*")] {
   se size ratio -1;
   end=strlen(ff);
   se te pngcai font "Helvetica, 15";
   set output "./output_change_ela/hfinal".ff[end-7:].".png"
   se xla "Easting (100 m)";
   se yla "Northing (100 m)";
   se key above;
   p ff u 1:2:($3/1000.) w p pt 5 pale ti ff[end-7:]
}

#v and a
reset;
do for [ff in system("ls ./output_change_ela/va*")] {
   end=strlen(ff);
   se te pngcai font "Helvetica, 15";
   set output "./output_change_ela/v".ff[end-7:].".png"
   se xla "t (yr)";
   se yla "volume (km^3)";
   se key above;
   p ff u 1:2 w l lw 2  ti ff[end-7:]
   set output "./output_change_ela/a".ff[end-7:].".png"
   se xla "t (yr)";
   se yla "area (km^2)";
   se key above;
   p ff u 1:3 w l lw 2  ti ff[end-7:]
}


#conservation
do for [ff in system("ls ./output_change_ela/consrv*")] {
   end=strlen(ff);
   se te pngcai font "Helvetica, 15";
   set output "./output_change_ela/consrv".ff[end-7:].".png"
   se ke le;
   se xla "time (yr)";
   se yla "km^3";
   se xra [0:1000]
   p ff u 1:($2*1e-5) w l lw 2 ti "cum. acc.",''  u 1:(($4-$3)*1e-5) w l lw 2 ti "cum. abl.", ''  u 1:($5*1e-5) w l lw 2 ti "vol";
}