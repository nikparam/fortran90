FILES=system("ls *.out")
plot for [i=1:words(FILES)] word(FILES,i) u 1:3 title word(FILES,i) w l
