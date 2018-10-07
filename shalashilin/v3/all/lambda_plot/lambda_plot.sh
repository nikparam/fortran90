set title 'lamda plot at time = 4.0000'
set xlabel '||S*dD/dt+i(H-i*tau)D||^2'
set ylabel '||dD/dt||^2'
FILES=system("ls *.out")
plot for [i=1:words(FILES)] word(FILES,i) u 3:2 every ::6440::6440 title word(FILES,i) 
