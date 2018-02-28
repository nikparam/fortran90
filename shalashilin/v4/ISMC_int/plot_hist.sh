#!/bin/bash
binwidth = 0.1
bin( x, width ) = width * floor( x / width )
#gauss(x) = a / ( sqrt( 2 * pi ) * sigma ) * exp( -(x-1)**2 / (2*sigma**2) )
sigma = 1; a = 10;
#fit gauss(x) 'out.out' via sigma, a
plot 'out.out' using ( bin( $1,binwidth) ) : (1.0) smooth freq with boxes #, gauss(x)
