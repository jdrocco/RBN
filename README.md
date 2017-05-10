Finds the dynamical attractors in a randomly initiated boolean 
network [1,2].

RBN.c - Chooses each strategy from all the 2^(2^K) 
possibilities, biased towards a preselected value of p (p = mean 
fraction of ones in the strategies of the various nodes). 
The p values of the strategies will have a normal distribution 
about the controlling p value for each network distribution, 
which itself is chosen at random on the interval [highp,lowp].

Versions which produce a linear interpolation for p, or delta 
functions for p and effective k, are easy to overconstrain.

Input file - RBNparam:

K - degree of network
N - number of nodes
RANDINITS - number of random initial states for each system
DIAGFLAG - most verbose logging
PRINTFLAG - somewhat verbose logging
SYSTEMS - number of system initializations (i.e. strategies & inputs, 
          the "wiring" of the network) 
LOWP - minimum value of p 
HIGHP - maximum value of p
MAXSTEPS - maximum iterations for each random initial state
SAVEINT - interval for saving network state
MAXATTRACTORS - maximum size of attractor library

Three levels of output detail:
1) RBNout - prints a summary of data obtained for each system realization
2) RBNdetail - prints the outcome of each random initial configuration:
	whether an attractor was found, which one, how long it took to find
	and how long the attractor itself is
   Activated by PRINTFLAG in the RBNparam file.
3) print to screen - prints all the inputs, strategies, and states for every
	iteration.  Only do for small simulations.
   Activated by DIAGFLAG in the RBNparam file.

Limitations:
max K = 6
Cannot do a distribution of K except by randomly drawing higher
digit strategies and waiting for one having a lower effective K.

[1] Kauffman, S. A. (1969). Metabolic stability and epigenesis in randomly 
constructed genetic nets. Journal of Theoretical Biology, 22:437-467.
[2] Luque, B. and Sole, R.V. (2000). Lyapunov exponents in random Boolean 
networks. Physica A, 284:33-45.
