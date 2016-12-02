# Planning of the Belgium gas transmission
#########################################################
# Sets
#########################################################

set TOWNS;                           # nodes
set ARCS within {TOWNS,TOWNS};       # pipelines 
set ARCSA within ARCS;               # active pipelines

#########################################################
# Parameters
#########################################################


########### town parameters #############################
param slo {TOWNS}; # inflow lower bound (mill M3 per day)
param sup {TOWNS}; # inflow upper bound (mill M3 per day)
param plo {TOWNS}; # pressure lower bound (bar)
param pup {TOWNS}; # pressure upper bound (bar)
param c {TOWNS};   # purchase price of gas ($ per MBTU)

########## arc parameters ###############################
param D {ARCS};          # diameter (mm)
param L {ARCS};          # length   (km)
param T;   # gas temperature (K)
param e;   # absolute rugosity (mm)
param den; # density of gas realtive to air (-) 
param z;   # compessibility factor (-)

param lambda {(i,j) in ARCS} = 
	1/(2*log10(3.7*D[i,j]/e))^2; 
	# lambda parameter
param c2 {(i,j) in ARCS} = 
	96.074830e-15*D[i,j]^5/lambda[i,j]/z/T/L[i,j]/den; 
	# C^2 parameter  


#########################################################
# Variables
#########################################################

var f {TOWNS,TOWNS};         # arc flow (mill M3 per day)
var s {TOWNS};               # supply - demand (mill M3 per day)
var pi {TOWNS};              # pressure^2 (bar^2)

#########################################################
# Objective function
#########################################################

# minimize the total cost of supplies
minimize total_cost: sum{j in TOWNS} c[j]*(s[j]*10^6/28.31685);

#########################################################
# Constraints
#########################################################

# flow balance
subject to Flow_Balance {i in TOWNS}: 
sum{j in TOWNS}f[i,j] = s[i] + sum{j in TOWNS}f[j,i];


# disallow flow on non-existing arcs
subject to Flow_Feasibility {(i,j) in {TOWNS,TOWNS} diff ARCS}:
f[i,j] = 0;

# passive arc equilibrium
subject to Passive_Equilibrium {(i,j) in ARCS diff ARCSA}:
abs(f[i,j])*f[i,j] = c2[i,j]*(pi[i]-pi[j]);

# active arc equilibrium
subject to Active_Equilibrium {(i,j) in ARCSA}:
abs(f[i,j])*f[i,j] >= c2[i,j]*(pi[i]-pi[j]);


# fulfill demands
subject to Demand {i in TOWNS}:
slo[i] <= s[i] <= sup[i];

# fulfill pressure requirements
subject to Pressure {i in TOWNS}:
plo[i]^2 <= pi[i] <= pup[i]^2;


# active arcs requirement
subject to Active_Arc {(i,j) in ARCSA}:
f[i,j] >= 0;










