import Traj_analyze_fracoccupied as main
import sys
N = int(sys.argv[1])
Rval = int(sys.argv[2])
Nreps = int(sys.argv[3])
main.frontend_cross_set_Rval(N, Rval, Nreps)
