module load lmm-tools
#need a gmx_mpi version
source /usr/local/gromacs-5.1.2_plumed-2.3b/bin/GMXRC.bash

pdb=model0.pdb
trj=center_traj.xtc

t_acc4 -f $trj -s $pdb -S -o acc.xvg

