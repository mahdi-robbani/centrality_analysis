 
# First step in the Pyinteraph pipeline: construction of the initial networks for each frame in the MD simulation.

model_path="model0.pdb"
trajectory_path="center_traj.xtc"

# Run Pyinteraph
# Use all residues except the glycines
pyinteraph -f -s $model_path -t $trajectory_path -r $model_path --hc-graph hc_graph.dat --hc-residues ALA,ARG,ASN,ASP,CYS,GLN,GLU,HIS,ILE,LEU,LYS,MET,PHE,PRO,SER,THR,TRP,TYR,VAL
