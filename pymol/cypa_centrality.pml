# load files
load ../wt/centrality/wt_betweenness.pdb
load ../wt/centrality/wt_closeness.pdb
load ../wt/centrality/wt_degree.pdb
load ../wt/centrality/wt_eigenvector.pdb
load ../s99t/centrality/s99t_betweenness.pdb
load ../s99t/centrality/s99t_closeness.pdb
load ../s99t/centrality/s99t_degree.pdb
load ../s99t/centrality/s99t_eigenvector.pdb

# Change thickness according to b factor value
hide everything
show cartoon
cartoon putty, all

# Change background color
bg_color white

# Selections
select active_site, resi 55+60+102+113+122+126
select mutants, resi 6+29+66+99
select imp_res, resi 61+63+98+101+103+109

# color selections
color red, active_site
color orange, mutants
color blue, imp_res

# label residues
#label (active_site and name ca), resi
#label (mutants and name ca), resi
#label (imp_res and name ca), resi
label name ca, resi

#make compatible
#set pse_export_version, 1.7