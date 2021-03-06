# load files
load ../wt/centrality/wt_betweenness.pdb
load ../wt/centrality/wt_closeness.pdb
load ../wt/centrality/wt_degree.pdb
load ../wt/centrality/wt_eigenvector.pdb
load ../s99t/centrality/s99t_betweenness.pdb
load ../s99t/centrality/s99t_closeness.pdb
load ../s99t/centrality/s99t_degree.pdb
load ../s99t/centrality/s99t_eigenvector.pdb
load ../wt/cc.pdb

# Align
align wt_betweenness, s99t_betweenness
align wt_closeness, s99t_closeness
align wt_degree, s99t_degree
align wt_eigenvector, s99t_eigenvector
align cc, s99t_betweenness

# Change thickness according to b factor value
hide everything
spectrum b, blue_red, wt_betweenness, 0, 0.20
spectrum b, blue_red, wt_closeness, 0, 0.07
spectrum b, blue_red, wt_degree, 0, 0.031
spectrum b, blue_red, wt_eigenvector, 0, 0.52

spectrum b, blue_red, s99t_betweenness, 0, 0.20
spectrum b, blue_red, s99t_closeness, 0, 0.07
spectrum b, blue_red, s99t_degree, 0, 0.031
spectrum b, blue_red, s99t_eigenvector, 0, 0.52

spectrum b, rainbow, cc
as cartoon
#show cartoon
#cartoon putty, all

# Change background color
bg_color white

# Selections
#select degree, resi 82+102+111 and (wt_degree or s99t_degree)
#select betweenness, resi 99+100+108+113+127 and (wt_betweenness or s99t_betweenness)
#select closeness, resi 66+84+125 and (wt_closeness or s99t_closeness)
#select eigenvector, resi 32+85+86+102+108+113+127 and (wt_eigenvector or s99t_eigenvector)

# new selections
select betweenness, resi 113+99+92+122+126+102+108+100+32+129+98+22+132+136+139+142+156+12 and (wt_betweenness or s99t_betweenness)
select closeness, resi 100+108+32 and (wt_closeness or s99t_closeness)
select degree, resi 38+32+92+142 and (wt_degree or s99t_degree)
select eigenvector, resi 92+119+128+122+99 and (wt_eigenvector or s99t_eigenvector)

select imp, resi 32+108+92
select active_site, resi 55+60+102+113+122+126

# show stick
#show line, active_site
show stick, degree
show stick, betweenness
show stick, closeness
show stick, eigenvector
show stick, imp


# color selections
#color green, degree
#color green, betweenness
#color green, closeness
#color green, eigenvector
#color green, active_site

# label residues
label (name ca and resi 10+30+50+70+90+110+130+150), resi

label (betweenness and name ca), resi
label (closeness and name ca), resi
label (degree and name ca), resi
label (eigenvector and name ca), resi

label (imp and name ca), "%s-%s" % (resn,resi)
label (active_site and name ca), "%s-%s" % (resn,resi)

# label color
#set label_color, green, (betweenness and name ca)
#set label_color, green, (closeness and name ca)
#set label_color, green, (degree and name ca)
#set label_color, green, (eigenvector and name ca)
set label_color, green, (imp and name ca)
set label_color, cyan, (active_site and name ca)


#label attributes
set label_size, 25
set label_font_id, 7

#make compatible
#set pse_export_version, 1.7