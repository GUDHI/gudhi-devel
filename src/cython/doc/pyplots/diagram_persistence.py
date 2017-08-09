import gudhi

rips_complex = gudhi.RipsComplex(off_file='tore3D_1307.off', max_edge_length=0.2)
simplex_tree = rips_complex.create_simplex_tree(max_dimension=3)
diag = simplex_tree.persistence()
plt = gudhi.plot_persistence_diagram(diag, band_boot=0.13)
plt.show()
