import gudhi

point_cloud = gudhi.read_off(off_file=gudhi.__root_source_dir__ + \
    '/data/points/tore3D_1307.off')
rips_complex = gudhi.RipsComplex(points=point_cloud, max_edge_length=0.3)
simplex_tree = rips_complex.create_simplex_tree(max_dimension=3)
diag = simplex_tree.persistence()
plt = gudhi.plot_persistence_diagram(diag)
plt.show()
