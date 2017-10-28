import gudhi

periodic_cc = gudhi.PeriodicCubicalComplex(perseus_file=gudhi.__root_source_dir__ + \
    '/data/bitmap/3d_torus.txt')
diag = periodic_cc.persistence()
plt = gudhi.plot_persistence_barcode(diag)
plt.show()
