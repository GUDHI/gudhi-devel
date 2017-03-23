import gudhi

periodic_cc = gudhi.PeriodicCubicalComplex(perseus_file='../3d_torus.txt')
diag = periodic_cc.persistence()
gudhi.plot_persistence_barcode(diag)
