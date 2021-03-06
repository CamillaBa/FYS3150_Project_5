# FYS3150_Project_5


This repository contains 11 header files,
1 cpp file and 1 python file. The remaining files are data files for plots.

The header files are:

* vector.h                          (A vector class made for simple 1-dim arrays)
* matrix.h                          (A matrix class made for simple 2-dim arrays)
* print_vector_matrix.h             (File containind various print functions)
* tridiag_solver.h                  (Solves a special tridiag matrix eq)

* one_dim_diffusion_eq.h            (A class to represent a spacial solution at a given time)
* one_dim_diffusion_eq_heat_prod.h  (Similar, but with functional heat prod argument)
* two_dim_diffusion_eq.h            (A class to represent a spacial solution at a given time)
* two_dim_diffusion_eq_heat_prod.h  (Similar, but with functional heat prod argument)

* physics.h                         (Various functions for the refertilization problem)
* unit_tests.h                      (Various unit tests)

Then there are "Header.h" and "Main.cpp". 
The file "Header.h" simply includes all of the above header files. 
"Header.h" is included to "Main.cpp".
In "Main.cpp", you will find code to generate the data used for plotting.
The plots can be generated by "plots.py".

When the file "plots.py" is executed as is, with
plot data contained in the folder named data,
the plots used in the report will be generated.

When the file "Main.cpp" is executed, it produce the data in the folder,
and output:

"Tridiagonal matrix solver passes unit test.
Completed iteration: 0
Completed iteration: 10000
Completed iteration: 0
Completed iteration: 10000
Completed iteration: 20000
Completed iteration: 30000
Completed iteration: 40000
Completed iteration: 50000
Completed iteration: 60000
Completed iteration: 70000
Completed iteration: 80000
Completed iteration: 90000
Completed iteration: 100000
Completed iteration: 110000
Completed iteration: 120000
Completed iteration: 130000
Completed iteration: 140000
Completed iteration: 150000
Completed iteration: 160000
Completed iteration: 170000
Completed iteration: 180000
Completed iteration: 190000
Completed iteration: 200000
Completed iteration: 210000
Completed iteration: 220000
Completed iteration: 230000
Completed iteration: 240000
Completed iteration: 250000
Completed iteration: 260000
Completed iteration: 270000
Completed iteration: 280000
Completed iteration: 290000
Completed iteration: 300000
Completed iteration: 310000
Completed iteration: 320000
Completed iteration: 330000
Completed iteration: 340000
Completed iteration: 350000
Completed iteration: 360000
Completed iteration: 370000
Completed iteration: 380000
Completed iteration: 390000
Completed iteration: 400000
Completed iteration: 410000
Completed iteration: 420000
Completed iteration: 430000
Completed iteration: 440000
Completed iteration: 450000
Completed iteration: 460000
Completed iteration: 470000
Completed iteration: 480000
Completed iteration: 490000
Completed iteration: 500000
Completed iteration: 510000
Completed iteration: 520000
Completed iteration: 530000
Completed iteration: 540000
Completed iteration: 550000
Completed iteration: 560000
Completed iteration: 570000
Completed iteration: 580000
Completed iteration: 590000
Completed iteration: 600000
Completed iteration: 610000
Completed iteration: 620000
Completed iteration: 630000
Completed iteration: 640000
Completed iteration: 650000
Completed iteration: 660000
Completed iteration: 670000
Completed iteration: 680000
Completed iteration: 690000
Completed iteration: 700000
Completed iteration: 710000
Completed iteration: 720000
Completed iteration: 730000
Completed iteration: 740000
Completed iteration: 750000
Completed iteration: 760000
Completed iteration: 770000
Completed iteration: 780000
Completed iteration: 790000
Completed iteration: 800000
Completed iteration: 810000
Completed iteration: 820000
Completed iteration: 830000
Completed iteration: 840000
Completed iteration: 850000
Completed iteration: 860000
Completed iteration: 870000
Completed iteration: 880000
Completed iteration: 890000
Completed iteration: 900000
Completed iteration: 910000
Completed iteration: 920000
Completed iteration: 930000
Completed iteration: 940000
Completed iteration: 950000
Completed iteration: 960000
Completed iteration: 970000
Completed iteration: 980000
Completed iteration: 990000
Steady state unit test passed: 1 dimensional heat equation with heat production
agrees with analytical steady state solution. Rel error:0.00066296
Completed iteration: 0
Completed iteration: 10000
Completed iteration: 20000
Completed iteration: 30000
Completed iteration: 40000
Completed iteration: 50000
Completed iteration: 60000
Completed iteration: 70000
Completed iteration: 80000
Completed iteration: 90000
Completed iteration: 100000
Completed iteration: 110000
Completed iteration: 120000
Completed iteration: 130000
Completed iteration: 140000
Completed iteration: 150000
Completed iteration: 160000
Completed iteration: 170000
Completed iteration: 180000
Completed iteration: 190000
Completed iteration: 200000
Completed iteration: 210000
Completed iteration: 220000
Completed iteration: 230000
Completed iteration: 240000
Completed iteration: 250000
Completed iteration: 260000
Completed iteration: 270000
Completed iteration: 280000
Completed iteration: 290000
Completed iteration: 300000
Completed iteration: 310000
Completed iteration: 320000
Completed iteration: 330000
Completed iteration: 340000
Completed iteration: 350000
Completed iteration: 360000
Completed iteration: 370000
Completed iteration: 380000
Completed iteration: 390000
Completed iteration: 400000
Completed iteration: 410000
Completed iteration: 420000
Completed iteration: 430000
Completed iteration: 440000
Completed iteration: 450000
Completed iteration: 460000
Completed iteration: 470000
Completed iteration: 480000
Completed iteration: 490000
Completed iteration: 500000
Completed iteration: 510000
Completed iteration: 520000
Completed iteration: 530000
Completed iteration: 540000
Completed iteration: 550000
Completed iteration: 560000
Completed iteration: 570000
Completed iteration: 580000
Completed iteration: 590000
Completed iteration: 600000
Completed iteration: 610000
Completed iteration: 620000
Completed iteration: 630000
Completed iteration: 640000
Completed iteration: 650000
Completed iteration: 660000
Completed iteration: 670000
Completed iteration: 680000
Completed iteration: 690000
Completed iteration: 700000
Completed iteration: 710000
Completed iteration: 720000
Completed iteration: 730000
Completed iteration: 740000
Completed iteration: 750000
Completed iteration: 760000
Completed iteration: 770000
Completed iteration: 780000
Completed iteration: 790000
Completed iteration: 800000
Completed iteration: 810000
Completed iteration: 820000
Completed iteration: 830000
Completed iteration: 840000
Completed iteration: 850000
Completed iteration: 860000
Completed iteration: 870000
Completed iteration: 880000
Completed iteration: 890000
Completed iteration: 900000
Completed iteration: 910000
Completed iteration: 920000
Completed iteration: 930000
Completed iteration: 940000
Completed iteration: 950000
Completed iteration: 960000
Completed iteration: 970000
Completed iteration: 980000
Completed iteration: 990000
Success!"
