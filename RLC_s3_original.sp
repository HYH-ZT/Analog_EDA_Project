* Example RLC-3 - Flat version
Vin 1 0 SIN 0 3 200 0

* First segment (original X1)
R1 1 5 3.5
L1 5 2 1.2e-3  
C1 2 0 7.3e-6

* Second segment (original X2)
R2 2 6 3.5
L2 6 3 1.2e-3
C2 3 0 7.3e-6

* Third segment (original X3)
R3 3 7 3.5
L3 7 4 1.2e-3
C3 4 0 7.3e-6

* Load capacitor
Cl 4 0 10e-6

.shooting 200
*.DC
*.TRAN 0.004e-3 20e-3
*.PROBE V(4)
*.print V(4)
.PLOTNV 4 2 1
.PLOTNV I(Vin)
.end