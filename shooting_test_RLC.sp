* Driven RLC (guaranteed oscillation)

VIN 1 0 SIN 0 1 0.1591549431 0
* FREQ = 1/(2π) ≈ 0.1591549431 Hz  -> ω = 2πF = 1 rad/s

R1 1 2 1
L1 2 3 1
C1 3 0 1

.TRAN 0.01 20

.PLOTNV 1
.PLOTNV 2
.PLOTNV 3
.PLOTNV I(VIN)