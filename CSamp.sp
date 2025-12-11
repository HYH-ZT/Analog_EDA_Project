* 简单共源放大器网络测试电路

* 电源定义
VDD 101 0 DC 5
VINN 102 0 SIN 0 3 0.159155 0
VINP 104 0 SIN 5 3 0.159155 0
VTN 106 0 DC 0
VTP 107 0 DC 0 

* MOS管
M1   103 102 106 n 1e-6 1e-6 2
M2   105 104 101 p 1e-6 1e-6 1

* 电阻网络
R1 101 103 1e5
R2 105 107 1e5

*  MOS模型
.MODEL 1 VT -0.75 MU 5e-2 COX 0.3e-4 LAMBDA 0.00 CJ0 4.0e-14
.MODEL 2 VT 0.83 MU 1.5e-1 COX 0.3e-4 LAMBDA 0.00 CJ0 4.0e-14

* 分析类型（根据需要选择）
.hb 1 5
.print hb V(103) V(105) 

