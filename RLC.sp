* Simple Resistor Network
* 简单电阻网络测试电路

* 电源定义
VDD 101 0 DC 5
* V1 101 0 SIN 5 1 0.159155 0

* 电阻网络
* R1 101 102 1
L1 101 102 0.5
L2 101 103 0.5
C1 102 103 1
M1   102 103 0   n 10e-6 0.35e-6 2
M2   103 102 0   n 10e-6 0.35e-6 2
* C1 102 103 1
* L1 103 0 1
* R1 102 0 1

.MODEL 1 VT -0.75 MU 5e-2 COX 0.3e-4 LAMBDA 0.05 CJ0 4.0e-14
.MODEL 2 VT 0.83 MU 1.5e-1 COX 0.3e-4 LAMBDA 0.05 CJ0 4.0e-14


* 分析类型（根据需要选择）
.hb 1 3
.print hb V(101) V(102)

