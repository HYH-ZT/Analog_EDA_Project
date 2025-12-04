* Simple Resistor Network
* 简单电阻网络测试电路

* 电源定义
VDD 101 0 DC 5

* 电阻网络
R1 101 102 100
C1 102 0 0.01
* L1 101 102 1
* R1 102 0 1



* 分析类型（根据需要选择）
.tran 0.01 5
.print tran V(102) I(VDD)

