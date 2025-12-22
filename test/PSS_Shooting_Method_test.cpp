// #include <iostream>
// #include <fstream>
// #include <string>
// #include <vector>
// #include "circuit.hpp"   // 包含 circuit/device/model 定义
// #include "parse_netlist.hpp"
// #include "solver.hpp"

// using namespace std;

// int main_pss_shooting_method_test(){
//     circuit ckt;
//     vector<string> netlist;
//     vector<analysis> analyses;
//     //parseNetlistFile("resistor_net.sp", ckt, analyses); //测试电阻网络
//     parseNetlistFile("buffer.sp", ckt, analyses); 

//     solver sol(ckt, analyses[0]);
//     double frequency = 10e6;
//     double period_T = 1.0 / frequency;
//     sol.PSS_solve_shooting(period_T, period_T / 100, 100, 1e-6);
//     //用稳态分析得到的节点电压值跑一次瞬态
//     cout << "Running transient simulation for one period to verify PSS result...\n";
//     sol.TRAN_solve(period_T, period_T / 100);
//     //Debug:展示节点电压结果
//     cout << "\n\n\n";
//     sol.print_node_voltages();
//     return 0;
// }