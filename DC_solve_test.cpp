#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "circuit.hpp"   // 包含 circuit/device/model 定义
#include "parse_netlist.hpp"
#include "solver.hpp"

using namespace std;

int main(){
    circuit ckt;
    vector<string> netlist;
    vector<analysis> analyses;
    //parseNetlistFile("resistor_net.sp", ckt, analyses); //测试电阻网络
    parseNetlistFile("buffer.sp", ckt, analyses); 
    // 使用节点名和电压值的映射设置初值
    std::map<std::string, double> initial_voltages;
    initial_voltages["101"] = 5.0;
    initial_voltages["102"] = 2.0;
    initial_voltages["103"] = 4.0;
    initial_voltages["104"] = 3.0;
    initial_voltages["105"] = 1.0;

    solver sol(ckt, analyses[0]);
    sol.setLinearSolverMethod(LinearSolverMethod::LU_DECOMPOSITION);
    // sol.DC_solve(initial_voltages);
    sol.DC_solve();
    return 0;
}