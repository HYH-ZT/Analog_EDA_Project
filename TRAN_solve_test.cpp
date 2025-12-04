#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "circuit.hpp"   // 包含 circuit/device/model 定义
#include "parse_netlist.hpp"
#include "solver.hpp"
#include <map>

using namespace std;

int main(){
    circuit ckt;
    vector<string> netlist;
    vector<analysis> analyses;
    //parseNetlistFile("resistor_net.sp", ckt, analyses); //测试电阻网络
    parseNetlistFile("buffer.sp", ckt, analyses); 
    // 使用节点名和电压值的映射设置初值

    solver sol(ckt, analyses[0]);
    sol.TRAN_solve();
}