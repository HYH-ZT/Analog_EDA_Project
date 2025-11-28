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

    solver sol(ckt, analyses[0]);
    sol.DC_solve();
}