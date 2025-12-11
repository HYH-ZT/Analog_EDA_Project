#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "circuit.hpp"   // 包含 circuit/device/model 定义
#include "parse_netlist.hpp"
#include "solver.hpp"
#include <map>

using namespace std;

int main_dc(){
    circuit ckt;
    vector<string> netlist;
    vector<analysis> analyses;
    //parseNetlistFile("resistor_net.sp", ckt, analyses); //测试电阻网络
    parseNetlistFile("dbmixer.sp", ckt, analyses); 
    //Debug输出解析结果
    std::cout << "Parsed Circuit:\n";
    std::cout << "Nodes:\n";
    for (const auto& pair : ckt.node_map) {
        std::cout << "  " << pair.first << " -> " << pair.second << "\n";
    }
    std::cout << "Devices:\n";
    for (const auto& dev : ckt.linear_devices) {
        std::cout << "  " << dev.name << " (" << dev.type << ") Nodes: ";
        for (const auto& node_name : dev.node_names) {
            std::cout << node_name << " ";
        }
        std::cout << " Parameters: ";
        for (const auto& param : dev.parameters) {
            std::cout << param.first << "=" << param.second << " ";
        }
        std::cout << "\n";
    }
    for (const auto& dev : ckt.nonlinear_devices) {
        std::cout << "  " << dev.name << " (" << dev.type << ") Nodes: ";
        for (const auto& node_name : dev.node_names) {
            std::cout << node_name << " ";
        }
        std::cout << " Parameters: ";
        for (const auto& param : dev.parameters) {
            std::cout << param.first << "=" << param.second << " ";
        }
        std::cout << "\n";
    }
    for (const auto& dev : ckt.sources) {
        std::cout << "  " << dev.name << " (" << dev.type << ") Nodes: ";
        for (const auto& node_name : dev.node_names) {
            std::cout << node_name << " ";
        }
        std::cout << " Parameters: ";
        for (const auto& param : dev.parameters) {
            std::cout << param.first << "=" << param.second << " ";
        }
        std::cout << "\n";
    }
    //输出打印变量
    std::cout << "Analyses:\n";
    for (const auto& analysis : analyses) {
        std::cout << "  Type: " << analysis.type << "\n";
        std::cout << "  Print Variables: ";
        for (const auto& var : analysis.print_variables) {
            std::cout << var << " ";
        }
        std::cout << "\n";
    }

    solver sol(ckt, analyses[0]);

    // //设置hb参数
    // sol.hb_params.fundamental_omega = 1; // 1MHz
    // sol.hb_params.num_harmonics = 5; // 30个谐波
    // sol.hb_params.max_iterations = 50;
    // sol.hb_params.tolerance = 1e-6;
    // //设置基频初始解 //如果要运行PSS，必须要用这个接口设置初始解，因为其中确定了Liner_Y的大小
    // sol.HB_set_initial_xw({{"102",1},{"103",-1}}); //空参数表示使用默认初始解
    //根据网表设定参数
    sol.hb_params.fundamental_omega = 2 * 3.14159265358979323846 * analyses[0].parameters["freq"]; // 基频角频率
    sol.hb_params.num_harmonics = static_cast<int>(analyses[0].parameters["harm"]); // 谐波数量
    sol.hb_params.max_iterations = 50; // 最大迭代次数
    sol.hb_params.tolerance = 1e-6; // 收敛容限
    sol.hb_params.relaxation_factor = 1; // 松弛因子

    ckt.extract_MOS_capacitances(); //提取MOS管寄生电容

    //Debug输出提取寄生电容后的线性器件
    std::cout << "Linear Devices after Extracting MOS Capacitances:\n";
    for (const auto& dev : ckt.linear_devices) {
        std::cout << "  " << dev.name << " (" << dev.type << ") Nodes: ";
        for (const auto& node_name : dev.node_names) {
            std::cout << node_name << " ";
        }
        std::cout << " Parameters: ";
        for (const auto& param : dev.parameters) {
            std::cout << param.first << "=" << param.second << " ";
        }
        std::cout << "\n";
    }

    sol.HB_set_initial_xw({}); //空参数表示使用默认初始解

    sol.PSS_solve_harmonic_balance();

    return 0;
}