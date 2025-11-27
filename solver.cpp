#include "solver.hpp"
#include <Eigen/Dense>
#include "circuit.hpp"
solver::solver(circuit& ckt_, analysis& analysis_):ckt(ckt_), analysis_type(analysis_){
    liner_Y.resize(ckt.node_map.size() - 1, ckt.node_map.size() - 1); //不包含地节点
    liner_Y.setZero();
    J.resize(ckt.node_map.size() - 1);
    J.setZero();
}

//直流分析
void solver::DC_solve(){
    //构建直流MNA矩阵
    //构建线性器件的MNA矩阵
    for (auto &dev : ckt.linear_devices){
        char c = toupper(dev.name[0]);
        //电阻
        if (c == 'R'){
            int n1 = dev.nodes[0];
            int n2 = dev.nodes[1];
            double value = dev.parameters["VALUE"];
            double conductance = 1.0 / value;
            if (n1 != 0){
                liner_Y(n1 - 1, n1 - 1) += conductance; //节点编号从1开始，矩阵索引从0开始
            }
            if (n2 != 0){
                liner_Y(n2 - 1, n2 - 1) += conductance;
            }
            if (n1 != 0 && n2 != 0){
                liner_Y(n1 - 1, n2 - 1) -= conductance;
                liner_Y(n2 - 1, n1 - 1) -= conductance;
            }
        }
        if (c == 'C'){
            //电容断路
        }
        if (c == 'L'){
            //电感短路，视为电压为0的电压源
            device source;
            source.name = "L_short_" + dev.name;
            source.type = "V_DC";
            source.nodes = dev.nodes;
            source.parameters["DC"] = 0.0;
            ckt.sources.push_back(source);
        }
    }
    MNA_Y = liner_Y;
    for (auto &dev : ckt.nonlinear_devices){
        char c = toupper(dev.name[0]);
        //MOS管
        if (c == 'M'){
            int d = dev.nodes[0];
            int g = dev.nodes[1];
            int s = dev.nodes[2];
            //简单的DC工作点假设：Vgs=2V, Vds=2V, kn=0.5mA/V^2, Vth=1V
            double Vg = (g == 0) ? 0.0 : node_voltages[g - 1];
            double Vs = (s == 0) ? 0.0 : node_voltages[s - 1];
            double Vd = (d == 0) ? 0.0 : node_voltages[d - 1];
            double Vgs = Vg - Vs;
            if ()
            double kn = dev.parameters["KP"]; //单位A/V^2
            double Vth = dev.parameters["VTO"]; //单位V
            double Ids = 0.0;
            if (Vgs <= Vth){
                Ids = 0.0; //截止区
            }
            else if (Vds < (Vgs - Vth)){
                Ids = kn * ((Vgs - Vth) * Vds - 0.5 * Vds * Vds); //线性区
            }
            else{
                Ids = 0.5 * kn * (Vgs - Vth) * (Vgs - Vth); //饱和区
            }
            if (d != 0){
                J(d - 1) -= Ids; //流出漏极为负
            }
            if (s != 0){
                J(s - 1) += Ids; //流入源极为正
            }
        }
    }
    for (auto &dev : ckt.sources){
        char c = toupper(dev.name[0]);
        //独立电压源
        if (c == 'V'){
            int n1 = dev.nodes[0];
            int n2 = dev.nodes[1];
            double value = dev.parameters["DC"];
            //引入辅助变量和KVL方程
            if (n1 == 0){
                //节点n2电压为-value
                J(n2 - 1) = -value;
                MNA_Y.row(n2 - 1).setZero();
                MNA_Y(n2 - 1, n2 - 1) = 1;
            }
            else if (n2 == 0){
                //节点n1电压为value
                J(n1 - 1) = value;
                MNA_Y.row(n1 - 1).setZero();
                MNA_Y(n1 - 1, n1 - 1) = 1;
            }
            else{
                //引入支路电流变量与新方程
                int new_var_index = MNA_Y.rows();
                MNA_Y.conservativeResize(new_var_index + 1, new_var_index + 1);
                MNA_Y.row(new_var_index).setZero();
                MNA_Y.col(new_var_index).setZero();
                J.conservativeResize(new_var_index + 1);
                J(new_var_index) = value;
                //KVL方程
                MNA_Y(new_var_index, n1 - 1) = 1;
                MNA_Y(new_var_index, n2 - 1) = -1;
                //支路电流对节点的贡献
                MNA_Y(n1 - 1, new_var_index) = 1;
                MNA_Y(n2 - 1, new_var_index) = -1;
            }
        }
        if (c == 'I'){
            //独立电流源
            int n1 = dev.nodes[0];
            int n2 = dev.nodes[1];
            double value = dev.parameters["DC"];
            if (n1 != 0){
                J(n1 - 1) -= value; //流出节点为负
            }
            if (n2 != 0){
                J(n2 - 1) += value; //流入节点为正
            }
        }
    }
    //求解直流MNA矩阵方程
}   