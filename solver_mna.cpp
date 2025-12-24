#include "solver.hpp"
#include "solver_internal.hpp"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

void solver::build_linear_MNA(bool in_tran){
    //构建线性MNA矩阵liner_Y和电流源向量J
    //构建线性器件的MNA矩阵
    //根据节点数量调整liner_Y和J的大小
    int nodeCount = (int)ckt.node_list.size() - 1; //不包含地节点
    liner_Y.resize(nodeCount, nodeCount);
    liner_Y.setZero();
    J.resize(nodeCount);
    J.setZero();
    for (auto &dev : ckt.linear_devices){
        char c = toupper(dev.name[0]);
        //电阻
        if (c == 'R'){
            int n1 = dev.nodes[0];
            int n2 = dev.nodes[1];
            double value = dev.parameters["value"];
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
        if (c == 'L' && in_tran == false){
            //电感短路，视为电压为0的电压源
            device source;
            source.name = "VL_short_" + dev.name;
            source.type = "V_DC";
            source.nodes = dev.nodes;
            source.parameters["DC"] = 0.0;
            source.original_device_name = dev.name; // 标记为电感的等效器件
            //保存指向该电感器件的指针，以便后续更新电流
            source.original_device_ptr = &dev;
            ckt.sources.push_back(source);
        }
    }
}

//MNA矩阵操作辅助函数

void solver::addToY(int rowNode, int colNode, double val) {
    if (rowNode == 0 || colNode == 0) return;
    int r = rowNode - 1;
    int c = colNode - 1;
    
    // 检查是否需要扩展矩阵大小
    int current_rows = MNA_Y.rows();
    int current_cols = MNA_Y.cols();
    int required_size = std::max(r + 1, c + 1);
    
    if (required_size > current_rows || required_size > current_cols) {
        // 需要扩展矩阵
        int new_size = std::max({current_rows, current_cols, required_size});
        
        // 保存原矩阵
        Eigen::MatrixXd old_matrix = MNA_Y;
        
        // 重新调整大小并初始化为零
        MNA_Y.resize(new_size, new_size);
        MNA_Y.setZero();
        
        // 复制原矩阵内容到新矩阵
        if (current_rows > 0 && current_cols > 0) {
            MNA_Y.block(0, 0, current_rows, current_cols) = old_matrix;
        }
    }
    
    MNA_Y(r, c) += val;
}


void solver::addToJ(int node, double val) {
    if (node == 0) return;
    int idx = node - 1;
    
    // 检查是否需要扩展向量大小
    int current_size = J.size();
    int required_size = idx + 1;
    
    if (required_size > current_size) {
        // 需要扩展向量
        // 保存原向量
        Eigen::VectorXd old_vector = J;
        
        // 重新调整大小并初始化为零
        J.resize(required_size);
        J.setZero();
        
        // 复制原向量内容到新向量
        if (current_size > 0) {
            J.head(current_size) = old_vector;
        }
    }
    
    J(idx) += val;
}


void solver::build_nonlinear_MNA() {
    // 假设 MNA_Y 已经被初始化为 liner_Y 的副本（见你 DC_solve 片段）
    // 假设 J 已预先被清零并且大小等于节点数（node_count）
    // node_voltages 存的是上一次迭代的节点电压，索引对应节点编号-1，节点0(ground)不在该向量中

    //const double gmin = 1e-12; // cutoff 时的小导纳容错
    int nodeCount = (int)ckt.node_list.size() - 1; // node_list 含地吗？根据你的mapping调整
    // 注意：如果 node_voltages 的长度不同，请确保索引安全
    // 对每一个 MOS 器件做线性化并 stamp
    for (const auto &dev : ckt.nonlinear_devices) {
        if (dev.type != "MOS") continue;

        // 节点编号（你的 parse 保证 nodes 顺序：0: T1, 1: G, 2: T2）
        int n1= dev.nodes.size() > 0 ? dev.nodes[0] : 0;
        int ng = dev.nodes.size() > 1 ? dev.nodes[1] : 0;
        int n2 = dev.nodes.size() > 2 ? dev.nodes[2] : 0;

        // 读取几何与类型
        double W = dev.parameters.count("W") ? dev.parameters.at("W") : 1e-6;
        double L = dev.parameters.count("L") ? dev.parameters.at("L") : 1e-6;
        double type = dev.parameters.count("TYPE") ? dev.parameters.at("TYPE") : 1.0; // 1: n, -1: p

        // 找到 model 并读取参数（KP, VTO, LAMBDA）
        const model* pmodel = ckt.findModelConst(dev.model);
        double MU;
        double COX;
        double VT;
        double LAMBDA;
        if (pmodel) {
            if (pmodel->parameters.count("MU")) MU = pmodel->parameters.at("MU");
            if (pmodel->parameters.count("VT")) VT = pmodel->parameters.at("VT");
            if (pmodel->parameters.count("COX")) COX = pmodel->parameters.at("COX");
            if (pmodel->parameters.count("LAMBDA")) LAMBDA = pmodel->parameters.at("LAMBDA");
        }
        double KP = MU * COX; // 过程跨导参数
        // beta = KP * (W / L)
        double beta = KP * (W / L);

        // 取得上一次迭代的节点电压（若为地，视为0）
        auto V_of = [&](int node)->double {
            if (node == 0) return 0.0;
            int idx = node - 1;
            if (idx < 0 || idx >= (int)node_voltages.size()) return 0.0;
            return node_voltages[idx];
        };

        double V1 = V_of(n1);
        double Vg0 = V_of(ng);
        double V2 = V_of(n2);

        // 对于 PMOS，把电压和阈值变换为针对模型的常用形式
        // 这里我们通过 TYPE 字段处理：如果 TYPE==-1（p），把电压极性翻转
        double sign = type; // n: +1, p: -1

        //根据电压高低确定 Drain 和 Source, 注意NMOS和PMOS的源漏定义
        int nd, ns;
        double Vd0, Vs0;

        if (V1 > V2 && type > 0 || V1 <= V2 && type < 0) { 
            nd = n1;
            Vd0 = V1;
            ns = n2; 
            Vs0 = V2;
        } else {
            nd = n2; 
            Vd0 = V2;
            ns = n1; 
            Vs0 = V1;
        }

        double Vgs = sign * (Vg0 - Vs0);
        double Vds = sign * (Vd0 - Vs0);
        double Vth = VT * sign; 


        // 计算 Id0, gm, gds, Ieq
        double Id0 = 0.0;
        double gm = 0.0;
        double gds = 0.0;
        double Ieq = 0.0;

        if (Vgs <= Vth) {
            // cutoff
            // 给一个很小的导纳以保证数值稳定性
            gm = 0.0;
            gds = 1e-12;
            Id0 = 0.0;
            Ieq = 0.0;
        } else {
            // 依据 Vds 与 Vgs-Vth 判定工作区
            double Vov = Vgs - Vth; // overdrive
            if (Vds < Vov) {
                // 线性区 (triode)
                // Id = beta * ( (Vov)*Vds - 0.5*Vds^2 )
                Id0 = beta * (Vov * Vds - 0.5 * Vds * Vds);
                // 导数计算
                // gm = ∂Id/∂Vg = beta * Vds
                gm = beta * Vds;
                // gds = ∂Id/∂Vd = beta * (Vov - Vds)
                gds = beta * (Vov - Vds);
            } else {
                // 饱和区 (saturation)
                // Id = 0.5 * beta * Vov^2 * (1 + lambda*Vds)
                // LAMBDA = 0;
                Id0 = 0.5 * beta * Vov * Vov * (1.0 + LAMBDA * Vds);
                // gm = ∂Id/∂Vg = beta * Vov * (1 + lambda*Vds)
                gm = beta * Vov * (1.0 + LAMBDA * Vds);
                // gds = ∂Id/∂Vd = 0.5 * lambda * beta * Vov^2
                gds = 1e-12 + 0.5 * LAMBDA * beta * Vov;
            }
            Ieq = Id0 - gm * Vgs - gds * Vds;
        }

        // 对 PMOS，我们计算的 Id0/gm/gds 是基于 sign * voltages 的正向定义。
        // 物理上上述 Id0 表示从 drain -> source（按 sign 方向）。 为简化，我们保持 Id0 表示流出 drain 到 source（即 drain→source）。
        // 若 TYPE==-1（p），sign=-1 已经在 Vgs/Vds 中处理好，Id0 的符号应是正确的。

        // 现在构造线性方程：小信号 i ≈ gds*(vd - vs) + gm*(vg - vs) + Iconst
        // 求常量 Iconst 使在工作点成立：
        // Iconst = Id0 - gds*Vd0 - gm*Vg0 + (gds+gm)*Vs0

        // --- 把 gm, gds, Iconst stamp 进 MNA_Y 和 J ---
        // 注意：节点编号为 0 表示地。MNA_Y 行列对应节点编号 1..N -> 索引 0..N-1

        //处理gm，受控电流源
        addToY(nd, ng, gm);
        addToY(ns, ng, -gm);
        addToY(nd, ns, -gm);
        addToY(ns, ns, gm);

        //处理gds，电导
        addToY(nd, nd, gds);
        addToY(nd, ns, -gds);
        addToY(ns, nd, -gds);
        addToY(ns, ns, gds);

        //处理Ieq，直流源
        Ieq *= sign;
        addToJ(nd, -Ieq);
        addToJ(ns, Ieq);
    } // end for each MOS
}


void solver::build_sources_MNA(){
    //清空打印支路电流索引
    ckt.print_branch_current_indices.clear();
    ckt.plot_branch_current_indices.clear();
    for (auto &dev : ckt.sources){
        char c = toupper(dev.name[0]);
        //独立电压源
        if (c == 'V'){
            int n1 = dev.nodes[0];
            int n2 = dev.nodes[1];
            double value = dev.parameters["DC"];

            //先引入支路电流变量，从n1流向n2
            int new_var_index = MNA_Y.rows();
            MNA_Y.conservativeResize(new_var_index + 1, new_var_index + 1);
            MNA_Y.row(new_var_index).setZero();
            MNA_Y.col(new_var_index).setZero();
            J.conservativeResize(new_var_index + 1);
            J(new_var_index) = value;
            //记录改电压源对应的支路电流变量索引，从序号0开始
            dev.branch_current_index = new_var_index- ((int)ckt.node_list.size() - 1);
            //如果需要打印支路电流，存在ckt.print_branch_current_indices,需要减去节点数偏移
            if(dev.printI){
                ckt.print_branch_current_indices.push_back(new_var_index - ((int)ckt.node_list.size() - 1));
            }
            if(dev.plotI){
                ckt.plot_branch_current_indices.push_back(new_var_index - ((int)ckt.node_list.size() - 1));
            }
            // 如果是瞬态等效电压源，更新动态器件映射
            if (!dev.original_device_name.empty()) {
                dynamic_device_current_map[dev.original_device_name] = new_var_index- ((int)ckt.node_list.size() - 1);
            }
            //如果是电感的等效电压源，保存指向该电感器件的指针，以便后续更新电流
            if (dev.original_device_ptr) {
                dev.original_device_ptr->branch_current_index = new_var_index- ((int)ckt.node_list.size() - 1);
            }
            //KVL方程
            if (n1 != 0){
                MNA_Y(new_var_index, n1 - 1) = 1;
                MNA_Y(n1 - 1, new_var_index) = 1;
            }
            if (n2 != 0){
                MNA_Y(new_var_index, n2 - 1) = -1;
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
}


void solver::build_sources_MNA_ramp(double alpha){
    // 清空打印支路电流索引
    ckt.print_branch_current_indices.clear();
    ckt.plot_branch_current_indices.clear();

    for (auto &dev : ckt.sources){
        char c = toupper(dev.name[0]);

        // 独立电压源
        if (c == 'V'){
            int n1 = dev.nodes[0];
            int n2 = dev.nodes[1];

            // 斜坡：把源值按 alpha 缩放
            double value = alpha * dev.parameters["DC"];

            // 先引入支路电流变量，从 n1 流向 n2
            int new_var_index = MNA_Y.rows();
            MNA_Y.conservativeResize(new_var_index + 1, new_var_index + 1);
            MNA_Y.row(new_var_index).setZero();
            MNA_Y.col(new_var_index).setZero();

            J.conservativeResize(new_var_index + 1);
            J(new_var_index) = value;

            // 记录该电压源对应的支路电流变量索引（从序号0开始）
            dev.branch_current_index = new_var_index - ((int)ckt.node_list.size() - 1);

            // 如果需要打印/绘制支路电流，索引要减去节点数偏移
            if (dev.printI){
                ckt.print_branch_current_indices.push_back(new_var_index - ((int)ckt.node_list.size() - 1));
            }
            if (dev.plotI){
                ckt.plot_branch_current_indices.push_back(new_var_index - ((int)ckt.node_list.size() - 1));
            }

            // 如果是瞬态等效电压源，更新动态器件映射
            if (!dev.original_device_name.empty()) {
                dynamic_device_current_map[dev.original_device_name] =
                    new_var_index - ((int)ckt.node_list.size() - 1);
            }

            // 如果是电感的等效电压源，保存指向该电感器件的指针，以便后续更新电流
            if (dev.original_device_ptr) {
                dev.original_device_ptr->branch_current_index =
                    new_var_index - ((int)ckt.node_list.size() - 1);
            }

            // KVL 方程
            if (n1 != 0){
                MNA_Y(new_var_index, n1 - 1) = 1;
                MNA_Y(n1 - 1, new_var_index) = 1;
            }
            if (n2 != 0){
                MNA_Y(new_var_index, n2 - 1) = -1;
                MNA_Y(n2 - 1, new_var_index) = -1;
            }
        }

        // 独立电流源
        if (c == 'I'){
            int n1 = dev.nodes[0];
            int n2 = dev.nodes[1];

            // 斜坡：把源值按 alpha 缩放
            double value = alpha * dev.parameters["DC"];

            if (n1 != 0){
                J(n1 - 1) -= value; // 流出节点为负
            }
            if (n2 != 0){
                J(n2 - 1) += value; // 流入节点为正
            }
        }
    }
}



void solver::build_sources_MNA(bool in_tran,double time){
    //清空打印支路电流索引
    ckt.print_branch_current_indices.clear();
    ckt.plot_branch_current_indices.clear();
    for (auto &dev : ckt.sources){
        char c = toupper(dev.name[0]);
        //独立电压源
        if (c == 'V'){
            int n1 = dev.nodes[0];
            int n2 = dev.nodes[1];
            double value = dev.parameters["DC"];
            if(in_tran) {
                //瞬态分析时，处理瞬态电压源
                if (dev.type == "V_SIN"){
                    //正弦波电压源
                    double A = dev.parameters["AMP"];
                    double FREQ = dev.parameters["FREQ"];
                    double PHASE = dev.parameters["PHASE"];
                    double DC = dev.parameters["DC"];
                    value = DC + A * sin(2 * M_PI * FREQ * time + PHASE * M_PI / 180.0);
                    
                    //调试
                    //std::cout << "Time: " << time << "s, Voltage Source " << dev.name << " Value: " << value << " V\n";
                    //system("pause");  //调试
                }
            }
            //先引入支路电流变量，从n1流向n2
            int new_var_index = MNA_Y.rows();
            MNA_Y.conservativeResize(new_var_index + 1, new_var_index + 1);
            MNA_Y.row(new_var_index).setZero();
            MNA_Y.col(new_var_index).setZero();
            J.conservativeResize(new_var_index + 1);
            J(new_var_index) = value;
            //记录改电压源对应的支路电流变量索引，从序号0开始
            dev.branch_current_index = new_var_index- ((int)ckt.node_list.size() - 1);
            //如果需要打印支路电流，存在ckt.print_branch_current_indices,需要减去节点数偏移
            if(dev.printI){
                ckt.print_branch_current_indices.push_back(new_var_index - ((int)ckt.node_list.size() - 1));
            }
            if(dev.plotI){
                ckt.plot_branch_current_indices.push_back(new_var_index - ((int)ckt.node_list.size() - 1));
            }
            // 如果是瞬态等效电压源，更新动态器件映射
            if (!dev.original_device_name.empty()) {
                dynamic_device_current_map[dev.original_device_name] = new_var_index- ((int)ckt.node_list.size() - 1);
            }
            //如果是电感的等效电压源，保存指向该电感器件的指针，以便后续更新电流
            if (dev.original_device_ptr) {
                dev.original_device_ptr->branch_current_index = new_var_index- ((int)ckt.node_list.size() - 1);
            }
            //KVL方程
            if (n1 != 0){
                MNA_Y(new_var_index, n1 - 1) = 1;
                MNA_Y(n1 - 1, new_var_index) = 1;
            }
            if (n2 != 0){
                MNA_Y(new_var_index, n2 - 1) = -1;
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
}
