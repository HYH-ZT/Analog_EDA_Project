#include "solver.hpp"
#include "solver_internal.hpp"
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

void solver::build_transient_ckt(double tstep){
    // 清除之前添加的瞬态等效器件
    auto it_linear = ckt.linear_devices.begin();
    while (it_linear != ckt.linear_devices.end()) {
        if (!it_linear->original_device_name.empty()) { // 使用属性判断而不是字符串查找
            it_linear = ckt.linear_devices.erase(it_linear);
        } else {
            ++it_linear;
        }
    }
    
    auto it_sources = ckt.sources.begin();
    while (it_sources != ckt.sources.end()) {
        if (!it_sources->original_device_name.empty()) { // 使用属性判断而不是字符串查找
            it_sources = ckt.sources.erase(it_sources);
        } else {
            ++it_sources;
        }
    }
    
    // 处理线性器件中的电容和电感
    // 先复制一份原始器件列表，避免在遍历过程中修改集合导致的问题
    std::vector<device> original_linear_devices;
    for (const auto &dev : ckt.linear_devices) {
        if (dev.original_device_name.empty()) { // 只处理原始器件，不处理等效器件
            original_linear_devices.push_back(dev);
        }
    }
    
    for (const auto &dev : original_linear_devices) {
        // 安全检查
        if (dev.name.empty() || dev.nodes.size() < 2 || dev.node_names.size() < 2) {
            continue;
        }
        
        char c = toupper(dev.name[0]);
        // 获取上一时刻的节点电压
        auto V_prev = [&](int node) -> double {
            if (node == 0) return 0.0;
            int idx = node - 1;
            if (idx < 0 || idx >= (int)node_voltages.size()) return 0.0;
            return node_voltages[idx];
        };
        if (c == 'C') {
            // 电容处理：基于梯形积分法，使用电压源串联电阻的形式
            // i_C(t) = C * dv/dt ≈ (C/tstep) * (v_n+1 - v_n)
            // 等效为电压源串联电阻：V_eq 串联 R_eq
            // 其中 R_eq = tstep/C, V_eq = -v_n
            
            int n1 = dev.nodes[0];
            int n2 = dev.nodes[1];
            if (dev.parameters.find("value") == dev.parameters.end()) {
                continue; // 跳过没有value参数的器件
            }
            double C = dev.parameters.at("value");
            if (C <= 0) {
                continue; // 跳过无效的电容值
            }

            //获取上一时刻支路电流
            //直接从映射中获取该电容对应的支路电流
            double I_prev = 0.0;
            auto it = dynamic_device_current_map.find(dev.name);
            if (it != dynamic_device_current_map.end() && 
                it->second >= 0 && it->second < branch_currents.size()) {
                I_prev = branch_currents[it->second];
            } else {
                // 如果找不到之前的支路电流，使用初始直流分析的值
                I_prev = 0.0;
            }

            double V1_prev = V_prev(n1);
            double V2_prev = V_prev(n2);
            
            // 梯形积分法：v_C(t+dt) = v_C(t) + (dt/C) * (i_C(t+dt) + i_C(t))/2
            // 重新整理：i_C(t+dt) = (2C/dt) * (v_C(t+dt) - v_C(t)) - i_C(t)
            // 等效电路：v_C(t+dt) = R_eq * i_C(t+dt) + V_eq
            // 其中：R_eq = dt/(2C), V_eq = v_C(t) - R_eq * i_C(t)
            //前向欧拉法等效参数
            double R_eq, V_eq;
            if(transient_method == TransientMethod::FORWARD_EULER){
                R_eq = 0.0;
                V_eq = (V1_prev - V2_prev) + (tstep / C) * I_prev;
            }
            //后退欧拉法等效参数
            else if(transient_method == TransientMethod::BACKWARD_EULER){
                R_eq = tstep / C;
                V_eq = (V1_prev - V2_prev);
            }
            else{
                //默认使用梯形法
                R_eq = tstep / (2.0 * C);
                V_eq = (V1_prev - V2_prev) + R_eq * I_prev;
            }
            
            //如果R_eq为0，表示纯电压源，直接添加电压源
            if (R_eq == 0.0){
                device equiv_voltage;
                equiv_voltage.name = "V" + dev.name + "_transient_V";
                equiv_voltage.type = "V";
                equiv_voltage.nodes = {n1, n2};
                equiv_voltage.node_names = {dev.node_names[0], dev.node_names[1]};
                equiv_voltage.parameters["DC"] = V_eq;
                equiv_voltage.original_device_name = dev.name; // 电容等效器件标识
                ckt.sources.push_back(equiv_voltage);
                continue;
            }
            else if (R_eq < 0.0){
                //无效电阻，跳过
                continue;
            }
            else{
                // 创建中间节点用于串联
                std::string mid_node_name = dev.name + "_transient_mid";
                int mid_node_id = ckt.getNodeID(mid_node_name);
                
                // 创建等效电阻 (从n1到中间节点)
                device equiv_resistor;
                equiv_resistor.name = "R" + dev.name + "_transient_R";
                equiv_resistor.type = "R";
                equiv_resistor.nodes = {n1, mid_node_id};
                equiv_resistor.node_names = {dev.node_names[0], mid_node_name};
                equiv_resistor.parameters["value"] = R_eq;
                equiv_resistor.original_device_name = dev.name; // 电容等效器件标识
                ckt.linear_devices.push_back(equiv_resistor);
                
                // 创建等效电压源 (从中间节点到n2)
                device equiv_voltage;
                equiv_voltage.name = "V" + dev.name + "_transient_V";
                equiv_voltage.type = "V";
                equiv_voltage.nodes = {mid_node_id, n2};
                equiv_voltage.node_names = {mid_node_name, dev.node_names[1]};
                equiv_voltage.parameters["DC"] = V_eq;
                equiv_voltage.original_device_name = dev.name; // 电容等效器件标识
                ckt.sources.push_back(equiv_voltage);
            }
        }
        
        if (c == 'L') {
            // 电感处理：基于梯形积分法，使用电压源串联电阻的形式
            // v_L(t) = L * di/dt ≈ (L/tstep) * (i_n+1 - i_n)
            // 等效为电压源串联电阻：V_eq 串联 R_eq
            // 其中 R_eq = L/tstep, V_eq = -L*i_n/tstep
            
            int n1 = dev.nodes[0];
            int n2 = dev.nodes[1];
            if (dev.parameters.find("value") == dev.parameters.end()) {
                continue; // 跳过没有value参数的器件
            }
            double L = dev.parameters.at("value");
            if (L <= 0) {
                continue; // 跳过无效的电感值
            }

            
            // 获取上一时刻的电感电流
            // 直接从映射中获取该电感对应的支路电流
            double I_prev = 0.0;
            //这里用map查找，可能会慢一些，改进方法是直接保存动态器件的支路电流索引
            auto it = dynamic_device_current_map.find(dev.name);
            if (it != dynamic_device_current_map.end() && 
                it->second >= 0 && it->second < branch_currents.size()) {
                I_prev = branch_currents[it->second];
            }
            // 如果找不到，使用初始值0

            double V1_prev = V_prev(n1);
            double V2_prev = V_prev(n2);

            double R_eq, V_eq,I_eq;
            if(transient_method == TransientMethod::FORWARD_EULER){
                R_eq = -1; //表示无穷大
                V_eq = 0;
                I_eq = I_prev + (tstep / L) * (V1_prev - V2_prev);
            }
            else if(transient_method == TransientMethod::BACKWARD_EULER){
                // //使用梯形法戴维南等效参数
                // R_eq = L / tstep;
                // V_eq = - L * I_prev / tstep;    
                // I_eq = 0;
                //使用梯形法诺顿等效参数
                R_eq = L / tstep;
                V_eq = 0;
                I_eq = I_prev;
            }
            else{
                // //使用梯形法戴维南等效参数
                // R_eq = 2 * L / tstep;
                // V_eq = -2 * L * I_prev / tstep - (V1_prev - V2_prev);
                // I_eq = 0;
                //使用梯形法诺顿等效参数
                R_eq = 2 * L / tstep;
                V_eq = 0;
                I_eq = I_prev + (tstep / (2.0 * L)) * (V1_prev - V2_prev);
            }
            //如果R_eq为-1，表示纯电流源，直接添加电流源
            // 创建中间节点用于串联电压源，用于记录电流
            std::string mid_node_name = dev.name + "_transient_mid";
            int mid_node_id = ckt.getNodeID(mid_node_name);
            
            // 创建等效电流源和与其并联电阻 (从n1到中间节点)
            if(I_eq != 0.0){
            device equiv_current;
            equiv_current.name = "I" + dev.name + "_transient_I";
            equiv_current.type = "I";
            equiv_current.nodes = {n1, mid_node_id};
            equiv_current.node_names = {dev.node_names[0], mid_node_name};
            equiv_current.parameters["DC"] = I_eq;
            equiv_current.original_device_name = dev.name; // 电感等效器件标识
            ckt.sources.push_back(equiv_current);
            }

            // 创建等效电阻 (从n1到中间节点)
            if(R_eq > 0.0){
                device equiv_resistor;
                equiv_resistor.name = "R" + dev.name + "_transient_R";
                equiv_resistor.type = "R";
                equiv_resistor.nodes = {n1, mid_node_id};
                equiv_resistor.node_names = {dev.node_names[0], mid_node_name};
                equiv_resistor.parameters["value"] = R_eq;
                equiv_resistor.original_device_name = dev.name; // 电感等效器件标识
                ckt.linear_devices.push_back(equiv_resistor);
            }
            // 创建等效电压源 (从中间节点到n2)
            device equiv_voltage;
            equiv_voltage.name = "V" + dev.name + "_transient_V";
            equiv_voltage.type = "V";
            equiv_voltage.nodes = {mid_node_id, n2};
            equiv_voltage.node_names = {mid_node_name, dev.node_names[1]};
            equiv_voltage.parameters["DC"] = V_eq;
            equiv_voltage.original_device_name = dev.name; // 电感等效器件标识
            ckt.sources.push_back(equiv_voltage);
        }
    }
}
    
//瞬态分析，参数可调节

void solver::TRAN_solve(double tstop, double tstep,int use_initial_dc){
    //提取电容信息
    ckt.extract_MOS_capacitances();
    //确定需要打印的变量
    parse_print_variables();

    //先进行直流分析，获得初始条件
    if(use_initial_dc == 2){
        DC_solve();
    }
    else{
        //零初值条件
        node_voltages = Eigen::VectorXd::Zero(ckt.node_map.size() - 1);
    }
    


    //进行前向欧拉瞬态分析
    // double tstop = analysis_type.parameters["tstop"];
    // double tstep = analysis_type.parameters["tstep"];
    int steps = static_cast<int>(tstop / tstep);
    tran_plot_data.clear();
    tran_print_data.clear();
    for (int step = 0; step <= steps; ++step){
        double time = step * tstep;
        //构建瞬态分析电路
        build_transient_ckt(tstep);
        

        //直接调用DC求解器
        //求解非线性MNA方程，以上次节点电压为初值
        DC_solve(node_voltages, true, time);


        //记录需要画图节点此时的电压
        for (auto plot_node_id : ckt.plot_node_ids){
            double v = 0.0;
            if (plot_node_id == 0) v = 0.0;
            else if (plot_node_id - 1 >= 0 && plot_node_id - 1 < node_voltages.size()) v = node_voltages[plot_node_id - 1];
            tran_plot_data[plot_node_id].push_back(std::make_pair(time, v));
        }
        //记录需要画图的支路电流
        for (auto plot_current_dev_index : ckt.plot_branch_current_indices){
            if (plot_current_dev_index >= 0 && plot_current_dev_index < branch_currents.size()){
                double i = branch_currents[plot_current_dev_index];
                tran_plot_data[-(plot_current_dev_index+1)].push_back(std::make_pair(time, i));
            }
        }

        //记录需要打印节点此时的电压
        for (auto print_node_id : ckt.print_node_ids){
            double v = 0.0;
            if (print_node_id == 0) v = 0.0;
            else if (print_node_id - 1 >= 0 && print_node_id - 1 < node_voltages.size()) v = node_voltages[print_node_id - 1];
            tran_print_data[print_node_id].push_back(std::make_pair(time, v));
        }
        //记录需要打印的支路电流
        for (auto print_current_dev_index : ckt.print_branch_current_indices){
            if (print_current_dev_index >= 0 && print_current_dev_index < branch_currents.size()){
                double i = branch_currents[print_current_dev_index];
                tran_print_data[-(print_current_dev_index+1)].push_back(std::make_pair(time, i));
            }
        }
        



    }
}

//打印瞬态结果
void solver::print_tran_results(){
    //输出tran_print_data到文件
    std::ofstream out("transient_print_data.txt", std::ios::out);
    //输出表头
    out << "Time(s)";
    for (const auto& pair : tran_print_data){
        int id = pair.first;
        if (id >= 0){
            std::string name = "NODE";
            if (id >= 0 && id < (int)ckt.node_list.size()) name = ckt.node_list[id];
            out << "\tV(" << name << ")";
        }
        else{
            int branch_index = -(id + 1);
            std::string dev_name = "BRANCH";
            if (branch_index >= 0 && branch_index < (int)ckt.sources.size()){
                dev_name = ckt.sources[branch_index].name;
            }
            out << "\tI(" << dev_name << ")";
        }
    }
    out << "\n";

    //假设所有数据长度相同，遍历时间点
    if (!tran_print_data.empty()){
        size_t data_size = tran_print_data.begin()->second.size();
        for (size_t i = 0; i < data_size; ++i){
            //输出时间
            out << tran_print_data.begin()->second[i].first;
            //输出各变量数据
            for (const auto& pair : tran_print_data){
                out << "\t" << pair.second[i].second;
            }
            out << "\n";
        }
    }
    out.close();
}


Eigen::MatrixXd solver::TRAN_solve_return(double tstop, double tstep,int use_initial_dc){

    //初始化返回矩阵
    int ori_node_size = ckt.node_map.size() - 1; //节点电压数量
    // int steps = static_cast<int>(tstop / tstep);
    int steps = static_cast<int>(std::round(tstop / tstep));
    Eigen::MatrixXd results = Eigen::MatrixXd::Zero(base_size, steps + 1);
    // std::cout << "Transient Analysis Results Matrix Size: " << results.rows() << " x " << results.cols() << "\n";
    //先进行直流分析，获得初始条件
    if(use_initial_dc == 2){
        DC_solve();
    }
    else{
        //零初值条件
        node_voltages = Eigen::VectorXd::Zero(ckt.node_map.size() - 1);
    }
    // std::cout << "Initial Node Voltages for Transient Analysis:\n" << node_voltages << "\n";
    for (int step = 0; step <= steps; ++step){
        double time = step * tstep;
        // std::cout << "Transient Analysis Time: " << time << " s\n";
        //构建瞬态分析电路
        build_transient_ckt(tstep);
        // std::cout << "Transient Circuit built for time " << time << " s\n";
        //直接调用DC求解器
        //求解非线性MNA方程，以上次节点电压为初值
        DC_solve(node_voltages, true, time);
        // clear_transient_equiv_devices();
        // std::cout << "Node Voltages at time " << time << " s:\n" << node_voltages << "\n";
        //记录当前时刻的节点电压到返回矩阵
        // std::cout << "node_map: " << ckt.node_map.size() << " \n";
        for (int i = 0; i < ori_node_size; ++i){
            // std::cout << i << ": " << node_voltages[i] << "\n";
            results(i, step) = node_voltages(i);
        }
        // std::cout << "Branch Currents at time " << time << " s:\n";
        // //电流节点也存入(还没有找到好方法存电流，因为瞬态的电路多了很多电压源，电流节点索引不对)
        // for (int i = ckt.node_map.size() - 1; i < base_size; ++i){
        //     results(i, step) = branch_currents[i - (ckt.node_map.size() - 1)];
        // }
        // std::cout << "complete" << "\n";
    }   

    // //Debug: 输出返回矩阵
    // std::cout << "Transient Analysis Results Matrix:\n" << results << "\n";
    return results;
}