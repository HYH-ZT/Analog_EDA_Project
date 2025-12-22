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
                R_eq = L / tstep;
                V_eq = 0;
                I_eq = I_prev;
            }
            else{
                // //使用梯形法戴维南等效参数
                // R_eq = 2 * L / tstep;
                // V_eq = -2 * L * I_prev / tstep - (V1_prev - V2_prev);
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

//瞬态分析

void solver::TRAN_solve(){
    //提取电容信息
    ckt.extract_MOS_capacitances();
    //确定需要打印的变量
    parse_print_variables();
    //确定瞬态分析参数
    double tstop = analysis_type.parameters["tstop"];
    double tstep = analysis_type.parameters["tstep"];

    //遍历print_variables

    //先进行直流分析，获得初始条件
    //DC_solve();
    //零初值条件
    node_voltages = Eigen::VectorXd::Zero(ckt.node_map.size() - 1);

    // //Debug展示初始节点电压结果
    // std::cout << "Initial Node Voltages for Transient Analysis:\n";
    // for (const auto& pair : ckt.node_map){
    //     const std::string& node_name = pair.first;
    //     int node_id = pair.second;
    //     if (node_id == 0){
    //         std::cout << "Node " << node_name << " (ID " << node_id << "): 0 V (Ground)\n";
    //     }
    //     else{
    //         std::cout << "Node " << node_name << " (ID " << node_id << "): " << node_voltages[node_id - 1] << " V\n";
    //     }
    // }


    int steps = static_cast<int>(tstop / tstep);
    tran_plot_data.clear();
    for (int step = 0; step <= steps; ++step){
        double time = step * tstep;
        // std::cout << "Transient Analysis Time: " << time << " s\n";
        //构建瞬态分析电路
        build_transient_ckt(tstep);
        
        // //Debug: 输出当前电路的线性器件和源列表
        // std::cout << "Linear Devices:\n";
        // for (const auto& dev : ckt.linear_devices) {
        //     std::cout << "  " << dev.name << " (" << dev.type << ") Nodes: ";
        //     for (const auto& node_name : dev.node_names) {
        //         std::cout << node_name << " ";
        //     }
        //     std::cout << " Parameters: ";
        //     for (const auto& param : dev.parameters) {
        //         std::cout << param.first << "=" << param.second << " ";
        //     }
        //     std::cout << "\n";
        // }
        // std::cout << "Sources:\n";
        // for (const auto& dev : ckt.sources) {
        //     std::cout << "  " << dev.name << " (" << dev.type << ") Nodes: ";
        //     for (const auto& node_name : dev.node_names) {
        //         std::cout << node_name << " ";
        //     }
        //     std::cout << " Parameters: ";
        //     for (const auto& param : dev.parameters) {
        //         std::cout << param.first << "=" << param.second << " ";
        //     }
        //     std::cout << "\n";
        // }
        // std::cout << "\n";

        //直接调用DC求解器
        //求解非线性MNA方程，以上次节点电压为初值
        DC_solve(node_voltages, true,time);

        // //Debug: 输出MNA矩阵和J向量(注意DC_solve中会行交换MNA)
        // std::cout << "MNA_Y at time " << time << " s:\n" << MNA_Y << "\n";
        // std::cout << "J at time " << time << " s:\n" << J << "\n";

        // //Debug:展示节点电压结果
        // std::cout << "Node Voltages at time " << time << " s:\n";
        // for (const auto& pair : ckt.node_map){
        //     const std::string& node_name = pair.first;
        //     int node_id = pair.second;
        //     if (node_id == 0){
        //         std::cout << "Node " << node_name << " (ID " << node_id << "): 0 V (Ground)\n";
        //     }
        //     else{
        //         std::cout << "Node " << node_name << " (ID " << node_id << "): " << node_voltages[node_id - 1] << " V\n";
        //     }
        // }

        //记录需要画图节点此时的电压
        for (auto plot_node_id : ckt.plot_node_ids){
            double v = 0.0;
            if (plot_node_id == 0) v = 0.0;
            else if (plot_node_id - 1 >= 0 && plot_node_id - 1 < node_voltages.size()) v = node_voltages[plot_node_id - 1];
            tran_plot_data[plot_node_id].push_back(std::make_pair(time, v));
            //输出调试信息
            //std::cout << "Plot Data - Time: " << time << " s, Node ID: " << plot_node_id << ", Voltage: " << v << " V\n";
        }
        //记录需要画图的支路电流
        for (auto plot_current_dev_index : ckt.plot_branch_current_indices){
            if (plot_current_dev_index >= 0 && plot_current_dev_index < branch_currents.size()){
                double i = branch_currents[plot_current_dev_index];
                tran_plot_data[-(plot_current_dev_index + 1)].push_back(std::make_pair(time, i));
                //输出调试信息
                //std::cout << "Plot Data - Time: " << time << " s, Branch Index: " << plot_current_dev_index << ", Current: " << i << " A\n";
            }
        }

        // 根据需要打印的变量，存到文件中
        {
            // 输出文件: transient_print.txt
            if (step == 0) {
                std::ofstream hdr("transient_print.txt", std::ios::out);
                hdr << "Time(s)";
                for (int node_id : ckt.print_node_ids){
                    std::string name = "NODE";
                    if (node_id >= 0 && node_id < (int)ckt.node_list.size()) name = ckt.node_list[node_id];
                    hdr << "\tV(" << name << ")";
                }
                //只需要遍历所有sources，按顺序输出支路电流表头
                for (const auto &d : ckt.sources){
                    if (d.printI) hdr << "\tI(" << d.name << ")";
                }
                //关闭
                hdr << "\n";
                hdr.close();
            }

            std::ofstream out("transient_print.txt", std::ios::app);
            out << time;
            for (int node_id : ckt.print_node_ids) {
                double v = 0.0;
                if (node_id == 0) v = 0.0;
                else if (node_id - 1 >= 0 && node_id - 1 < node_voltages.size()) v = node_voltages[node_id - 1];
                out << "\t" << v;
            }
            for (int current_dev_index : ckt.print_branch_current_indices) {
                if (current_dev_index >= 0 && current_dev_index < branch_currents.size()){
                    out << "\t" << branch_currents[current_dev_index];
                }
            }

            out << "\n";
            out.close();
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
    

    // //Debug展示初始节点电压结果
    // std::cout << "Initial Node Voltages for Transient Analysis:\n";
    // for (const auto& pair : ckt.node_map){
    //     const std::string& node_name = pair.first;
    //     int node_id = pair.second;
    //     if (node_id == 0){
    //         std::cout << "Node " << node_name << " (ID " << node_id << "): 0 V (Ground)\n";
    //     }
    //     else{
    //         std::cout << "Node " << node_name << " (ID " << node_id << "): " << node_voltages[node_id - 1] << " V\n";
    //     }
    // }

    //进行前向欧拉瞬态分析
    // double tstop = analysis_type.parameters["tstop"];
    // double tstep = analysis_type.parameters["tstep"];
    int steps = static_cast<int>(tstop / tstep);
    tran_plot_data.clear();
    for (int step = 0; step <= steps; ++step){
        double time = step * tstep;
        std::cout << "Transient Analysis Time: " << time << " s\n";
        //构建瞬态分析电路
        build_transient_ckt(tstep);
        
        // //Debug: 输出当前电路的线性器件和源列表
        // std::cout << "Linear Devices:\n";
        // for (const auto& dev : ckt.linear_devices) {
        //     std::cout << "  " << dev.name << " (" << dev.type << ") Nodes: ";
        //     for (const auto& node_name : dev.node_names) {
        //         std::cout << node_name << " ";
        //     }
        //     std::cout << " Parameters: ";
        //     for (const auto& param : dev.parameters) {
        //         std::cout << param.first << "=" << param.second << " ";
        //     }
        //     std::cout << "\n";
        // }
        // std::cout << "Sources:\n";
        // for (const auto& dev : ckt.sources) {
        //     std::cout << "  " << dev.name << " (" << dev.type << ") Nodes: ";
        //     for (const auto& node_name : dev.node_names) {
        //         std::cout << node_name << " ";
        //     }
        //     std::cout << " Parameters: ";
        //     for (const auto& param : dev.parameters) {
        //         std::cout << param.first << "=" << param.second << " ";
        //     }
        //     std::cout << "\n";
        // }

        //直接调用DC求解器
        //求解非线性MNA方程，以上次节点电压为初值
        DC_solve(node_voltages, true, time);

        // //Debug:展示节点电压结果
        // std::cout << "Node Voltages at time " << time << " s:\n";
        // for (const auto& pair : ckt.node_map){
        //     const std::string& node_name = pair.first;
        //     int node_id = pair.second;
        //     if (node_id == 0){
        //         std::cout << "Node " << node_name << " (ID " << node_id << "): 0 V (Ground)\n";
        //     }
        //     else{
        //         std::cout << "Node " << node_name << " (ID " << node_id << "): " << node_voltages[node_id - 1] << " V\n";
        //     }
        // }

        //记录需要画图节点此时的电压
        for (auto plot_node_id : ckt.plot_node_ids){
            double v = 0.0;
            if (plot_node_id == 0) v = 0.0;
            else if (plot_node_id - 1 >= 0 && plot_node_id - 1 < node_voltages.size()) v = node_voltages[plot_node_id - 1];
            tran_plot_data[plot_node_id].push_back(std::make_pair(time, v));
            //输出调试信息
            //std::cout << "Plot Data - Time: " << time << " s, Node ID: " << plot_node_id << ", Voltage: " << v << " V\n";
        }
        //记录需要画图的支路电流
        for (auto plot_current_dev_index : ckt.plot_branch_current_indices){
            if (plot_current_dev_index >= 0 && plot_current_dev_index < branch_currents.size()){
                double i = branch_currents[plot_current_dev_index];
                tran_plot_data[-(plot_current_dev_index+1)].push_back(std::make_pair(time, i));
                // //输出调试信息
                // std::cout << "save current index" <<plot_current_dev_index+1 << "\n";
                //std::cout << "Plot Data - Time: " << time << " s, Branch Index: " << plot_current_dev_index << ", Current: " << i << " A\n";
            }
        }

        //根据需要打印的变量，存到文件中
        {
            // 输出文件: transient_print.txt
            if (step == 0) {
                std::ofstream hdr("transient_print.txt", std::ios::out);
                hdr << "Time(s)";
                for (int node_id : ckt.print_node_ids) {
                    std::string name = "NODE";
                    if (node_id >= 0 && node_id < (int)ckt.node_list.size()) name = ckt.node_list[node_id];
                    hdr << "\tV(" << name << ")";
                }
                //只需要遍历所有sources，按顺序输出支路电流表头
                for (const auto &d : ckt.sources){
                    if (d.printI) hdr << "\tI(" << d.name << ")";
                }
                //关闭
                hdr << "\n";
                hdr.close();
            }

            std::ofstream out("transient_print.txt", std::ios::app);
            out << time;
            for (int node_id : ckt.print_node_ids) {
                double v = 0.0;
                if (node_id == 0) v = 0.0;
                else if (node_id - 1 >= 0 && node_id - 1 < node_voltages.size()) v = node_voltages[node_id - 1];
                out << "\t" << v;
            }
            for (int current_dev_index : ckt.print_branch_current_indices) {
                if (current_dev_index >= 0 && current_dev_index < branch_currents.size()){
                    out << "\t" << branch_currents[current_dev_index];
                }
            }

            out << "\n";
            out.close();
        }

    }
}


//瞬态分析，初值可调节

void solver::TRAN_solve_with_initial_value(double tstop, double tstep){
    int steps = static_cast<int>(tstop / tstep);
    tran_plot_data.clear();
    for (int step = 0; step <= steps; ++step){
        double time = step * tstep;
        // std::cout << "Transient Analysis Time: " << time << " s\n";

        //直接调用DC求解器
        //求解非线性MNA方程，以上次节点电压为初值
        DC_solve(node_voltages, true, time);

        //记录需要画图节点此时的电压
        for (auto plot_node_id : ckt.plot_node_ids){
            double v = 0.0;
            if (plot_node_id == 0) v = 0.0;
            else if (plot_node_id - 1 >= 0 && plot_node_id - 1 < node_voltages.size()) v = node_voltages[plot_node_id - 1];
            tran_plot_data[plot_node_id].push_back(std::make_pair(time, v));
            //输出调试信息
            //std::cout << "Plot Data - Time: " << time << " s, Node ID: " << plot_node_id << ", Voltage: " << v << " V\n";
        }
        //记录需要画图的支路电流
        for (auto plot_current_dev_index : ckt.plot_branch_current_indices){
            if(plot_current_dev_index >=0 && plot_current_dev_index < ckt.sources.size()){
                double i = branch_currents[plot_current_dev_index];
                tran_plot_data[-(plot_current_dev_index+1)].push_back(std::make_pair(time, i));
                // //输出调试信息
                // std::cout << "save current index" <<plot_current_dev_index+1 << "\n";
                //std::cout << "Plot Data - Time: " << time << " s, Branch Index: " << plot_current_dev_index << ", Current: " << i << " A\n";
            }
        }
    }
}

// HB分析
