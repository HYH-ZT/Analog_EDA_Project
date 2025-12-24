#include "solver.hpp"
#include "solver_internal.hpp"
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <string>

void solver::DC_solve() {
    const int maxNewtonIter = 500;
    const double tol = 1e-9;
    //确定需要打印的节点电压和支路电流
    parse_print_variables();

    // 1. 只构建一次线性矩阵
    build_linear_MNA();

    int nodeCount = (int)ckt.node_list.size() - 1;
    node_voltages = Eigen::VectorXd::Zero(nodeCount);

    int iter;
    for (iter = 0; iter < maxNewtonIter; iter++) {

        // 2. 每次迭代重新构造 MNA
        MNA_Y = liner_Y;
        J = Eigen::VectorXd::Zero(MNA_Y.rows());

        // 3. MOS stamp
        build_nonlinear_MNA();

        // 4. 电源 stamp
        build_sources_MNA();

        // 5. 求解线性方程
        //保存旧节点电压用于收敛性检查
        Eigen::VectorXd old_node_voltages = node_voltages;
        solve_linear_MNA();

        // 6. 检查收敛性
        double max_diff = (node_voltages - old_node_voltages).cwiseAbs().maxCoeff();
        if (max_diff < tol) {
            std::cout << "Converged after " << iter + 1 << " iterations.\n";
            break;
        }
    }

    if(iter == maxNewtonIter){ 
        std::cout << "Warning: DC did NOT converge.\n";
    }
}



void solver::DC_solve_ramp() {
    const int    maxNewtonIter = 500;
    const double tol           = 1e-9;

    // 自适应步长参数（可按电路规模调整）
    const double alphaMinStep  = 1e-4;   // 最小步长：再小就认为失败
    const double alphaMaxStep  = 0.2;    // 最大步长：避免一下子跨太大
    double       alphaStep     = 0.05;   // 初始步长（建议 0.02~0.1）

    // 收敛快则加速（阈值/倍率可调）
    const int    fastIterThres = 8;      // 小于该迭代次数视为“很快”
    const double growFactor    = 1.5;    // 步长放大倍数
    const double shrinkFactor  = 0.5;    // 步长缩小倍数

    parse_print_variables();
    build_linear_MNA();

    int nodeCount = (int)ckt.node_list.size() - 1;
    node_voltages = Eigen::VectorXd::Zero(nodeCount);

    // 从 alpha=0 开始往 1 推进
    double alpha = 0.0;

    // 为了能在“某一步失败”时回退，保存上一次成功解
    Eigen::VectorXd last_good_solution = node_voltages;

    while (alpha < 1.0 - 1e-15) {
        // 这一步的目标 alpha
        double next_alpha = alpha + alphaStep;
        if (next_alpha > 1.0) next_alpha = 1.0;

        // 用上一次成功解作为本步初值
        node_voltages = last_good_solution;

        bool converged = false;
        int  iter      = 0;

        for (iter = 0; iter < maxNewtonIter; ++iter) {
            MNA_Y = liner_Y;
            J     = Eigen::VectorXd::Zero(MNA_Y.rows());

            build_nonlinear_MNA();

            // 源按 next_alpha 加载
            build_sources_MNA_ramp(next_alpha);

            Eigen::VectorXd old_node_voltages = node_voltages;
            solve_linear_MNA();

            double max_diff = (node_voltages - old_node_voltages).cwiseAbs().maxCoeff();
            if (max_diff < tol) {
                converged = true;
                break;
            }
        }

        if (converged) {
            // 接受本步
            alpha = next_alpha;
            last_good_solution = node_voltages;

            // 如果收敛很快，尝试增大步长加速
            if (iter < fastIterThres) {
                alphaStep = std::min(alphaStep * growFactor, alphaMaxStep);
            }

        } else {
            // 本步失败：缩小步长并重试
            alphaStep *= shrinkFactor;

            // 回退到上一次成功解（确保状态一致）
            node_voltages = last_good_solution;

            if (alphaStep < alphaMinStep) {
                std::cout << "Warning: DC did NOT converge during ramp stepping.\n";
                std::cout << "Stopped at alpha=" << alpha
                          << " (next_alpha=" << next_alpha << "), alphaStep=" << alphaStep << "\n";
                break;
            }
        }
    }

    if (alpha >= 1.0 - 1e-15) {
        std::cout << "Ramp DC converged to alpha=1.\n";
    }
}


// 使用节点名和电压值的映射来设置初值
//可能删除的版本
void solver::DC_solve(const std::map<std::string, double>& node_voltage_map, bool in_tran) {
    const int maxNewtonIter = 500;
    const double tol = 1e-9;

    // 1. 只构建一次线性矩阵
    build_linear_MNA(in_tran);

    int nodeCount = (int)ckt.node_list.size() - 1;
    node_voltages = Eigen::VectorXd::Zero(nodeCount);

    // 根据输入的节点名和电压值设置初值
    for (const auto& entry : node_voltage_map) {
        const std::string& node_name = entry.first;
        double voltage = entry.second;
        
        // 检查节点是否存在
        auto it = ckt.node_map.find(node_name);
        if (it != ckt.node_map.end()) {
            int node_id = it->second;
            if (node_id > 0 && node_id <= nodeCount) {
                node_voltages(node_id - 1) = voltage;
                std::cout << "Set initial voltage for node " << node_name 
                         << " (ID " << node_id << "): " << voltage << " V\n";
            }
        } else {
            std::cout << "Warning: Node " << node_name << " not found in circuit\n";
        }
    }
    int iter;
    for (iter = 0; iter < maxNewtonIter; iter++) {

        // 2. 每次迭代重新构造 MNA
        MNA_Y = liner_Y;
        J = Eigen::VectorXd::Zero(MNA_Y.rows());

        // 3. MOS stamp
        build_nonlinear_MNA();

        // 4. 电源 stamp
        build_sources_MNA();

        // 5. 求解线性方程
        //保存旧节点电压用于收敛性检查
        Eigen::VectorXd old_node_voltages = node_voltages;
        solve_linear_MNA();

        // 或许可以用（
        // Eigen::VectorXd delta_node_voltages = node_voltages - old_node_voltages;
        // node_voltages = old_node_voltages + 0.8 * delta_node_voltages;

        // // Debug: 输出当前迭代的节点电压
        // std::cout << "Node Voltages:\n" << node_voltages << std::endl << std::endl; 

        // 6. 检查收敛性
        double max_diff = (node_voltages - old_node_voltages).cwiseAbs().maxCoeff();
        // std::cout << "Max voltage change: " << max_diff << "\n";
        if (max_diff < tol) {
            std::cout << "Converged after " << iter + 1 << " iterations.\n";
            break;
        }
    }
    if(iter == maxNewtonIter){ 
        std::cout << "Warning: DC did NOT converge.\n";
    }

    //展示节点电压结果
    std::cout << "DC Analysis Node Voltages:\n";
    for (const auto& pair : ckt.node_map){
        const std::string& node_name = pair.first;
        int node_id = pair.second;
        if (node_id == 0){
            std::cout << "Node " << node_name << " (ID " << node_id << "): 0 V (Ground)\n";
        }
        else{
            std::cout << "Node " << node_name << " (ID " << node_id << "): " << node_voltages[node_id - 1] << " V\n";
        }
    }
}


//根据给定的初始节点电压进行直流求解

void solver::DC_solve(const Eigen::VectorXd& initial_node_voltages, bool in_tran, double time) {
    const int maxNewtonIter = 1000;
    const double tol = 1e-9;

    // 1. 只构建一次线性矩阵
    build_linear_MNA(in_tran);

    int nodeCount = (int)ckt.node_list.size() - 1;
    node_voltages = initial_node_voltages;
    
    int iter;
    for (iter = 0; iter < maxNewtonIter; iter++) {

        // 2. 每次迭代重新构造 MNA
        MNA_Y = liner_Y;
        J = Eigen::VectorXd::Zero(MNA_Y.rows());

        // 3. MOS stamp
        build_nonlinear_MNA();

        // 4. 电源 stamp
        build_sources_MNA(in_tran,time);

        // 5. 求解线性方程
        //保存旧节点电压用于收敛性检查
        Eigen::VectorXd old_node_voltages = node_voltages;
        solve_linear_MNA();

        // 6. 检查收敛性
        int original_node_count = std::min((int)old_node_voltages.size(), (int)node_voltages.size());
        double max_diff = 0.0;
        if (original_node_count > 0) {
            max_diff = (node_voltages.head(original_node_count) - old_node_voltages.head(original_node_count)).cwiseAbs().maxCoeff();
        }

        if (max_diff < tol) {
            break;
        }
    }

    if(iter == maxNewtonIter){ 
        std::cout << "Warning: DC did NOT converge.\n";
    }
}

void solver::print_dc_results(){
    std::ofstream out("dc_print_results.txt");
    out << "DC";
    for (int node_id : ckt.print_node_ids) {
        std::string name = "NODE";
        if (node_id >= 0 && node_id < (int)ckt.node_list.size()) {
            name = ckt.node_list[node_id];
        }
        out << "\tV(" << name << ")";
    }
    for (const auto &d : ckt.sources){
        if (d.printI) out << "\tI(" << d.name << ")";
    }
    out << "\n";

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
