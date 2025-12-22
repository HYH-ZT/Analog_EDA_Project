#include "solver.hpp"
#include "solver_internal.hpp"
#include <chrono>
#include <iostream>
#include <map>
#include <string>

void solver::DC_solve() {
    const int maxNewtonIter = 500;
    const double tol = 1e-9;

    // //debug：打印所有器件信息
    // std::cout << "Circuit Devices:\n";
    // for (const auto& dev : ckt.linear_devices) {
    //     std::cout << "Linear Device: " << dev.name << ", Type: " << dev.type << ", Nodes: ";
    //     for (const auto& n : dev.nodes) std::cout << n << " ";
    //     std::cout << "\n";
    // }
    // for (const auto& dev : ckt.nonlinear_devices) {
    //     std::cout << "Nonlinear Device: " << dev.name << ", Type: " << dev.type << ", Nodes: ";
    //     for (const auto& n : dev.nodes) std::cout << n << " ";
    //     std::cout << "\n";
    // }
    // for (const auto& dev : ckt.sources) {
    //     std::cout << "Source: " << dev.name << ", Type: " << dev.type << ", Nodes: ";
    //     for (const auto& n : dev.nodes) std::cout << n << " ";
    //     std::cout << "\n";
    // }

    //确定需要打印的节点电压和支路电流
    parse_print_variables();

    // 1. 只构建一次线性矩阵
    build_linear_MNA();

    // //debug：打印线性MNA矩阵
    // std::cout << "Linear MNA Matrix (liner_Y):\n" << liner_Y << "\n";

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

        // // Debug: 输出当前迭代的 MNA 矩阵和 J 向量
        // std::cout << "Iteration " << iter + 1 << ":\n";
        // //输出节点序号与名称对照
        // for (const auto& pair : ckt.node_map) {
        //     std::cout << "Node ID " << pair.second << ": " << pair.first << "\n";
        // }
        // std::cout << "MNA_Y:\n" << MNA_Y << "\n";
        // std::cout << "J:\n" << J << "\n";

        // 5. 求解线性方程
        //保存旧节点电压用于收敛性检查
        Eigen::VectorXd old_node_voltages = node_voltages;
        solve_linear_MNA();

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

    //std::cout << "DC did NOT converge.\n";

        // //Debug:根据需要打印的变量，输出结果
        // //根据ckt.print_node_ids,ckt.print_branch_current_indices，输出结果
        // for (int node_id : ckt.print_node_ids) {
        //     double voltage = (node_id == 0) ? 0.0 : node_voltages[node_id - 1];
        //     std::cout << "  Voltage at node ID " << node_id << ": " << voltage << " V\n";
        // }
        // for (int branch_index : ckt.print_branch_current_indices) {
        //     if (branch_index >= 0 && branch_index < branch_currents.size()) {
        //         double current = branch_currents[branch_index];
        //         std::cout << "  Current through branch " << branch_index << ": " << current << " A\n";
        //     }
        // }

    // //Debug:展示所有节点电压结果
    // std::cout << "\n";
    // std::cout << "DC Analysis Node Voltages:\n";
    // print_node_voltages();
}



void solver::DC_solve_ramp() {
    const int    maxNewtonIter = 500;
    const double tol           = 1e-9;

    // 自适应步长参数（可按作业/电路规模调整）
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
        // 这一步的目标 alpha（不要超过 1）
        double next_alpha = alpha + alphaStep;
        if (next_alpha > 1.0) next_alpha = 1.0;

        // 用上一次成功解作为本步初值（非常关键）
        node_voltages = last_good_solution;

        bool converged = false;
        int  iter      = 0;

        for (iter = 0; iter < maxNewtonIter; ++iter) {
            MNA_Y = liner_Y;
            J     = Eigen::VectorXd::Zero(MNA_Y.rows());

            build_nonlinear_MNA();

            // 源按 next_alpha 加载（本步目标）
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

            // （可选）打印进度
            // std::cout << "alpha=" << alpha << " converged in " << (iter+1) << " iters, step=" << alphaStep << "\n";
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

            // （可选）打印失败信息
            // std::cout << "alpha step failed at next_alpha=" << next_alpha << ", shrinking step to " << alphaStep << "\n";
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

        // // Debug: 输出当前迭代的 MNA 矩阵和 J 向量
        // std::cout << "Iteration " << iter + 1 << ":\n";
        // std::cout << "MNA_Y:\n" << MNA_Y << "\n";
        // std::cout << "J:\n" << J << "\n";

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

    //std::cout << "DC did NOT converge.\n";

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

        // // // Debug: 输出当前迭代的 MNA 矩阵和 J 向量

        // std::cout << "liner_Y:\n" << liner_Y << "\n";
        // std::cout << "J:\n" << J << "\n";

    int nodeCount = (int)ckt.node_list.size() - 1;
    node_voltages = initial_node_voltages;
    
    int iter;
    for (iter = 0; iter < maxNewtonIter; iter++) {

        // 2. 每次迭代重新构造 MNA
        MNA_Y = liner_Y;
        J = Eigen::VectorXd::Zero(MNA_Y.rows());

        // //Debug: 输出当前MNA
        // std::cout << "Iteration " << iter + 1 << ":\n";
        // std::cout << "MNA_Y before nonlinear stamp:\n" << MNA_Y << "\n";
        // std::cout << "J before nonlinear stamp:\n" << J << "\n";

        // 3. MOS stamp
        build_nonlinear_MNA();

        // // Debug: 输出当前迭代的 MNA 矩阵和 J 向量
        // std::cout << "Iteration " << iter + 1 << ":\n";
        // std::cout << "MNA_Y:\n" << MNA_Y << "\n";
        // std::cout << "J:\n" << J << "\n";

        // 4. 电源 stamp
        build_sources_MNA(in_tran,time);

        // // Debug: 输出当前迭代的 MNA 矩阵和 J 向量
        // std::cout << "Source " << iter + 1 << ":\n";
        // std::cout << "Before_solve_MNA_Y:\n" << MNA_Y << "\n";
        // std::cout << "J:\n" << J << "\n";


        // 5. 求解线性方程
        //保存旧节点电压用于收敛性检查
        Eigen::VectorXd old_node_voltages = node_voltages;
        solve_linear_MNA();

        // //Debug: 输出当前迭代的解
        // std::cout << "after_solve Node Voltages:\n" << node_voltages << std::endl;
        // std::cout << "currents" << branch_currents << std::endl;


        // 6. 检查收敛性
        // double max_diff = (node_voltages - old_node_voltages).cwiseAbs().maxCoeff();
        int original_node_count = std::min((int)old_node_voltages.size(), (int)node_voltages.size());
        double max_diff = 0.0;
        if (original_node_count > 0) {
            max_diff = (node_voltages.head(original_node_count) - old_node_voltages.head(original_node_count)).cwiseAbs().maxCoeff();
        }
        // std::cout << "Max voltage change: " << max_diff << "\n";
        if (max_diff < tol) {
            // std::cout << "Converged after " << iter + 1 << " iterations.\n";
            break;
        }
    }

    if(iter == maxNewtonIter){ 
        std::cout << "Warning: DC did NOT converge.\n";
    }

    // //Debug展示节点电压结果
    // std::cout << "DC Analysis Node Voltages:\n";
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
}


//基于梯形欧拉法构建瞬态分析电路，修改ckt，而不是MNA矩阵
