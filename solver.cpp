#include "solver.hpp"
#include <Eigen/Dense>
#include <iostream>
#include "circuit.hpp"
#include <map>
#include <string>
#include <vector>
#include <fstream>
#include <chrono>
#include <cmath>

#include <cstdlib>

#define OOB_ABORT(msg) do { \
  std::cerr << "\n[OOB_ABORT] " << msg << "\n"; \
  std::abort(); \
} while(0)

#define CHECK_IDX(idx, n, where) do { \
  if ((idx) < 0 || (idx) >= (n)) { \
    std::cerr << "\n[OOB] " << where \
              << " idx=" << (idx) << " size=" << (n) << "\n"; \
    std::abort(); \
  } \
} while(0)

#ifndef M_PI
constexpr double M_PI = 3.141592653589793238462643383279502884;
#endif



solver::solver(circuit& ckt_, analysis& analysis_, 
               LinearSolverMethod lsm, 
               TransientMethod tm,
               SteadyStateMethod ssm)
    : ckt(ckt_), analysis_type(analysis_),
      linear_solver_method(lsm), transient_method(tm), steady_state_method(ssm) {
    liner_Y.resize(ckt.node_map.size() - 1, ckt.node_map.size() - 1); //不包含地节点
    liner_Y.setZero();
    J.resize(ckt.node_map.size() - 1);
    J.setZero();
    
    // // 设置默认求解方法
    // linear_solver_method = LinearSolverMethod::LU_DECOMPOSITION;
    // transient_method = TransientMethod::TRAPEZOIDAL;
    // steady_state_method = SteadyStateMethod::SHOOTING;
}

void solver::print_branch_current(int branch_index){
    if (branch_index >= 0 && branch_index < branch_currents.size()){
        std::cout << "I(" << ckt.sources[branch_index].name << ") " << branch_currents[branch_index] << " A\n";
    }
    else{
        std::cout << "Error: Branch current index " << branch_index << " is out of range.\n";
    }
}

void solver::set_initial_node_voltages(std::string node_name, double voltage) {
    int node_id = ckt.getNodeID(node_name);
    if (node_id != -1) {
        if (node_id == 0) {
            std::cout << "Warning: Cannot set voltage for ground node.\n";
        } else {
            node_voltages[node_id - 1] = voltage;
        }
    } else {
        std::cout << "Error: Node " << node_name << " not found in the circuit.\n";
    }
}

void solver::print_node_voltages() {
    std::cout << "Node Voltages:\n";
    for (const auto& pair : ckt.node_map) {
        const std::string& node_name = pair.first;
        int node_id = pair.second;
        if (node_id == 0) {
            std::cout << "Node " << node_name << " (ID " << node_id << "): 0 V (Ground)\n";
        } else {
            std::cout << "Node " << node_name << " (ID " << node_id << "): " << node_voltages[node_id - 1] << " V\n";
        }
    }
}

void solver::print_node_voltage(const int node_id) {
    std::string node_name = ckt.node_list[node_id];
    if (node_id == 0) {
        std::cout << "Node " << node_name << " (ID " << node_id << "): 0 V (Ground)\n";
    } else if (node_id > 0 && node_id < ckt.node_map.size()) {
        std::cout << "Node " << node_name << " (ID " << node_id << "): " << node_voltages[node_id - 1] << " V\n";
    } else {
        std::cout << "Error: Node ID " << node_id << " is out of range.\n";
    }
}

void solver::print_node_voltage(const std::string& node_name) {
    int node_id = ckt.getNodeID(node_name);
    if (node_id == -1) {
        std::cout << "Error: Node " << node_name << " not found in the circuit.\n";
    } else {
        print_node_voltage(node_id);
    }
}

//解析打印变量，填充ckt.print_node_ids，以及sources中器件的printI标志
void solver::parse_print_variables() {
    ckt.print_node_ids.clear();
    ckt.plot_node_ids.clear();
    //Debug: 输出要解析的打印变量
    std::cout << analysis_type.print_variables.size() << " print variables to parse:\n";
    for (const auto& var : analysis_type.print_variables) {
        std::cout << "  " << var << "\n";
    }
    std::cout << analysis_type.plot_variables.size() << " plot variables to parse:\n";
    for (const auto& var : analysis_type.plot_variables) {
        std::cout << "  " << var << "\n";
    }

    for (const auto& var : analysis_type.print_variables) {
        if (var[0] == 'V' && var[1] == '(' && var.back() == ')') {
            //节点电压变量
            std::string node_name = var.substr(2, var.size() - 3);
            int node_id = ckt.getNodeID(node_name);
            if (node_id != -1) {
                ckt.print_node_ids.push_back(node_id);
            } else {
                std::cout << "Warning: Node " << node_name << " not found for printing voltage.\n";
            }
        } else if (var[0] == 'I' && var[1] == '(' && var.back() == ')') {
            //支路电流变量
            std::string device_name = var.substr(2, var.size() - 3);
            bool found = false;
            for (auto& src : ckt.sources) {
                if (src.name == device_name) {
                    src.printI = true; //标记该器件需要打印电流
                    found = true;
                    break;
                }
            }
            if (!found) {
                std::cout << "Warning: Source " << device_name << " not found for printing current.\n";
            }
        } else {
            std::cout << "Warning: Unrecognized print variable format: " << var << "\n";
        }
    }
    for (const auto& var : analysis_type.plot_variables) {
        if (var[0] == 'I' && var[1] == '(' && var.back() == ')') {
            //支路电流变量
            std::string device_name = var.substr(2, var.size() - 3);
            bool found = false;
            for (auto& src : ckt.sources) {
                if (src.name == device_name) {
                    src.plotI = true; //标记该器件需要打印电流
                    found = true;
                    break;
                }
            }
            if (!found) {
                std::cout << "Warning: Source " << device_name << " not found for printing current.\n";
            }
        }
        else{
            //节点电压变量
            std::string node_name = var;
            int node_id = ckt.getNodeID(node_name);
            if (node_id != -1) {
                ckt.plot_node_ids.push_back(node_id);
            } else {
                std::cout << "Warning: Node " << node_name << " not found for printing voltage.\n";
            }
        }
    }
}

//高斯消去法线性MNA方程求解
void solver::solve_linear_MNA_Gauss(){
    // //展示MNA矩阵和J向量
    // std::cout << "MNA Matrix:\n" << MNA_Y << "\n";
    // std::cout << "J Vector:\n" << J << "\n";
    int n = MNA_Y.rows();
    Eigen::VectorXd b = J;
    for (int k = 0; k < n; ++k){
        //部分主元选择
        int max_row = k;
        for (int i = k + 1; i < n; ++i){
            if (std::abs(MNA_Y(i, k)) > std::abs(MNA_Y(max_row, k))){
                max_row = i;
            }
        }
        //如果需要，交换行
        if (max_row != k){
            MNA_Y.row(k).swap(MNA_Y.row(max_row));
            std::swap(b(k), b(max_row));
        }
        // //消元过程
        // //Debug:输出交换后的MNA_Y
        // std::cout << MNA_Y << std::endl;
        // std::cout << b << std::endl;
        for (int i = k + 1; i < n; ++i){
            double factor = MNA_Y(i, k) / MNA_Y(k, k);
            MNA_Y.row(i) -= factor * MNA_Y.row(k);
            b(i) -= factor * b(k);
        }
    }
    //回代过程
    Eigen::VectorXd x(n);
    for (int i = n - 1; i >= 0; --i){
        double sum = b(i);
        //把已知的变量代入，移项到右边
        for (int j = i + 1; j < n; ++j){
            sum -= MNA_Y(i, j) * x(j);
        }
        x(i) = sum / MNA_Y(i, i);
    }
    //存储节点电压结果
    node_voltages.resize(ckt.node_map.size() - 1);
    for (int i = 0; i < ckt.node_map.size() - 1; ++i){
        node_voltages[i] = x(i);
    }
    //如果有额外的变量（如支路电流），也存储起来
    if (x.size() > ckt.node_map.size() - 1) {
        branch_currents.resize(x.size() - (ckt.node_map.size() - 1));
        for (int i = ckt.node_map.size() - 1; i < x.size(); ++i) {
            branch_currents[i - (ckt.node_map.size() - 1)] = x(i);
        }
        // std::cout << "Branch Currents:\n";
        // for (int i = 0; i < branch_currents.size(); ++i) {
        //     std::cout << "I" << i << ": " << branch_currents[i] << " A\n";
        // }
    }
}

//LU分解法求解MNA方程
void solver::solve_linear_MNA_LU(){
    // //展示MNA矩阵和J向量
    // std::cout << "MNA Matrix:\n" << MNA_Y << "\n";
    // std::cout << "J Vector:\n" << J << "\n";
    
    int n = MNA_Y.rows();
    if (n == 0) {
        std::cout << "Error: Empty MNA matrix\n";
        return;
    }
    
    //使用Eigen的LU分解求解器
    Eigen::PartialPivLU<Eigen::MatrixXd> lu(MNA_Y);
    
    //检查矩阵是否可逆
    if (lu.determinant() == 0.0) {
        std::cout << "Warning: Matrix is singular or near-singular\n";
    }
    
    //求解线性方程组 MNA_Y * x = J
    Eigen::VectorXd x = lu.solve(J);
    
    // //Debug检查解的精度
    // double relative_error = (MNA_Y * x - J).norm() / J.norm();
    // std::cout << "LU Decomposition Relative Error: " << relative_error << "\n";
    
    //存储节点电压结果
    node_voltages.resize(ckt.node_map.size() - 1);
    for (int i = 0; i < std::min((int)(ckt.node_map.size() - 1), (int)x.size()); ++i){
        node_voltages[i] = x(i);
    }
    
    //如果有额外的变量（如支路电流），也存储起来
    if (x.size() > ckt.node_map.size() - 1) {
        branch_currents.resize(x.size() - (ckt.node_map.size() - 1));
        for (int i = ckt.node_map.size() - 1; i < x.size(); ++i) {
            branch_currents[i - (ckt.node_map.size() - 1)] = x(i);
        }
        // std::cout << "Branch Currents:\n";
        // for (int i = 0; i < branch_currents.size(); ++i) {
        //     std::cout << "I" << i << ": " << branch_currents[i] << " A\n";
        // }
    }
}

//完整LU分解法（手动实现）
void solver::get_linear_MNA_LU_manual(){
    // std::cout << "Manual LU Decomposition:\n";
    // std::cout << "Original MNA Matrix:\n" << MNA_Y << "\n";
    // std::cout << "Original J Vector:\n" << J << "\n";
    
    int n = MNA_Y.rows();
    if (n == 0) {
        std::cout << "Error: Empty MNA matrix\n";
        return;
    }
    
    //初始化L和U矩阵
    L = Eigen::MatrixXd::Identity(n, n);  //L矩阵初始为单位矩阵
    U = Eigen::MatrixXd::Zero(n, n);      //U矩阵初始为零矩阵
    
    //复制原矩阵到临时矩阵A进行操作
    Eigen::MatrixXd A = MNA_Y;
    //同样需要复制J向量，因为主元交换时需要同步交换
    Eigen::VectorXd b = J;
    
    //部分主元选择的置换向量
    std::vector<int> pivot_order(n);
    for (int i = 0; i < n; ++i) pivot_order[i] = i;
    
    //进行LU分解
    for (int k = 0; k < n; ++k) {
        //部分主元选择
        int max_row = k;
        for (int i = k + 1; i < n; ++i) {
            if (std::abs(A(i, k)) > std::abs(A(max_row, k))) {
                max_row = i;
            }
        }
        
        //如果需要，交换行
        if (max_row != k) {
            // 只交换A矩阵和b向量的行，不交换L矩阵
            A.row(k).swap(A.row(max_row));
            std::swap(b(k), b(max_row));
            std::swap(pivot_order[k], pivot_order[max_row]);
            
            // 交换L矩阵中已计算的部分（第k列之前的元素）
            for (int j = 0; j < k; ++j) {
                std::swap(L(k, j), L(max_row, j));
            }
            
            // std::cout << "Swapping rows " << k << " and " << max_row << " (including J vector)\n";
        }
        
        //检查主元是否接近零
        if (std::abs(A(k, k)) < 1e-12) {
            std::cout << "Warning: Near-zero pivot at position (" << k << "," << k << "): " << A(k, k) << "\n";
        }
        
        //计算U矩阵的第k行
        for (int j = k; j < n; ++j) {
            U(k, j) = A(k, j);
        }
        
        //计算L矩阵的第k列并进行消元
        for (int i = k + 1; i < n; ++i) {
            if (std::abs(U(k, k)) > 1e-12) {
                L(i, k) = A(i, k) / U(k, k);
                //消元：第i行减去L(i,k)*第k行
                for (int j = k; j < n; ++j) {
                    A(i, j) -= L(i, k) * U(k, j);
                }
            } else {
                std::cout << "Warning: Division by near-zero pivot\n";
                L(i, k) = 0.0;
            }
        }
    }
    
    //存储经过置换的J向量，供后续求解使用
    J_permuted = b;
    
    // //验证LU分解的正确性
    // Eigen::MatrixXd product = L * U;
    // Eigen::MatrixXd error_matrix = product - MNA_Y;
    // double max_error = error_matrix.cwiseAbs().maxCoeff();
    
    // std::cout << "LU Decomposition Results:\n";
    // std::cout << "L Matrix:\n" << L << "\n";
    // std::cout << "U Matrix:\n" << U << "\n";
    // std::cout << "L*U Product:\n" << product << "\n";
    // std::cout << "Permuted J Vector:\n" << J_permuted << "\n";
    // std::cout << "Maximum decomposition error: " << max_error << "\n";
    
    // if (max_error < 1e-10) {
    //     std::cout << "LU decomposition successful!\n";
    // } else {
    //     std::cout << "Warning: LU decomposition may have numerical issues\n";
    // }
}

//使用已计算的L和U矩阵求解MNA方程
void solver::solve_with_LU_matrices(){
    if (L.rows() == 0 || U.rows() == 0) {
        std::cout << "Error: L or U matrix not computed. Call get_linear_MNA_LU_manual() first.\n";
        return;
    }
    
    if (J_permuted.size() == 0) {
        std::cout << "Error: Permuted J vector not available. Call get_linear_MNA_LU_manual() first.\n";
        return;
    }
    
    int n = L.rows();
    // std::cout << "Solving using pre-computed L and U matrices:\n";
    // std::cout << "Permuted J Vector:\n" << J_permuted << "\n";
    
    //第一步：前向替代求解 L*y = J_permuted
    Eigen::VectorXd y(n);
    for (int i = 0; i < n; ++i) {
        double sum = J_permuted(i);  //使用经过置换的J向量
        for (int j = 0; j < i; ++j) {
            sum -= L(i, j) * y(j);
        }
        y(i) = sum / L(i, i);
    }
    
    // std::cout << "Intermediate vector y (L*y = J_permuted):\n" << y << "\n";
    
    //第二步：后向替代求解 U*x = y
    Eigen::VectorXd x(n);
    for (int i = n - 1; i >= 0; --i) {
        double sum = y(i);
        for (int j = i + 1; j < n; ++j) {
            sum -= U(i, j) * x(j);
        }
        x(i) = sum / U(i, i);
    }
    
    // std::cout << "Solution vector x (U*x = y):\n" << x << "\n";
    
    //验证解的正确性（使用原始的MNA矩阵和J向量）
    Eigen::VectorXd residual = MNA_Y * x - J;
    double residual_norm = residual.norm();
    // std::cout << "Residual norm ||MNA_Y*x - J||: " << residual_norm << "\n";
    
    //存储节点电压结果
    node_voltages.resize(std::min((int)(ckt.node_map.size() - 1), (int)x.size()));
    for (int i = 0; i < node_voltages.size(); ++i){
        node_voltages[i] = x(i);
    }
    // //Debug: 输出节点电压
    // std::cout << "Node Voltages:\n";
    // for (int i = 0; i < node_voltages.size(); ++i) {
    //     std::cout << "V" << i+1 << ": " << node_voltages[i] << " V\n";
    // }
    
    //如果有额外的变量（如支路电流），也存储起来
    if (x.size() > ckt.node_map.size() - 1) {
        branch_currents.resize(x.size() - (ckt.node_map.size() - 1));
        for (int i = ckt.node_map.size() - 1; i < x.size(); ++i) {
            branch_currents[i - (ckt.node_map.size() - 1)] = x(i);
        }
        // std::cout << "Branch Currents:\n";
        // for (int i = 0; i < branch_currents.size(); ++i) {
        //     std::cout << "I" << i << ": " << branch_currents[i] << " A\n";
        // }
    }
}

//Gauss-Jacobi迭代法求解线性MNA方程
void solver::solve_linear_MNA_Gauss_Jacobi(){
    // //Debug：打印MNA矩阵和J向量
    // std::cout << "MNA_Y:\n" << MNA_Y << "\n";
    // std::cout << "J:\n" << J << "\n";
    //先进行行交换，确保对角线元素不为零
    int n = MNA_Y.rows();
    // for (int k = 0; k < n; ++k){
    //     if (std::abs(MNA_Y(k, k)) < 1e-6){
    //         //寻找所有行中绝对值最大的元素进行交换
    //         int max_row = k;
    //         for (int i = 0; i < n; ++i){
    //             if (std::abs(MNA_Y(i, k)) > std::abs(MNA_Y(max_row, k))){
    //                 max_row = i;
    //             }
    //         }
    //         if (max_row != k){
    //             MNA_Y.row(k).swap(MNA_Y.row(max_row));
    //             std::swap(J(k), J(max_row));
    //             // std::cout << "Swapping rows " << k << " and " << max_row << " to avoid zero diagonal\n";
    //         }
    //     }
    // }

    // //Debug:用LU分解法先得到解
    // //使用Eigen的LU分解求解器
    // Eigen::PartialPivLU<Eigen::MatrixXd> lu(MNA_Y);
    
    // //检查矩阵是否可逆
    // if (lu.determinant() == 0.0) {
    //     std::cout << "Warning: Matrix is singular or near-singular\n";
    // }
    
    // //求解线性方程组 MNA_Y * x = J
    // Eigen::VectorXd x = lu.solve(J);

    // std::cout << "LU Decomposition Solution x:\n" << x << "\n";
    // /////////////////////////////////////////////////////////////

    // //使用对角占优化处理，增强收敛性
    // Eigen::MatrixXd inversePerm = diagMax();
    // //Debug: 输出MNA矩阵和J向量在主元交换后的形式
    // std::cout << "MNA_Y Matrix before pivoting:\n" << MNA_Y << "\n";

    //对于MNA(i,i)为零的情况，必然是电压源，寻找上方为1或者-1的行，哪个行对应的对角元小，与那行交换
    for (int i = 0; i < n; ++i){
        if (std::abs(MNA_Y(i, i)) < 1e-12){
            //寻找上方为1或者-1的行
            int row1 = -1;
            int row2 = -1;
            for (int j = 0; j < n; ++j){
                if (j != i && (std::abs(MNA_Y(j, i) - 1.0) < 1e-12 || std::abs(MNA_Y(j, i) + 1.0) < 1e-12)){
                    if (row1 == -1){
                        row1 = j;
                    } else {
                        row2 = j;
                        break;
                    }
                }
            }
            //选择对角元较小的行进行交换
            int swap_row = -1;
            if (row1 != -1 && row2 != -1){
                if (std::abs(MNA_Y(row1, row1)) < std::abs(MNA_Y(row2, row2))){
                    swap_row = row1;
                } else {
                    swap_row = row2;
                }
            } else if (row1 != -1){
                swap_row = row1;
            } else if (row2 != -1){
                swap_row = row2;
            }
            if (swap_row != -1){
                MNA_Y.row(i).swap(MNA_Y.row(swap_row));
                std::swap(J(i), J(swap_row));
                // std::cout << "Swapping rows " << i << " and " << swap_row << " to avoid zero diagonal\n";
            }
        }
    }


    // //Debug: 输出调整后的MNA矩阵和J向量
    // std::cout << "Adjusted MNA_Y Matrix:\n" << MNA_Y << "\n";
    // std::cout << "Adjusted J Vector:\n" << J << "\n";

    //构造迭代矩阵，判断收敛性
    Eigen::MatrixXd D = MNA_Y.diagonal().asDiagonal();
    Eigen::MatrixXd D_inv = D.inverse();
    Eigen::MatrixXd L_plus_U = MNA_Y - D;
    Eigen::MatrixXd iteration_matrix = -D_inv * L_plus_U;
    // //Debug: 输出D矩阵和D的逆
    // std::cout << "D Matrix:\n" << D << "\n";
    // std::cout << "D Inverse Matrix:\n" << D_inv << "\n";
    // //Debug: 输出L+U矩阵
    // std::cout << "L + U Matrix:\n" << L_plus_U << "\n";
    // //Debug: 输出迭代矩阵
    // std::cout << "Iteration Matrix:\n" << iteration_matrix << "\n";
    //计算特征值判断收敛性
    Eigen::EigenSolver<Eigen::MatrixXd> iteration_es(iteration_matrix);
    // std::cout << "Eigenvalues of iteration matrix:\n" << iteration_es.eigenvalues() << "\n";
    double max_eigenvalue = 0.0;
    for (int i = 0; i < iteration_es.eigenvalues().size();i++){
        double abs_eigenvalue = std::abs(iteration_es.eigenvalues()(i));
        if (abs_eigenvalue > max_eigenvalue){
            max_eigenvalue = abs_eigenvalue;
        }
    }
    if(max_eigenvalue >= 1.0){
        std::cout << "Warning: Gauss-Jacobi iteration may not converge (max eigenvalue = " << max_eigenvalue << ")\n";
    } else {
        std::cout << "Gauss-Jacobi iteration expected to converge (max eigenvalue = " << max_eigenvalue << ")\n";
    }

    // // std::cout << "Solving MNA equations using Gauss-Jacobi Iteration:\n";
    // // //展示特征值
    // Eigen::EigenSolver<Eigen::MatrixXd> es(MNA_Y);
    // // std::cout << "Eigenvalues of MNA_Y:\n" << es.eigenvalues() << "\n";
    // // //如果存在特征值绝对值大于1，则可能不收敛
    // for (int i = 0; i < es.eigenvalues().size(); ++i) {
    //     if (std::abs(es.eigenvalues()(i)) > 1.0) {
    //         std::cout << "Warning: Eigenvalue " << es.eigenvalues()(i) << " may indicate non-convergence\n";
    //         break;
    //     }
    // }

    Eigen::VectorXd x_old = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd x_new = Eigen::VectorXd::Zero(n);
    const int max_iterations = 1000;
    const double tolerance = 1e-10;
    
    for (int iter = 0; iter < max_iterations; ++iter) {

        for (int i = 0; i < n; ++i) {
            double sum = J(i);
            for (int j = 0; j < n; ++j) {
                if (j != i) {
                    sum -= MNA_Y(i, j) * x_old(j);
                }
            }
            if (std::abs(MNA_Y(i, i)) < 1e-12) {
                // std::cout << "Warning: Near-zero diagonal element at row " << i << "\n";
                x_new(i) = 0.0;
            } else {
                x_new(i) = sum / MNA_Y(i, i);
            }
        }
        
        // //Debug: 输出每次迭代的结果
        // std::cout << "Iteration " << iter + 1 << ": x = " << x_new << "\n";

        //检查收敛性
        double error = (x_new - x_old).norm();
        if (error < tolerance) {
            std::cout << "Converged in " << iter + 1 << " iterations with error: " << error << "\n";
            break;
        }
        
        x_old = x_new;
        
        if (iter == max_iterations - 1) {
            std::cout << "Warning: Did not converge within the maximum number of iterations\n";
        }
    }
    // //交换回来
    // x_new = inversePerm * x_new;
    // 3. 恢复解的顺序
    
    // //Debug: 输出最终解
    // std::cout << "Gauss-Jacobi Solution x:\n" << x_new << "\n";
    
    // //存储节点电压结果
    // node_voltages.resize(ckt.node_map.size() - 1);
    // for (int i = 0; i < ckt.node_map.size() - 1; ++i){
    //     node_voltages[i] = x_new(i);
    // }
    // //如果有额外的变量（如支路电流），也存储起来
    // if (x_new.size() > ckt.node_map.size() - 1) {
    //     branch_currents.resize(x_new.size() - (ckt.node_map.size() - 1));
    //     for (int i = ckt.node_map.size() - 1; i < x_new.size(); ++i) {
    //         branch_currents[i - (ckt.node_map.size() - 1)] = x_new(i);
    //     }
    // }
   //////////////////////////////////////
    //存储节点电压结果
    node_voltages.resize(std::min((int)(ckt.node_map.size() - 1), (int)x_new.size()));
    for (int i = 0; i < node_voltages.size(); ++i){
        node_voltages[i] = x_new(i);
    }
    // //Debug: 输出节点电压
    // std::cout << "Node Voltages:\n";
    // for (int i = 0; i < node_voltages.size(); ++i) {
    //     std::cout << "V" << i+1 << ": " << node_voltages[i] << " V\n";
    // }
    
    //如果有额外的变量（如支路电流），也存储起来
    if (x_new.size() > ckt.node_map.size() - 1) {
        branch_currents.resize(x_new.size() - (ckt.node_map.size() - 1));
        for (int i = ckt.node_map.size() - 1; i < x_new.size(); ++i) {
            branch_currents[i - (ckt.node_map.size() - 1)] = x_new(i);
        }
        // std::cout << "Branch Currents:\n";
        // for (int i = 0; i < branch_currents.size(); ++i) {
        //     std::cout << "I" << i << ": " << branch_currents[i] << " A\n";
        // }
    }
}

//Gauss-Seidel迭代法求解线性MNA方程
void solver::solve_linear_MNA_Gauss_Seidel(){
    //先进行行交换，确保对角线元素不为零
    int n = MNA_Y.rows();

    //对于MNA(i,i)为零的情况，必然是电压源，寻找上方为1或者-1的行，哪个行对应的对角元小，与那行交换
    for (int i = 0; i < n; ++i){
        if (std::abs(MNA_Y(i, i)) < 1e-12){
            //寻找上方为1或者-1的行
            int row1 = -1;
            int row2 = -1;
            for (int j = 0; j < n; ++j){
                if (j != i && (std::abs(MNA_Y(j, i) - 1.0) < 1e-12 || std::abs(MNA_Y(j, i) + 1.0) < 1e-12)){
                    if (row1 == -1){
                        row1 = j;
                    } else {
                        row2 = j;
                        break;
                    }
                }
            }
            //选择对角元较小的行进行交换
            int swap_row = -1;
            if (row1 != -1 && row2 != -1){
                if (std::abs(MNA_Y(row1, row1)) < std::abs(MNA_Y(row2, row2))){
                    swap_row = row1;
                } else {
                    swap_row = row2;
                }
            } else if (row1 != -1){
                swap_row = row1;
            } else if (row2 != -1){
                swap_row = row2;
            }
            if (swap_row != -1){
                MNA_Y.row(i).swap(MNA_Y.row(swap_row));
                std::swap(J(i), J(swap_row));
                // std::cout << "Swapping rows " << i << " and " << swap_row << " to avoid zero diagonal\n";
            }
        }
    }
    Eigen::VectorXd x_old = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd x_new = Eigen::VectorXd::Zero(n);
    const int max_iterations = 5000;
    const double tolerance = 1e-9;
    for (int iter = 0; iter < max_iterations; ++iter) {
        x_new = x_old; //每次迭代开始时，先复制旧值
        for (int i = 0; i < n; ++i) {
            double sum = J(i);
            for (int j = 0; j < n; ++j) {
                if (j != i) {
                    sum -= MNA_Y(i, j) * x_new(j); //使用最新的x_new值
                }
            }
            if (std::abs(MNA_Y(i, i)) < 1e-12) {
                // std::cout << "Warning: Near-zero diagonal element at row " << i << "\n";
                x_new(i) = 0.0;
            } else {
                x_new(i) = sum / MNA_Y(i, i);
            }
        }
        //检查收敛性
        double error = (x_new - x_old).norm();
        if (error < tolerance) {
            std::cout << "Converged in " << iter + 1 << " iterations with error: " << error << "\n";
            break;
        }
        x_old = x_new;
        if (iter == max_iterations - 1) {
            std::cout << "Warning: Did not converge within the maximum number of iterations\n";
        }
    }
    //存储节点电压结果
    node_voltages.resize(std::min((int)(ckt.node_map.size() - 1), (int)x_new.size()));
    for (int i = 0; i < node_voltages.size(); ++i){
        node_voltages[i] = x_new(i);
    }
    // //Debug: 输出节点电压
    // std::cout << "Node Voltages:\n";
    // for (int i = 0; i < node_voltages.size(); ++i) {
    //     std::cout << "V" << i+1 << ": " << node_voltages[i] << " V\n";
    // }
    
    //如果有额外的变量（如支路电流），也存储起来
    if (x_new.size() > ckt.node_map.size() - 1) {
        branch_currents.resize(x_new.size() - (ckt.node_map.size() - 1));
        for (int i = ckt.node_map.size() - 1; i < x_new.size(); ++i) {
            branch_currents[i - (ckt.node_map.size() - 1)] = x_new(i);
        }
        // std::cout << "Branch Currents:\n";
        // for (int i = 0; i < branch_currents.size(); ++i) {
        //     std::cout << "I" << i << ": " << branch_currents[i] << " A\n";
        // }
    }
}

//根据设置的方法选择线性方程求解算法
void solver::solve_linear_MNA(){
    switch (linear_solver_method) {
        case LinearSolverMethod::GAUSS_ELIMINATION:
            solve_linear_MNA_Gauss();
            break;
        case LinearSolverMethod::LU_DECOMPOSITION:
            solve_linear_MNA_LU();
            break;
        case LinearSolverMethod::MANUAL_LU:
            get_linear_MNA_LU_manual();
            solve_with_LU_matrices();
            break;
        case LinearSolverMethod::GAUSS_JACOBI:
            solve_linear_MNA_Gauss_Jacobi();
            break;
        case LinearSolverMethod::GAUSS_SEIDEL:
            solve_linear_MNA_Gauss_Seidel();
            break;
        default:
            std::cout << "Error: Unknown linear solver method\n";
            solve_linear_MNA_LU(); // 默认使用LU分解
            break;
    }
}

//method: 0-高斯消去法，1-LU分解法，2-手动LU分解法，3-Gauss-Jacobi迭代法
void solver::solve_linear_MNA(int method){
    switch (method) {
        case 0:
            solve_linear_MNA_Gauss();
            break;
        case 1:
            solve_linear_MNA_LU();
            break;
        case 2:
            get_linear_MNA_LU_manual();
            solve_with_LU_matrices();
            break;
        case 3:
            solve_linear_MNA_Gauss_Jacobi();
            break;
        default:
            std::cout << "Error: Unknown method " << method << "\n";
            break;
    }
}

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


            // //引入辅助变量和KVL方程
            // if (n1 == 0){
            //     //节点n2电压为-value
            //     J(n2 - 1) = -value;
            //     MNA_Y.row(n2 - 1).setZero();
            //     MNA_Y(n2 - 1, n2 - 1) = 1;
            //     //引入支路电流变量与新方程

            // }
            // else if (n2 == 0){
            //     //节点n1电压为value
            //     J(n1 - 1) = value;
            //     MNA_Y.row(n1 - 1).setZero();
            //     MNA_Y(n1 - 1, n1 - 1) = 1;
            // }
            // else{
            //     //引入支路电流变量与新方程
            //     int new_var_index = MNA_Y.rows();
            //     MNA_Y.conservativeResize(new_var_index + 1, new_var_index + 1);
            //     MNA_Y.row(new_var_index).setZero();
            //     MNA_Y.col(new_var_index).setZero();
            //     J.conservativeResize(new_var_index + 1);
            //     J(new_var_index) = value;
            //     //KVL方程
            //     MNA_Y(new_var_index, n1 - 1) = 1;
            //     MNA_Y(new_var_index, n2 - 1) = -1;
            //     //支路电流对节点的贡献
            //     MNA_Y(n1 - 1, new_var_index) = 1;
            //     MNA_Y(n2 - 1, new_var_index) = -1;
            // }
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


            // //引入辅助变量和KVL方程
            // if (n1 == 0){
            //     //节点n2电压为-value
            //     J(n2 - 1) = -value;
            //     MNA_Y.row(n2 - 1).setZero();
            //     MNA_Y(n2 - 1, n2 - 1) = 1;
            //     //引入支路电流变量与新方程

            // }
            // else if (n2 == 0){
            //     //节点n1电压为value
            //     J(n1 - 1) = value;
            //     MNA_Y.row(n1 - 1).setZero();
            //     MNA_Y(n1 - 1, n1 - 1) = 1;
            // }
            // else{
            //     //引入支路电流变量与新方程
            //     int new_var_index = MNA_Y.rows();
            //     MNA_Y.conservativeResize(new_var_index + 1, new_var_index + 1);
            //     MNA_Y.row(new_var_index).setZero();
            //     MNA_Y.col(new_var_index).setZero();
            //     J.conservativeResize(new_var_index + 1);
            //     J(new_var_index) = value;
            //     //KVL方程
            //     MNA_Y(new_var_index, n1 - 1) = 1;
            //     MNA_Y(new_var_index, n2 - 1) = -1;
            //     //支路电流对节点的贡献
            //     MNA_Y(n1 - 1, new_var_index) = 1;
            //     MNA_Y(n2 - 1, new_var_index) = -1;
            // }
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

//直流分析
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
            if(plot_current_dev_index >=0 && plot_current_dev_index < ckt.sources.size()){
                double i = branch_currents[plot_current_dev_index];
                tran_plot_data[node_voltages.size() + plot_current_dev_index].push_back(std::make_pair(time, i));
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
                if(current_dev_index >=0 && current_dev_index < ckt.sources.size()){
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
            if(plot_current_dev_index >=0 && plot_current_dev_index < ckt.sources.size()){
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
                if(current_dev_index >=0 && current_dev_index < ckt.sources.size()){
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
void solver::hb_build_linear_MNA(){
    //先对直流点构建线性MNA矩阵,只能进行一次，电感会贴出来很多个电压源
    // build_linear_MNA(false);
    //确保MNA_Y大小正确
    MNA_Y = Eigen::MatrixXd::Zero(liner_Y.rows(), liner_Y.cols());
    MNA_Y = liner_Y;
    build_sources_MNA();
    //这样就得到的MNA_Y大小就是每个频率点的矩阵大小
    base_size = MNA_Y.rows();
    //初始化多频率点的线性MNA矩阵,包含正负频率点
    //每个频率点的矩阵大小为base_size * (num_harmonics + 1)
    //这里的num_harmonics是指正频率点的数量，不包括直流点
    hb_liner_Y = Eigen::MatrixXcd::Zero(base_size * (2*hb_params.num_harmonics+1), base_size * (2*hb_params.num_harmonics+1));
    //初始化多频率点的线性J向量
    hb_J = Eigen::VectorXcd::Zero(base_size * (2*hb_params.num_harmonics+1));
    int node_list_size = ckt.node_list.size() -1; //不包括地节点
    //把直流点的线性MNA矩阵放入多频率点矩阵的是中央位置
    hb_liner_Y.block(hb_params.num_harmonics * base_size, hb_params.num_harmonics * base_size, base_size, base_size) = MNA_Y.cast<std::complex<double> >();
    hb_J.segment(hb_params.num_harmonics * base_size, base_size) = J.cast<std::complex<double> >();

    //对各个频率点进行线性MNA矩阵构建
    for(int h = -hb_params.num_harmonics; h <= hb_params.num_harmonics; ++h){
        if(h == 0) continue; //直流点已经处理过了
        double omega;
        omega = h * hb_params.fundamental_omega;
        std::complex<double> jw(0.0, omega);
        //初始化子MNA为0
        Eigen::MatrixXcd Y_h = Eigen::MatrixXcd::Zero(base_size,base_size);
        //对电容和电感进行频率域修正
        for (const auto &dev : ckt.linear_devices) {
            char c = toupper(dev.name[0]);
            if(c == 'R'){
                int n1 = dev.nodes[0];
                int n2 = dev.nodes[1];
                if (dev.parameters.find("value") == dev.parameters.end()) {
                    continue; // 跳过没有value参数的器件
                }
                double R = dev.parameters.at("value");
                if (R <= 0) {
                    continue; // 跳过无效的电阻值
                }
                std::complex<double> Yr = 1.0 / R; // 电阻的导纳
                if (n1 != 0) Y_h(n1 - 1, n1 - 1) += Yr;
                if (n2 != 0) Y_h(n2 - 1, n2 - 1) += Yr;
                if (n1 != 0 && n2 != 0) {
                    Y_h(n1 - 1, n2 - 1) -= Yr;
                    Y_h(n2 - 1, n1 - 1) -= Yr;
                }
            }
            else if (c == 'C') {
                int n1 = dev.nodes[0];
                int n2 = dev.nodes[1];
                if (dev.parameters.find("value") == dev.parameters.end()) {
                    continue; // 跳过没有value参数的器件
                }
                double C = dev.parameters.at("value");
                if (C <= 0) {
                    continue; // 跳过无效的电容值
                }
                std::complex<double> Yc = jw * C; // 电容的导纳
                if (n1 != 0) Y_h(n1 - 1, n1 - 1) += Yc;
                if (n2 != 0) Y_h(n2 - 1, n2 - 1) += Yc;
                if (n1 != 0 && n2 != 0) {
                    Y_h(n1 - 1, n2 - 1) -= Yc;
                    Y_h(n2 - 1, n1 - 1) -= Yc;
                }
            }
            else if (c == 'L') {
                int n1 = dev.nodes[0];
                int n2 = dev.nodes[1];
                if (dev.parameters.find("value") == dev.parameters.end()) {
                    continue; // 跳过没有value参数的器件
                }
                double L = dev.parameters.at("value");
                if (L <= 0) {
                    continue; // 跳过无效的电感值
                }
                //电感引入电流支路，确定支路电流的行列号
                int branch_index = dev.branch_current_index + node_list_size; //支路电流在MNA矩阵中的行列号
                if (n1 != 0){
                    Y_h(n1 - 1, branch_index) = 1;
                    Y_h(branch_index, n1 - 1) = 1;
                } 
                if (n2 != 0){
                    Y_h(n2 - 1, branch_index) = -1;
                    Y_h(branch_index, n2 - 1) = -1;
                }
                Y_h(branch_index, branch_index) = -jw * L;
            }
        }
        //将该频率点的线性MNA矩阵放入多频率点矩阵中
        hb_liner_Y.block((h + hb_params.num_harmonics) * base_size, (h + hb_params.num_harmonics) * base_size, base_size, base_size) = Y_h;
    }
    // //Debug: 输出多频率点的线性MNA矩阵和J向量
    // //输出到文件中查看
    // std::ofstream out("hb_linear_MNA.txt");
    // if (out.is_open()) {
    //     out << "Harmonic Balance Linear MNA Matrix (Y):\n" << hb_liner_Y << "\n";
    //     out << "Harmonic Balance Linear J Vector:\n" << hb_J << "\n";
    //     out.close();
    // } else {
    //     std::cerr << "Error opening file for writing: hb_linear_MNA.txt\n";
    // }
    // std::cout << "Harmonic Balance Linear J Vector:\n" << hb_J << "\n";
    // std::cout << "Harmonic Balance Linear MNA Matrix (Y):\n" << hb_liner_Y << "\n";
}

//构建普通的DFT和IDFT变换矩阵
void solver::hb_initialize_DFT_matrices() {
    int num_harmonics = hb_params.num_harmonics;
    int N = 2 * num_harmonics + 1; // 总频率点数 = 时域点数
    double f0 = 1.0; // 基频（归一化）
    double T = 1.0 / f0; // 周期
    
    hb_DFT_matrix = Eigen::MatrixXcd::Zero(N, N);
    hb_iDFT_matrix = Eigen::MatrixXcd::Zero(N, N);
    
    // ============= DFT矩阵 (时域→频域) =============
    // 公式：DFT(k, n) = exp(-j·2π·f_k·t_n) / N
    // 其中 f_k = (k - num_harmonics) * f0
    //      t_n = n * T / N
    
    for(int k = 0; k < N; ++k) {
        // 频率索引 m = -num_harmonics, ..., 0, ..., +num_harmonics
        int m = k - num_harmonics;
        double freq = m * f0;
        
        for(int n = 0; n < N; ++n) {
            // 时间点 t = n * T / N
            double t = n * T / N;
            double theta = -2.0 * M_PI * freq * t;
            hb_DFT_matrix(k, n) = std::exp(std::complex<double>(0.0, theta)) / static_cast<double>(N);
        }
    }
    
    // ============= IDFT矩阵 (频域→时域) =============
    // 公式：IDFT(n, k) = exp(+j·2π·f_k·t_n)
    // 注意：DFT和IDFT通常是互逆的，所以IDFT不需要除以N
    
    for(int n = 0; n < N; ++n) {
        double t = n * T / N;
        
        for(int k = 0; k < N; ++k) {
            int m = k - num_harmonics;
            double freq = m * f0;
            double theta = 2.0 * M_PI * freq * t;
            hb_iDFT_matrix(n, k) = std::exp(std::complex<double>(0.0, theta));
        }
    }
    
}

//构建FT变换矩阵
void solver::hb_build_TF_matrix(){

    int N = 2 * hb_params.num_harmonics + 1;
    //初始化TF矩阵和iTF矩阵
    hb_T2F_matrix = Eigen::MatrixXcd::Zero(base_size * N, base_size * N);
    hb_F2T_matrix = Eigen::MatrixXcd::Zero(base_size * N, base_size * N);
    //构建TF矩阵
    for(int i = 0; i < base_size; ++i){
        for(int k = 0; k < N; ++k){
            for(int n = 0; n < N; ++n){
                hb_T2F_matrix(i + k * base_size, i + n * base_size) = hb_DFT_matrix(k, n);
                hb_F2T_matrix(i + n * base_size, i + k * base_size) = hb_iDFT_matrix(n, k);
            }
        }
    }
}

//实现将时域的解进行DFT变换
Eigen::VectorXcd solver::hb_DFT(Eigen::VectorXcd xt){
    //xt长度为base_size * (2*hb_params.num_harmonics+1)，且同一个时刻的变量放在一起，同一变量不同时刻的值相距base_size
    int N = 2 * hb_params.num_harmonics + 1;
    Eigen::VectorXcd result(base_size * N);
    for(int i = 0; i < base_size; ++i){
        //提取第i个变量的时域序列
        Eigen::VectorXcd x_t(N);
        for(int n = 0; n < N; ++n){
            x_t(n) = xt(i + n * base_size);
        }
        //进行DFT变换
        Eigen::VectorXcd x_w = hb_DFT_matrix * x_t;
        //放回结果向量
        for(int k = 0; k < N; ++k){
            result(i + k * base_size) = x_w(k);
        }
    }
    return result;
}

//实现将HB的解进行逆DFT变换
Eigen::VectorXcd solver::hb_iDFT(Eigen::VectorXcd xw)
{
    int N = 2 * hb_params.num_harmonics + 1;
    Eigen::VectorXcd result(base_size * N);
    for(int i = 0; i < base_size; ++i){
        //提取第i个变量的频域序列
        Eigen::VectorXcd x_w(N);
        for(int k = 0; k < N; ++k){
            x_w(k) = xw(i + k * base_size);
        }
        //进行IDFT变换
        Eigen::VectorXcd x_t = hb_iDFT_matrix * x_w;
        //放回结果向量
        for(int n = 0; n < N; ++n){
            result(i + n * base_size) = x_t(n);
        }
    }
    return result;
}

//加入sources到MNA和J中
void solver::hb_build_sources_MNA(){

    //确保hb_MNA_Y大小正确,此时hb_liner_Y已经构建完成
    hb_MNA_Y = hb_liner_Y;
    for (auto &dev : ckt.sources){
        char c = toupper(dev.name[0]);
        //独立电压源
        if (c == 'V'){
            int n1 = dev.nodes[0];
            int n2 = dev.nodes[1];
            //判断是否是sin源
            if(dev.type == "V_SIN"){
                //根据其频率，计算对应的频率点索引
                double freq = dev.parameters["FREQ"];
                int harmonic_index = static_cast<int>(std::round(freq / (hb_params.fundamental_omega / (2.0 * M_PI))));
                if(harmonic_index < -hb_params.num_harmonics || harmonic_index > hb_params.num_harmonics){
                    //频率点超出范围，跳过
                    continue;
                }
                double vdc = dev.parameters["DC"];
                double amplitude = dev.parameters["AMP"];
                double phase_deg = dev.parameters["PHASE"];
                double phase_rad = phase_deg * M_PI / 180.0 - M_PI/2.0; //转换为弧度，并减去90度，变为正弦波初相位
                //在对应的正负频率点加入电压源贡献
                int pos_index = (harmonic_index + hb_params.num_harmonics) * base_size;
                int neg_index = (-harmonic_index + hb_params.num_harmonics) * base_size;
                //在hb_MNA_Y中加入电压源支路,在DC分析中，已经得到了支路电流变量索引
                int branch_index = dev.branch_current_index + (ckt.node_list.size() - 1); //支路电流在MNA矩阵中的行列号
                //对所有频率遍历，加入KVL方程
                for(int h = -hb_params.num_harmonics; h <= hb_params.num_harmonics; ++h){
                    int row_index = (h + hb_params.num_harmonics) * base_size + branch_index;
                    //直流点单独处理
                    //同时对于Jacobian矩阵也加入
                    if(h == 0){
                        if (n1 != 0){
                            hb_MNA_Y(row_index, (h + hb_params.num_harmonics) * base_size + n1 - 1) = 1;
                            hb_MNA_Y((h + hb_params.num_harmonics) * base_size + n1 - 1, row_index) = 1;
                            hb_jacobian(row_index, (h + hb_params.num_harmonics) * base_size + n1 - 1) = 1;
                            hb_jacobian((h + hb_params.num_harmonics) * base_size + n1 - 1, row_index) = 1;
                        }
                        if (n2 != 0){
                            hb_MNA_Y(row_index, (h + hb_params.num_harmonics) * base_size + n2 - 1) = -1;
                            hb_MNA_Y((h + hb_params.num_harmonics) * base_size + n2 - 1, row_index) = -1;
                            hb_jacobian(row_index, (h + hb_params.num_harmonics) * base_size + n2 - 1) = -1;
                            hb_jacobian((h + hb_params.num_harmonics) * base_size + n2 - 1, row_index) = -1;
                        }
                        //在J向量中加入电压源值
                        hb_J(row_index) = vdc;
                        continue;
                    }
                    if (n1 != 0){
                        hb_MNA_Y(row_index, (h + hb_params.num_harmonics) * base_size + n1 - 1) = 1;
                        hb_MNA_Y((h + hb_params.num_harmonics) * base_size + n1 - 1, row_index) = 1;
                        hb_jacobian(row_index, (h + hb_params.num_harmonics) * base_size + n1 - 1) = 1;
                        hb_jacobian((h + hb_params.num_harmonics) * base_size + n1 - 1, row_index) = 1;
                    }
                    if (n2 != 0){
                        hb_MNA_Y(row_index, (h + hb_params.num_harmonics) * base_size + n2 - 1) = -1;
                        hb_MNA_Y((h + hb_params.num_harmonics) * base_size + n2 - 1, row_index) = -1;
                        hb_jacobian(row_index, (h + hb_params.num_harmonics) * base_size + n2 - 1) = -1;
                        hb_jacobian((h + hb_params.num_harmonics) * base_size + n2 - 1, row_index) = -1;
                    }
                    //在J向量中加入电压源值，只在对应频率点加入，其余为0
                    if(h == harmonic_index){
                        hb_J(row_index) = std::polar(amplitude / 2.0, phase_rad); //正频率点
                    }
                    else if(h == -harmonic_index){
                        hb_J(row_index) = std::polar(amplitude / 2.0, -phase_rad); //负频率点
                    }
                    else{
                        hb_J(row_index) = 0;
                    }
                }
                
            }
            else if(dev.type == "V_DC"){
                //直流电压源，放在直流点
                double value = dev.parameters["DC"];
                int branch_index = dev.branch_current_index + (ckt.node_list.size() - 1); //支路电流在MNA矩阵中的行列号
                //在hb_MNA_Y中加入电压源支路,对所有频率遍历，加入KVL方程
                for(int h = -hb_params.num_harmonics; h <= hb_params.num_harmonics; ++h){
                    int row_index = (h + hb_params.num_harmonics) * base_size + branch_index;
                    if (n1 != 0){
                        hb_MNA_Y(row_index, (h + hb_params.num_harmonics) * base_size + n1 - 1) = 1;
                        hb_MNA_Y((h + hb_params.num_harmonics) * base_size + n1 - 1, row_index) = 1;
                        hb_jacobian(row_index, (h + hb_params.num_harmonics) * base_size + n1 - 1) = 1;
                        hb_jacobian((h + hb_params.num_harmonics) * base_size + n1 - 1, row_index) = 1;
                    }
                    if (n2 != 0){
                        hb_MNA_Y(row_index, (h + hb_params.num_harmonics) * base_size + n2 - 1) = -1;
                        hb_MNA_Y((h + hb_params.num_harmonics) * base_size + n2 - 1, row_index) = -1;
                        hb_jacobian(row_index, (h + hb_params.num_harmonics) * base_size + n2 - 1) = -1;
                        hb_jacobian((h + hb_params.num_harmonics) * base_size + n2 - 1, row_index) = -1;
                    }
                    //在J向量中加入电压源值，只在直流点加入，其余为0
                    if(h == 0){
                        hb_J(row_index) = value;
                    }
                    else{
                        hb_J(row_index) = 0;
                    }
                }
            }
        }
        if (c == 'I'){
            //独立电流源
            int n1 = dev.nodes[0];
            int n2 = dev.nodes[1];
            double value = dev.parameters["DC"];
            //在hb_J中加入电流源贡献，只在直流点加入
            hb_J(hb_params.num_harmonics * base_size + n1 - 1) -= value;
            hb_J(hb_params.num_harmonics * base_size + n2 - 1) += value;
        }
    }
}

//非线性器件的HB贡献
void solver::hb_build_nonlinear_MNA(){
    // 从上一次的频域解 hb_xw 计算时域解 hb_xt。
    // 注意：不要在这里重置 hb_xw（它在 PSS 入口处初始化一次），
    // 否则会改变 hb_xw 的尺寸导致后续比较/收敛判断出错。
    hb_xt = hb_iDFT(hb_xw);
    //初始化时域雅可比矩阵
    t_jacobian = Eigen::MatrixXd::Zero(base_size * (2 * hb_params.num_harmonics + 1), base_size * (2 * hb_params.num_harmonics + 1));
    //初始化时域J向量和频域J向量
    hb_J = Eigen::VectorXcd::Zero(base_size * (2 * hb_params.num_harmonics + 1));
    //对非线性器件进行多频率点的贡献构建
    Eigen::VectorXd I_dev_time = Eigen::VectorXd::Zero(base_size * (2 * hb_params.num_harmonics + 1));
    for (const auto &dev : ckt.nonlinear_devices) {
        char c = toupper(dev.name[0]);
        //初始化电流向量
        if(c == 'M'){
            // //Debug: 输出器件信息
            // std::cout << "Processing Nonlinear Device: " << dev.name << " Type: " << dev.type << "\n";
            //读取器件参数，不考虑body节点
            int n1 = dev.nodes[0];
            int ng0 = dev.nodes[1];
            int n2 = dev.nodes[2];
            double W = dev.parameters.at("W");
            double L = dev.parameters.at("L");
            double type = dev.parameters.at("TYPE"); //1 for NMOS, -1 for PMOS
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
            //提取时域节点电压
            Eigen::VectorXcd V1 = Eigen::VectorXcd::Zero(2 * hb_params.num_harmonics + 1);
            Eigen::VectorXcd Vg0 = Eigen::VectorXcd::Zero(2 * hb_params.num_harmonics + 1);
            Eigen::VectorXcd V2 = Eigen::VectorXcd::Zero(2 * hb_params.num_harmonics + 1);
            // //debug: 输出节点信息
            // std::cout << "Start trans" << "\n";
            // 从 hb_xt 中提取各节点的时域序列时要注意：
            // - 节点号为 0 表示接地，应直接赋 0
            // - hb_xt 长度为 base_size * N (N = 2*num_harmonics+1)，索引必须在 [0, hb_xt.size()-1] 范围内
            int N = 2 * hb_params.num_harmonics + 1;
            for(int n = 0; n < N; ++n){
                // 节点 n1
                if (n1 != 0) {
                    int idx1 = (n1 - 1) + n * base_size;
                    if (idx1 >= 0 && idx1 < hb_xt.size()) {
                        V1(n) = hb_xt(idx1);
                    } else {
                        V1(n) = std::complex<double>(0.0, 0.0);
                    }
                } else {
                    V1(n) = std::complex<double>(0.0, 0.0);
                }

                // 栅极节点 ng0
                if (ng0 != 0) {
                    int idxg = (ng0 - 1) + n * base_size;
                    if (idxg >= 0 && idxg < hb_xt.size()) {
                        Vg0(n) = hb_xt(idxg);
                    } else {
                        Vg0(n) = std::complex<double>(0.0, 0.0);
                    }
                } else {
                    Vg0(n) = std::complex<double>(0.0, 0.0);
                }

                // 节点 n2
                if (n2 != 0) {
                    int idx2 = (n2 - 1) + n * base_size;
                    if (idx2 >= 0 && idx2 < hb_xt.size()) {
                        V2(n) = hb_xt(idx2);
                    } else {
                        V2(n) = std::complex<double>(0.0, 0.0);
                    }
                } else {
                    V2(n) = std::complex<double>(0.0, 0.0);
                }
            }
            // //取实部
            // std::cout << "Start real part extraction" << "\n";
            Eigen::VectorXd V1_real = V1.real();
            Eigen::VectorXd Vg0_real = Vg0.real();
            Eigen::VectorXd V2_real = V2.real();
            //根据电压高低确定 Drain 和 Source, 注意NMOS和PMOS的源漏定义
            //初始化漏极电流 Ids
            Eigen::VectorXd Ids_time = Eigen::VectorXd::Zero(2 * hb_params.num_harmonics + 1);
            int nd, ns;
            double Vd0, Vs0;
            //遍历所有时域点
            for(int t = 0; t < (2 * hb_params.num_harmonics + 1); ++t){
                // //Debug: 输出时域点信息
                // std::cout << "  Time Point " << t << "\n";
                if (V1_real(t) > V2_real(t) && type > 0 || V1_real(t) <= V2_real(t) && type < 0) { 
                    nd = n1;
                    Vd0 = V1_real(t);
                    ns = n2; 
                    Vs0 = V2_real(t);
                } else {
                    nd = n2; 
                    Vd0 = V2_real(t);
                    ns = n1; 
                    Vs0 = V1_real(t);
                }

                double Vgs = type * (Vg0_real(t) - Vs0);
                double Vds = type * (Vd0 - Vs0);
                double Vth = VT * type; 
                //计算漏极电流 Ids
                // 计算 Id0, gm, gds, Ieq
                double Id0 = 0.0;
                double gm = 0.0;
                double gds = 0.0;
                double Ieq = 0.0;

                if (Vgs <= Vth) {
                    // cutoff
                    continue;
                } 
                else {
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
                        Id0 = 0.5 * beta * Vov * Vov * (1.0 + LAMBDA * Vds);
                        // gm = ∂Id/∂Vg = beta * Vov * (1 + lambda*Vds)
                        gm = beta * Vov * (1.0 + LAMBDA * Vds);
                        // gds = ∂Id/∂Vd = 0.5 * lambda * beta * Vov^2
                        gds = 0.5 * LAMBDA * beta * Vov;
                    }
                    //规定Id流出源极
                    Id0 = type * Id0;
                    // Ids_time(t) = Id0;
                    //把电流加入到时域电流向量中
                    // I_dev_time(t * base_size + nd - 1) -= Id0;
                    // I_dev_time(t * base_size + ns - 1) += Id0;
                    //gm和gds添加到时域雅可比矩阵中
                    //对漏极节点加入贡献
                    // //debug
                    // std::cout << "add Jacobian contributions at time point " << t << "\n";
                    if (nd != 0) {
                        int row_index = t * base_size + nd - 1;
                        //电流
                        I_dev_time(row_index) -= Id0;
                        // 对漏极节点的电压导数贡献
                        t_jacobian(row_index, t * base_size + nd - 1) += gds;
                        // 对栅极节点的电压导数贡献
                        if (ng0 != 0) {
                            t_jacobian(row_index, t * base_size + ng0 - 1) += gm;
                        }
                        // 对源极节点的电压导数贡献
                        if (ns != 0) {
                            t_jacobian(row_index, t * base_size + ns - 1) += - (gm + gds);
                        }
                    }
                    //对源极节点加入贡献
                    if (ns != 0) {
                        int row_index = t * base_size + ns - 1;
                        //电流
                        I_dev_time(row_index) += Id0;
                        // 对漏极节点的电压导数贡献
                        if(nd != 0)
                        t_jacobian(row_index, t * base_size + nd - 1) += -gds;
                        // 对栅极节点的电压导数贡献
                        if (ng0 != 0) {
                            t_jacobian(row_index, t * base_size + ng0 - 1) += -gm;
                        }
                        // 对源极节点的电压导数贡献
                        if (ns != 0) {
                            t_jacobian(row_index, t * base_size + ns - 1) += (gm + gds);
                        }
                    }
                }
                //Debug:
                // std::cout << "Time point " << t << ": Vd0=" << Vd0 << ", Vs0=" << Vs0 << ", Vgs=" << Vgs << ", Vds=" << Vds << ", Ids=" << Ids_time(t) << "\n";
            }
            
        }
        //  //Debug:
        // std::cout << "MOS Device " << dev.name << " processed in HB.\n";

    }

    // //Debug:输出时域电流向量
    // std::cout << "Nonlinear Device Time-Domain Current Vector (I_dev_time):\n" << I_dev_time << "\n";

    //将时域电流向量进行DFT变换，得到频域分量
    Eigen::VectorXcd I_dev_freq = hb_DFT(I_dev_time.cast<std::complex<double> >());

    // //Debug:输出频域电流向量
    // std::cout << "Nonlinear Device Frequency-Domain Current Vector (I_dev_freq):\n" << I_dev_freq << "\n";

    //将I_dev_freq加入到hb_J中
    for(int k = 0; k < (2 * hb_params.num_harmonics + 1); ++k){
        int row_index;
        //对各个节点加入电流源贡献
        for(int i = 0; i < base_size; ++i){
            row_index = (k) * base_size + i;
            hb_J(row_index) += I_dev_freq(row_index);
        }
        // //Debug:
        // std::cout << "Frequency point " << k - hb_params.num_harmonics << ": Current Contribution added to hb_J.\n";
    }

    // //线性元件对于频域雅可比矩阵的贡献
    // hb_jacobian = hb_T2F_matrix * t_jacobian * hb_F2T_matrix;

    //分块来计算矩阵乘法
    //分块法
    // 变换与矩阵乘法（按块实现）用时开始
    auto start_transform = std::chrono::high_resolution_clock::now();

    // 说明：hb_T2F_matrix 和 hb_F2T_matrix 是由每个变量重复的 DFT/IDFT 小块构成。
    // 可以对 t_jacobian 按变量对 (i,j) 提取 NxN 小块 B，计算 hb_DFT_matrix * B * hb_iDFT_matrix，
    // 再把结果写回 hb_jacobian 的相应位置。这样复杂度从 O((base_size*N)^3) 降为 O(base_size^2 * N^3)。

    int N = 2 * hb_params.num_harmonics + 1;
    int total_size = base_size * N;
    hb_jacobian = Eigen::MatrixXcd::Zero(total_size, total_size);

    // 预分配临时矩阵以减少重复分配开销
    Eigen::MatrixXcd block(N, N);
    Eigen::MatrixXcd trans(N, N);

    // 如果希望启用并行化，请在编译时加入 -fopenmp，并取消下面的注释（确保 Eigen 支持线程安全写入）
    // #pragma omp parallel for collapse(2) private(block, trans)
    for(int i = 0; i < base_size; ++i){
        for(int j = 0; j < base_size; ++j){
            // 提取 N x N 小块：block(k,n) = t_jacobian(i + k*base_size, j + n*base_size)
            for(int k = 0; k < N; ++k){
                for(int n = 0; n < N; ++n){
                    block(k, n) = t_jacobian(i + k * base_size, j + n * base_size);
                }
            }

            // 进行小块的频域变换
            trans.noalias() = hb_DFT_matrix * block * hb_iDFT_matrix;

            // 写回结果到 hb_jacobian
            for(int k = 0; k < N; ++k){
                for(int n = 0; n < N; ++n){
                    hb_jacobian(i + k * base_size, j + n * base_size) = trans(k, n);
                }
            }
        }
    }

    // 变换与矩阵乘法用时结束
    auto end_transform = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> transform_duration = end_transform - start_transform;
    std::cout << "Time-Frequency Transform Time: " << transform_duration.count() << " seconds.\n";


    // //Debug: 输出时域雅可比矩阵
    // std::cout << "Time-Domain Jacobian Matrix (t_jacobian):\n" << t_jacobian << "\n";
    // //Debug: 输出频域雅可比矩阵
    // std::cout << "MOS Frequency-Domain Jacobian Matrix (hb_jacobian):\n" << hb_jacobian << "\n";

    hb_jacobian += hb_liner_Y;
}


void solver::hb_solve_linear_MNA(){
    //求解多频率点的线性MNA方程
    Eigen::VectorXcd hb_x = hb_MNA_Y.fullPivLu().solve(hb_J);
    // //Debug: 输出多频率点的解向量
    // std::cout << "Harmonic Balance Linear MNA Solution (x):\n" << hb_x << "\n";
    //进行IDFT变换，得到时域解
    // Eigen::VectorXcd hb_xt = hb_iDFT(hb_x);
    hb_xt = hb_F2T_matrix * hb_x;
    // //Debug: 输出时域解向量
    // std::cout << "Harmonic Balance Time-Domain Solution (xt):\n" << hb_xt << "\n";
}

//设置初始频域解
void solver::HB_set_initial_xw(const std::map<std::string, double>& node_voltage_map){
    //运行一篇直流分析，得到base_size
    //先对直流点构建线性MNA矩阵
    build_linear_MNA(false);
    MNA_Y = Eigen::MatrixXd::Zero(liner_Y.rows(), liner_Y.cols());
    MNA_Y = liner_Y;
    build_sources_MNA();
    //这样就得到的MNA_Y大小就是每个频率点的矩阵大小
    base_size = MNA_Y.rows();

    int N = 2 * hb_params.num_harmonics + 1;
    hb_xw = Eigen::VectorXcd::Zero(base_size * N);
    //遍历节点电压映射，设置初始频域解
    for(const auto& pair : node_voltage_map){
        const std::string& node_name = pair.first;
        double voltage = pair.second;
        //找到节点索引
        int node_index = ckt.getNodeID(node_name);
        if(node_index == -1){
            std::cerr << "Warning: Node " << node_name << " not found in circuit. Skipping initial condition setting.\n";
            continue;
        }
        //设置基频分量
        hb_xw(node_index - 1 + (hb_params.num_harmonics - 1) * base_size) = std::complex<double>(voltage, 0.0);
        hb_xw(node_index - 1 + (hb_params.num_harmonics + 1) * base_size) = std::complex<double>(voltage, 0.0);
    }

    // //Debug: 输出初始频域解
    // std::cout << "Initial Harmonic Balance Frequency-Domain Solution (xw):\n" << hb_xw << "\n";
}

void solver::PSS_solve_harmonic_balance(){

    //构建多频率点的线性MNA矩阵
    hb_build_linear_MNA();
    // //Debug: 输出多频率点的线性MNA矩阵
    // std::cout << "Harmonic Balance Linear MNA Matrix (Y):\n" << hb_liner_Y << "\n";

    //初始化DFT和IDFT矩阵
    hb_initialize_DFT_matrices();
    hb_build_TF_matrix();

    // //确定需要打印的节点
    // parse_print_variables();
    // 初始化频域解向量
    hb_xt = Eigen::VectorXcd::Zero(base_size * (2 * hb_params.num_harmonics + 1));
    hb_xt = hb_iDFT(hb_xw);

    //进行迭代求解
    Eigen::VectorXcd hb_xw_old = Eigen::VectorXcd::Zero(base_size * (2 * hb_params.num_harmonics + 1));
    for(int iter = 0; iter < hb_params.max_iterations; ++iter){
        
        //每次迭代用时
        auto start_iter = std::chrono::high_resolution_clock::now();

        std::cout << "Harmonic Balance Iteration " << iter + 1 << ":\n";
        // //Debug: 输出当前频域解
        // std::cout << "Current Frequency-Domain Solution (xw):\n" << hb_xw << "\n";
        // //Debug: 输出时域解
        // std::cout << "Current Time-Domain Solution (xt):\n" << hb_xt << "\n";

        //保存上一次的频域解
        hb_xw_old = hb_xw;
        //构建非线性器件的HB贡献

        //对构建的时间进行计时
        auto start_nonlinear = std::chrono::high_resolution_clock::now();
        
        hb_build_nonlinear_MNA();

        // 结束时间点
        auto end_nonlinear = std::chrono::high_resolution_clock::now();
        auto nonlinear_seconds = std::chrono::duration<double>(end_nonlinear - start_nonlinear).count();
        std::cout << "构建非线性器件贡献耗时: " << nonlinear_seconds << " 秒" << std::endl;

        // std::cout << "Nonlinear MNA contribution built.\n";
        //加入sources到多频率点MNA矩阵和J向量中

        //对构建的时间进行计时
        auto start_sources = std::chrono::high_resolution_clock::now();

        hb_build_sources_MNA();

        // 结束时间点
        auto end_sources = std::chrono::high_resolution_clock::now();
        auto sources_seconds = std::chrono::duration<double>(end_sources - start_sources).count();
        std::cout << "构建sources贡献耗时: " << sources_seconds << " 秒" << std::endl;
        // std::cout << "Sources contribution built.\n";
        //直接求解法求解多频率点的线性MNA方程
        // hb_solve_linear_MNA();
        //已经得到新的频域解
        //初值调节迭代



        // //Debug: 输出当前MNA矩阵和J向量
        // std::cout << "Harmonic Balance MNA Matrix (Y):\n" << hb_MNA_Y << "\n";
        // std::cout << "Harmonic Balance Jacobian Matrix:\n" << hb_jacobian << "\n";
        // std::cout << "Harmonic Balance J Vector:\n" << hb_J << "\n";

        //迭代求解
        Eigen::VectorXcd delta_F = hb_J - (hb_MNA_Y * hb_xw);

    //Debug: 计时开始
    //开始时间点
    auto start_solveMNA = std::chrono::high_resolution_clock::now();
    
    // Eigen::VectorXcd delta_xw = hb_jacobian.fullPivLu().solve(delta_F);
    //使用Eigen的LU分解求解器
    Eigen::PartialPivLU<Eigen::MatrixXcd> lu(hb_jacobian);
    
    //求解线性方程组 MNA_Y * x = J
    Eigen::VectorXcd delta_xw = lu.solve(delta_F);

    // Debug 或者直接获取秒
    // 结束时间点
    auto end_solveMNA = std::chrono::high_resolution_clock::now();
    auto seconds = std::chrono::duration<double>(end_solveMNA - start_solveMNA).count();
    std::cout << "求解jacobian矩阵耗时: " << seconds << " 秒" << std::endl;

    // //Debug:展示增量
    // std::cout << "Delta_xw:\n" << delta_xw << "\n";


        hb_xw += delta_xw;
        
        //检查收敛性
        double norm = (hb_xw - hb_xw_old).norm();
        //Debug: 输出收敛性指标
        std::cout << "Convergence Norm: " << norm << "\n";

        if (norm < hb_params.tolerance) {
            std::cout << "Converged after " << iter + 1 << " iterations.\n";
            break;
        }
        hb_xt = hb_iDFT(hb_xw);

        //每次迭代用时
        auto end_iter = std::chrono::high_resolution_clock::now();
        auto iter_seconds = std::chrono::duration<double>(end_iter - start_iter).count();
        std::cout << "Iteration " << iter + 1 << " completed in " << iter_seconds << " seconds.\n";

    }
    //最终求解得到时域解
    hb_xt = hb_iDFT(hb_xw);
    //给出稳态条件
    node_voltages = Eigen::VectorXd(ckt.node_list.size() -1);
    for(int i = 0; i < (ckt.node_list.size() -1); ++i){
        node_voltages(i) = hb_xt(i).real();
        //std::cout << "Steady-State Voltage at Node " << ckt.node_list[i +1] << ": " << node_voltages(i) << " V\n";
    }
    // //Debug: 输出DFT和IDFT矩阵
    // std::cout << "DFT Matrix:\n" << hb_DFT_matrix << "\n";
    // std::cout << "IDFT Matrix:\n" << hb_iDFT_matrix << "\n";
    // // Debug：自己构建一个频域解向量，测试IDFT变换
    // Eigen::VectorXcd test_xw = Eigen::VectorXcd::Zero(base_size * (2*hb_params.num_harmonics+1));
    // int i = 0;
    //     for(int k = 0; k < (2*hb_params.num_harmonics+1); ++k){
    //         if(k == 0 || k == 2*hb_params.num_harmonics)
    //             test_xw(i + k * base_size) = std::complex<double>(0.5, 0.0); //最高频点余弦波 
    //         else
    //         test_xw(i + k * base_size) = std::complex<double>(0.0, 0.0); //简单测试值
    //     }
    //     i = 1;
    //     for(int k = 0; k < (2*hb_params.num_harmonics+1); ++k){
    //         if(k == hb_params.num_harmonics)
    //             test_xw(i + k * base_size) = std::complex<double>(1.0, 0.0); //直流分量
    //         else
    //         test_xw(i + k * base_size) = std::complex<double>(0.0, 0.0); //简单测试值
    //     }
    //     std::cout << "Test DFT Input:\n" << test_xw << "\n";
    // Eigen::VectorXcd test_xt = hb_iDFT(test_xw);
    // std::cout << "Test IDFT Result:\n" << test_xt << "\n";
    // //测试DFT变换
    // Eigen::VectorXcd test_xw_back = hb_DFT(test_xt);
    // std::cout << "Test DFT Result:\n" << test_xw_back << "\n";


    // //Debug: 输出多频率点的线性MNA矩阵和J向量
    // //输出到文件中查看
    // std::ofstream file_Y("hb_MNA_Y.txt");
    // if (file_Y.is_open()) {
    //     file_Y << hb_MNA_Y << std::endl;
    //     file_Y << hb_J << std::endl;
    //     file_Y.close();
    // } else {
    //     std::cerr << "无法打开文件 hb_MNA_Y.txt 进行写入。" << std::endl;
    // }

    //打印需要打印的节点
        //根据需要打印的变量，存到文件中
        {
            // 输出文件: hb_print.txt
                std::ofstream hdr("hb_print.txt", std::ios::out);
                hdr << "Time(s)";
                for (int node_id : ckt.print_node_ids) {
                    std::string name = "NODE";
                    if (node_id >= 0 && node_id < (int)ckt.node_list.size()) name = ckt.node_list[node_id];
                    hdr << "\tV(" << name << ")";
                }
                // //只需要遍历所有sources，按顺序输出支路电流表头
                for (const auto &d : ckt.sources){
                    if (d.printI) hdr << "\tI(" << d.name << ")";
                }
                //关闭
                hdr << "\n";
                hdr.close();
            

            std::ofstream out("hb_print.txt", std::ios::app);
            //遍历所有时域点，输出需要打印的节点电压和支路电流
            int N = 2 * hb_params.num_harmonics + 1;
            double T = 1.0 / (hb_params.fundamental_omega / (2.0 * M_PI)); //周期
            double time = 0.0;
            for(int n = 0; n < N; ++n){
                time = n * T / N;
                //提取节点电压
                Eigen::VectorXd hb_node_voltages(ckt.node_list.size() -1);
                for(int i = 0; i < (ckt.node_list.size() -1); ++i){
                    hb_node_voltages(i) = hb_xt(i + n * base_size).real();
                }
                //1215提取支路电流
                Eigen::VectorXd branch_currents(ckt.sources.size());
                for(int i = 0; i < ckt.sources.size(); ++i){
                    int branch_index = ckt.sources[i].branch_current_index + (ckt.node_list.size() -1);
                    branch_currents(i) = hb_xt(branch_index + n * base_size).real();
                }
                out << time;
                for (int node_id : ckt.print_node_ids) {
                    double v = 0.0;
                    if (node_id == 0) v = 0.0;
                    else if (node_id - 1 >= 0 && node_id - 1 < hb_node_voltages.size()) v = hb_node_voltages[node_id - 1];
                    out << "\t" << v;
                }
                out << "\n";

            //1215 电流
            for (int current_dev_index : ckt.print_branch_current_indices) {
                if(current_dev_index >=0 && current_dev_index < ckt.sources.size()){
                    out << "\t" << branch_currents[current_dev_index];
                }
            }
            }

            //打印频域结果
            out << "\nFrequency Domain Results:\n";
            out << "Harmonic\tFrequency(Hz)";
            for (int node_id : ckt.print_node_ids) {
                std::string name = "NODE";
                if (node_id >= 0 && node_id < (int)ckt.node_list.size()) name = ckt.node_list[node_id];
                out << "\tV(" << name << ")";
            }
            out << "\n";
            for (int h = 0; h < N; ++h) {
                Eigen::VectorXcd hb_node_vw(ckt.node_list.size() -1);
                for(int i = 0; i < (ckt.node_list.size() -1); ++i){
                    hb_node_vw(i) = hb_xw(i + h * base_size);
                }
                out << h - hb_params.num_harmonics << "\t" << ((h - hb_params.num_harmonics) * (hb_params.fundamental_omega / (2.0 * M_PI)));
                for (int node_id : ckt.print_node_ids) {
                    std::complex<double> v = 0.0;
                    if (node_id == 0) v = 0.0;
                    else if (node_id - 1 >= 0 && node_id - 1 < hb_xw.size()) v = hb_node_vw[node_id - 1];
                    out << "\t" << v;
                }
                out << "\n";
            }

            out << "\n";
            out.close();
        }
}

void solver::PSS_solve_harmonic_balance(analysis& analysis, int max_iters, double tol, double relaxation_factor){
    //根据网表设定参数
    hb_params.fundamental_omega = 2 * 3.14159265358979323846 * analysis.parameters["freq"]; // 基频角频率
    hb_params.num_harmonics = static_cast<int>(analysis.parameters["harm"]); // 谐波数量
    hb_params.max_iterations = max_iters; // 最大迭代次数
    hb_params.tolerance = tol; // 收敛容限
    hb_params.relaxation_factor = relaxation_factor; // 松弛因子
    ckt.extract_MOS_capacitances(); //提取MOS管寄生电容
        //确定需要打印的节点
    parse_print_variables();
    HB_set_initial_xw({}); //空参数表示使用默认初始解
    PSS_solve_harmonic_balance();
    std::cout << "Harmonic Balance Analysis Completed.\n";
}

void solver::print_hb_time_domain_results(){
    //输出到文件中查看,把std::cout改为文件输出
    std::ofstream out("hb_time_domain_results.txt");
    out << "Time(s)";
    for (int node_id : ckt.print_node_ids) {
        std::string name = "NODE";
        if (node_id >= 0 && node_id < (int)ckt.node_list.size()) name = ckt.node_list[node_id];
        out << "\tV(" << name << ")";
    }
    for(const auto &d : ckt.sources){
        if (d.printI) out << "\tI(" << d.name << ")";
    }
    out << "\n";

    int N = 2 * hb_params.num_harmonics + 1;
    double T = 1.0 / (hb_params.fundamental_omega / (2.0 * M_PI)); //周期
    double time = 0.0;
    for(int n = 0; n < N; ++n){
        time = n * T / N;
        //提取节点电压

        Eigen::VectorXd hb_node_voltages(ckt.node_list.size() -1);
        for(int i = 0; i < (ckt.node_list.size() - 1); ++i){
            hb_node_voltages(i) = hb_xt(i + n * base_size).real();
        }
        //提取支路电流
        Eigen::VectorXd branch_currents(ckt.sources.size());
        for(int i = 0; i < ckt.sources.size(); ++i){
            int branch_index = ckt.sources[i].branch_current_index + (ckt.node_list.size() -1);
            branch_currents(i) = hb_xt(branch_index + n * base_size).real();
        }
        out << time;
        for (int node_id : ckt.print_node_ids) {
            double v = 0.0;
            if (node_id == 0) v = 0.0;
            else if (node_id - 1 >= 0 && node_id - 1 < hb_node_voltages.size()) v = hb_node_voltages[node_id - 1];
            out << "\t" << v;
        }
        //打印支路电流
        for (int current_dev_index : ckt.print_branch_current_indices) {
            // std::cout << current_dev_index << "\n";
            if(current_dev_index >=0 && current_dev_index < ckt.sources.size()){
                out << "\t" << branch_currents[current_dev_index];
            }
        }
        out << "\n";
    }
    out.close();
}

void solver::print_hb_frequency_domain_results(){
    //打印频域结果到文件中
    int N = 2 * hb_params.num_harmonics + 1;
    std::ofstream out("hb_frequency_domain_results.txt");
    out << "\nFrequency Domain Results:\n";
    out << "Harmonic\tFrequency(Hz)";
    for (int node_id : ckt.print_node_ids) {
        std::string name = "NODE";
        if (node_id >= 0 && node_id < (int)ckt.node_list.size()) name = ckt.node_list[node_id];
        out << "\tV(" << name << ")";
    }
    for(const auto &d : ckt.sources){
        if (d.printI) out << "\tI(" << d.name << ")";
    }
    out << "\n";
    for (int h = 0; h < N; ++h) {
        Eigen::VectorXcd hb_node_vw(base_size);
        for(int i = 0; i < base_size; ++i){
            hb_node_vw(i) = hb_xw(i + h * base_size);
        }
        //提取支路电流
        Eigen::VectorXcd branch_currents_w(ckt.sources.size());
        for(int i = 0; i < ckt.sources.size(); ++i){
            int branch_index = ckt.sources[i].branch_current_index + (ckt.node_list.size() -1);
            branch_currents_w(i) = hb_xw(branch_index + h * base_size);
        }
        out << h - hb_params.num_harmonics << "\t" << ((h - hb_params.num_harmonics) * (hb_params.fundamental_omega / (2.0 * M_PI)));
        for (int node_id : ckt.print_node_ids) {
            std::complex<double> v = 0.0;
            if (node_id == 0) v = 0.0;
            else if (node_id - 1 >= 0 && node_id - 1 < hb_node_vw.size()) v = hb_node_vw[node_id - 1];
            out << "\t" << v;
        }
        // //打印支路电流
        // for (int current_dev_index : ckt.print_branch_current_indices) {
        //     if(current_dev_index >=0 && current_dev_index < ckt.sources.size()){
        //         std::complex<double> i = hb_node_vw[current_dev_index + (ckt.node_list.size() -1)];
        //         out << "\t" << i;
        //     }
        // }
        //打印支路电流
        for (int current_dev_index : ckt.print_branch_current_indices) {
            // //Debug:
            // std::cout << current_dev_index << "\n";
            if(current_dev_index >=0 && current_dev_index < ckt.sources.size()){
                out << "\t" << branch_currents_w[current_dev_index];
            }
        }
        out << "\n";
    }
    out.close();
}

// void solver::plot_hb_time_domain_results() {
//     // 时间采样点数
//     int N = 2 * hb_params.num_harmonics + 1;
//     double T = 1.0 / (hb_params.fundamental_omega / (2.0 * M_PI)); //周期

//     // 构造时间轴
//     std::vector<double> time_points;
//     for (int n = 0; n < N; ++n) {
//         time_points.push_back(n * T / N);
//     }

//     // 创建一个图
//     plt::figure();

//     for (int node_id : get_plot_node_ids()) {
//         std::string name = "NODE";
//         if (node_id >= 0 && node_id < (int)ckt.node_list.size())
//             name = ckt.node_list[node_id];
//         // 提取节点电压
//         std::vector<double> voltages;
//         for (int n = 0; n < N; ++n) {
//             voltages.push_back(hb_xt((node_id - 1) + n * base_size).real());
//             std::cout << time_points[n] << "\t" << voltages[n] << "\n";
//         }
//         // 绘制节点电压波形
//         plt::plot(time_points, voltages, {{"label", "Node " + name}});
//         std::cout << "Plotted Node " << name << " voltage.\n";
//     }
//     plt::legend();
//     plt::title("HB Time-Domain Node Voltages");
//     plt::xlabel("Time (s)");
//     plt::ylabel("Voltage (V)");
//     plt::grid(true);
//     plt::show();

//     std::cout << "Plotted HB time-domain voltages for all nodes.\n";
// }



//神秘实验内容：对角线最优化
    // 找到最大匹配的BFS函数（用于匈牙利算法）
    bool solver::bfs_match(const std::vector<std::vector<int>>& graph, 
                   std::vector<int>& pairU, 
                   std::vector<int>& pairV, 
                   std::vector<int>& dist, 
                   int N) {
        std::vector<int> q;
        for (int u = 0; u < N; u++) {
            if (pairU[u] == -1) {
                dist[u] = 0;
                q.push_back(u);
            } else {
                dist[u] = -1;
            }
        }
        
        bool found = false;
        int qpos = 0;
        
        while (qpos < q.size()) {
            int u = q[qpos++];
            for (int v : graph[u]) {
                if (pairV[v] == -1) {
                    found = true;
                } else if (dist[pairV[v]] == -1) {
                    dist[pairV[v]] = dist[u] + 1;
                    q.push_back(pairV[v]);
                }
            }
        }
        
        return found;
    }

    // DFS寻找增广路
    bool solver::dfs_match(int u, const std::vector<std::vector<int>>& graph,
                   std::vector<int>& pairU, std::vector<int>& pairV,
                   std::vector<int>& dist) {
        for (int v : graph[u]) {
            if (pairV[v] == -1 || (dist[pairV[v]] == dist[u] + 1 && 
                dfs_match(pairV[v], graph, pairU, pairV, dist))) {
                pairU[u] = v;
                pairV[v] = u;
                return true;
            }
        }
        dist[u] = -1;
        return false;
    }

    // 匈牙利算法找到完美匹配
    std::vector<int> solver::hungarianMatch(const Eigen::MatrixXd& matrix) {
        int n = matrix.rows();
        std::vector<std::vector<int>> graph(n);
        
        // 构建二分图：对于每行，选择绝对值最大的k个元素
        // 这里选择前min(3, n)个最大元素，保证匹配存在
        for (int i = 0; i < n; i++) {
            std::vector<std::pair<double, int>> elements;
            for (int j = 0; j < n; j++) {
                elements.emplace_back(std::abs(matrix(i, j)), j);
            }
            // 按绝对值降序排序
            std::sort(elements.begin(), elements.end(), 
                     [](const auto& a, const auto& b) {
                         return a.first > b.first;
                     });
            
            // 取前k个作为候选边
            int k = std::min(3, n);
            for (int idx = 0; idx < k; idx++) {
                graph[i].push_back(elements[idx].second);
            }
        }
        
        std::vector<int> pairU(n, -1);  // 行到列的匹配
        std::vector<int> pairV(n, -1);  // 列到行的匹配
        std::vector<int> dist(n);
        
        int matching = 0;
        while (bfs_match(graph, pairU, pairV, dist, n)) {
            for (int u = 0; u < n; u++) {
                if (pairU[u] == -1 && dfs_match(u, graph, pairU, pairV, dist)) {
                    matching++;
                }
            }
        }
        
        // 如果找到完美匹配，返回列置换
        if (matching == n) {
            return pairU;  // pairU[i] 表示第i行匹配的列
        }
        
        // 如果没找到完美匹配，使用简单的贪心策略
        std::vector<int> colMatch(n, -1);
        std::vector<int> rowMatch(n, -1);
        
        for (int i = 0; i < n; i++) {
            // 找到当前行中绝对值最大的可用列
            int bestCol = -1;
            double maxVal = -1;
            
            for (int j = 0; j < n; j++) {
                if (colMatch[j] == -1) {
                    double val = std::abs(matrix(i, j));
                    if (val > maxVal) {
                        maxVal = val;
                        bestCol = j;
                    }
                }
            }
            
            if (bestCol != -1) {
                colMatch[bestCol] = i;
                rowMatch[i] = bestCol;
            }
        }
        
        // 处理未匹配的行
        for (int i = 0; i < n; i++) {
            if (rowMatch[i] == -1) {
                // 找到未使用的列
                for (int j = 0; j < n; j++) {
                    if (colMatch[j] == -1) {
                        rowMatch[i] = j;
                        colMatch[j] = i;
                        break;
                    }
                }
            }
        }
        
        return rowMatch;
    }


    Eigen::MatrixXd solver::diagMax() {
        int n = MNA_Y.rows();
        
        // 1. 找到列置换（通过最大匹配）
        std::vector<int> colPerm = hungarianMatch(MNA_Y);
        
        // 2. 应用列置换到矩阵
        Eigen::MatrixXd Y_perm_cols = Eigen::MatrixXd::Zero(n, n);
        for (int i = 0; i < n; i++) {
            Y_perm_cols.col(i) = MNA_Y.col(colPerm[i]);
        }
        
        // 3. 应用行置换：将每行最大元素放到对角线
        std::vector<int> rowPerm(n);
        for (int i = 0; i < n; i++) {
            // 找到当前行绝对值最大的元素
            int maxIdx = 0;
            double maxVal = std::abs(Y_perm_cols(i, 0));
            for (int j = 1; j < n; j++) {
                double val = std::abs(Y_perm_cols(i, j));
                if (val > maxVal) {
                    maxVal = val;
                    maxIdx = j;
                }
            }
            // 将最大元素所在的行换到对应列的位置
            rowPerm[maxIdx] = i;
        }
        
        // 检查rowPerm是否是有效排列
        std::vector<bool> used(n, false);
        for (int i = 0; i < n; i++) {
            if (rowPerm[i] < 0 || rowPerm[i] >= n || used[rowPerm[i]]) {
                // 如果排列无效，使用顺序排列
                for (int j = 0; j < n; j++) rowPerm[j] = j;
                break;
            }
            used[rowPerm[i]] = true;
        }
        
        // 4. 应用行置换到矩阵和激励向量
        Eigen::MatrixXd Y_perm = Eigen::MatrixXd::Zero(n, n);
        Eigen::VectorXd J_perm(n);
        for (int i = 0; i < n; i++) {
            Y_perm.row(i) = Y_perm_cols.row(rowPerm[i]);
            J_perm(i) = J(rowPerm[i]);
        }
        
        // 5. 更新原始矩阵和激励向量
        MNA_Y = Y_perm;
        J = J_perm;
        
        // 6. 构建逆置换矩阵
        // 逆置换顺序是：先逆行置换，再逆列置换
        Eigen::MatrixXd inversePerm = Eigen::MatrixXd::Identity(n, n);
        
        // 逆列置换矩阵
        Eigen::MatrixXd colPermMat = Eigen::MatrixXd::Zero(n, n);
        for (int i = 0; i < n; i++) {
            colPermMat(i, colPerm[i]) = 1.0;
        }
        
        // 逆行置换矩阵
        Eigen::MatrixXd rowPermMat = Eigen::MatrixXd::Zero(n, n);
        for (int i = 0; i < n; i++) {
            rowPermMat(rowPerm[i], i) = 1.0;  // 注意：这里是逆置换
        }
        
        // 整体逆置换矩阵 = 逆列置换 * 逆行置换
        // 因为原始变换是：Y_new = P_row * Y_old * P_col
        // 所以恢复时需要：Y_old = P_row^T * Y_new * P_col^T
        inversePerm = colPermMat.transpose() * rowPermMat.transpose();
        
        return inversePerm;
    }


    /**
     * @brief 使用逆置换矩阵恢复解的顺序
     * @param solution 置换后系统求得的解
     * @param inversePerm 由diagMax()返回的逆置换矩阵
     * @return 恢复原始顺序的解
     */
    Eigen::VectorXd solver::restoreOrder(const Eigen::VectorXd& solution, 
                                        const Eigen::MatrixXd& inversePerm) {
        return inversePerm * solution;
    }
    
    /**
     * @brief 检查矩阵是否对角占优
     * @param threshold 严格对角占优的阈值
     * @return 对角占优的行数
     */
    int solver::checkDiagonalDominance(double threshold) const {
        int n = MNA_Y.rows();
        int count = 0;
        
        for (int i = 0; i < n; i++) {
            double diag = std::abs(MNA_Y(i, i));
            double sum = 0.0;
            
            for (int j = 0; j < n; j++) {
                if (j != i) {
                    sum += std::abs(MNA_Y(i, j));
                }
            }
            
            if (diag > sum + threshold) {
                count++;
            }
        }
        
        return count;
    }

// 运行一次瞬态，用给定初始条件 init_x，返回 T 处的节点电压
std::pair<Eigen::VectorXd, Eigen::VectorXd> solver::run_transient_once(double T, double tstep, const Eigen::VectorXd& init_V, const Eigen::VectorXd& init_I)
{
    // 设置初始节点电压
    node_voltages = init_V;
    branch_currents = init_I;
    // 计算步数
    int steps = static_cast<int>(T / tstep);

    tran_plot_data.clear();
    //std::cout << branch_currents[0] << "\n\n";
    //system("pause");
    for (int step = 0; step <= steps + 1; step++) {
        double time = step * tstep;

        // 构建瞬态分析电路
        build_transient_ckt(tstep);
        
        // 以当前 node_voltages 为初值求解非线性 MNA
        DC_solve(node_voltages, true, time);

        // node_voltages 会在 DC_solve 中被更新，无需额外处理
        //std::cout << node_voltages.size() << std::endl;

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
                //输出调试信息
                //std::cout << "Plot Data - Time: " << time << " s, Branch Index: " << plot_current_dev_index << ", Current: " << i << " A\n";
                //std::cout << time << "\t" << i << "\n";
            }
        }
    }
    //system("pause");
    return std::make_pair(node_voltages, branch_currents);   // 即 v(T) 和 i(T)
}

void solver::PSS_solve_shooting(double period_T, double tstep, int max_iters, double tol){
    //确定需要打印的变量
    parse_print_variables();

    //提取电容信息
    ckt.extract_MOS_capacitances();

    // 节点个数
    int N = ckt.node_map.size() - 1;

    // ---- Step 0：初始化初始条件 X0 ----
    // 你可以用 DC 解，或者直接用 0
    DC_solve();
    Eigen::VectorXd V0 = node_voltages; // 初始节点电压
    Eigen::VectorXd I0 = branch_currents; // 初始支路电流

    std::cout << "Shooting Method Start: N = " << N << "\n";
    int iter = 0;
    for (iter = 0; iter < max_iters; iter++)
    {
        std::cout << "=== Shooting Iteration " << iter << " ===\n";

        // ---- Step1：从 X0 出发运行瞬态，得到周期末 v(T) ----
        std::pair<Eigen::VectorXd, Eigen::VectorXd> XT = run_transient_once(period_T, tstep, V0, I0);

        Eigen::VectorXd VT = XT.first;
        Eigen::VectorXd IT = XT.second;

        // ---- Step2：误差 F = XT - X0 ----
        if (VT.size() != V0.size()) {
            V0 = Eigen::VectorXd::Zero(VT.size());
        }
        if (IT.size() != I0.size()) {
            I0 = Eigen::VectorXd::Zero(IT.size());
        }
        Eigen::VectorXd F_V = VT - V0;
        double err_V = F_V.norm();

        Eigen::VectorXd F_I = IT - I0;
        double err_I = F_I.norm();

        std::cout << "Error norm = " << err_V << " " << err_I << "\n";

        // ---- Step3：检查收敛 ----
        if (err_V < tol && err_I < tol) {
            std::cout << "Shooting method converged.\n";
            node_voltages = VT;    // 最终稳态
            branch_currents = IT;
            break;
        }

        // ---- Step4：更新初始条件 ----
        // 松弛法：X0 ← X0 + α*(XT - X0)
        double alpha = 0.5;   // 可调，0.3~0.8 之间效果较好
        V0 = V0 + alpha * (VT - V0); // 临时存储当前初始条件
        I0 = I0 + alpha * (IT - I0);
    }
    if (iter == max_iters){
        std::cout << "Warning: Shooting method did NOT converge within max_iters.\n";
        node_voltages = V0;
        branch_currents = I0;
    }
    // //展示节点电压结果
    // std::cout << "PSS Analysis Node Voltages:\n";
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


void solver::PSS_solve_shooting_exact_jacobian(double period_T, double tstep, int max_iters, double tol) {
    //确定需要打印的变量
    parse_print_variables();

    //提取电容信息
    ckt.extract_MOS_capacitances();
    // 节点个数
    int N = ckt.node_map.size() - 1;
    DC_solve();
    Eigen::VectorXd V0 = node_voltages; // 初始节点电压
    Eigen::VectorXd I0 = branch_currents; // 初始支路电流
    Eigen::VectorXd X0 = Eigen::VectorXd::Zero(V0.size() + I0.size());
    X0 << V0, I0;
    int n = X0.size();
    
    Eigen::MatrixXd J(n, n);
    Eigen::VectorXd F(n);
    double prev_error = 1e20;
    for (int iter = 0; iter < max_iters; iter++) {
        // Step 1: 计算 F = XT - X0
        std::pair<Eigen::VectorXd, Eigen::VectorXd> XT_pair = run_transient_once(period_T, tstep, V0, I0);
        Eigen::VectorXd VT = XT_pair.first;
        Eigen::VectorXd IT = XT_pair.second;
        Eigen::VectorXd XT(n);
        XT << VT, IT;
        F = XT - X0;
        
        double error = F.norm();
        std::cout << "Iter " << iter << ": error = " << error << std::endl;
        
        if (error < tol) {
            std::cout << "Converged!" << std::endl;
            node_voltages = VT;
            branch_currents = IT;
            return;
        }
        
        // Step 2: 计算精确雅可比（只在需要时）
        bool recompute_jacobian = true;
        if (iter > 0) {
            // 检查收敛速度，决定是否重新计算雅可比
            double error_reduction = error / prev_error;
            recompute_jacobian = (error_reduction > 0.5) || (iter % 3 == 0);
        }
        
        if (recompute_jacobian) {
            std::cout << "  Computing exact Jacobian..." << std::endl;
            // 对每个分量进行扰动
            for (int j = 0; j < n; j++) {
                Eigen::VectorXd X0_perturbed = X0;
                //把X0_perturbed拆成电压和电流部分
                Eigen::VectorXd V0_perturbed = X0_perturbed.head(V0.size());
                Eigen::VectorXd I0_perturbed = X0_perturbed.tail(I0.size());
                X0_perturbed(j) += 1e-10; // 小扰动
                
                std::pair<Eigen::VectorXd, Eigen::VectorXd> XT_perturbed_pair = run_transient_once(period_T, tstep, V0_perturbed, I0_perturbed);
                Eigen::VectorXd VT_perturbed = XT_perturbed_pair.first;
                Eigen::VectorXd IT_perturbed = XT_perturbed_pair.second;
                Eigen::VectorXd XT_perturbed(n);
                XT_perturbed << VT_perturbed, IT_perturbed;
                Eigen::VectorXd F_perturbed = XT_perturbed - X0_perturbed;
                
                // 计算列 j 的导数：dF/dX0_j
                J.col(j) = (F_perturbed - F) / 1e-10;
            }
        }
        
        // Step 3: 牛顿步
        Eigen::VectorXd delta_X = -J.lu().solve(F);
        
        // Step 4: 线搜索
        double lambda = 1.0;
        for (int ls = 0; ls < 10; ls++) {
            Eigen::VectorXd X_new = X0 + lambda * delta_X;
            Eigen::VectorXd V_new = X_new.head(V0.size());
            Eigen::VectorXd I_new = X_new.tail(I0.size());
            std::pair<Eigen::VectorXd, Eigen::VectorXd> XT_new_pair = run_transient_once(period_T, tstep, V_new, I_new);
            Eigen::VectorXd XT_new = XT_new_pair.first;
            Eigen::VectorXd F_new = XT_new - X_new;
            
            if (F_new.norm() < F.norm() * (1.0 - 0.1 * lambda)) {
                X0 = X_new;
                break;
            }
            lambda *= 0.5;
        }
        
        prev_error = error;
    }
    
    std::cout << "Warning: Not converged" << std::endl;
}

void solver::PSS_solve_shooting_new(double period_T, double tstep, int max_iters, double tol){
    //确定需要打印的变量
    parse_print_variables();

    //提取电容信息
    ckt.extract_MOS_capacitances();

    // ---- Step 0：初始化初始条件 X0 ----
    // 你可以用 DC 解，或者直接用 0
    DC_solve();
    int N;
    int iter = 0;
    for (iter = 0; iter < max_iters; iter++)
    {
        std::cout << "=== Shooting Iteration " << iter << " ===\n";

        Eigen::VectorXd V0 = node_voltages; // 初始节点电压

        //构建瞬态分析电路
        build_transient_ckt(tstep);

        // ---- Step1：从 X0 出发运行瞬态，得到周期末 v(T) ----
        //基准transient求解
        TRAN_solve_with_initial_value(period_T, tstep);

        Eigen::VectorXd VT = node_voltages;

        if (iter == 0){
            N = VT.size();
            V0 = Eigen::VectorXd::Zero(N);
        }

        // ---- Step2：误差 F = XT - X0 ----
        //std::cout << VT << "\n\n" << V0 << "\n";
        Eigen::VectorXd F = VT - V0;
        double err = F.norm();

        std::cout << "Error norm = " << err << "\n";

        // ---- Step3：检查收敛 ----
        if (err < tol) {
            std::cout << "Shooting method converged.\n";
            break;
        }

        // ---- Step4：更新初始条件 ----

        // // 松弛法：X0 ← X0 + α*(XT - X0)
        // double alpha = 0.5;   // 可调，0.3~0.8 之间效果较好
        // node_voltages = V0 + alpha * (VT - V0); // 临时存储当前初始条件

        // 牛顿法：X0 ← X0 - J^{-1} * F
        Eigen::MatrixXd J(N, N);
        double eps = 1e-9;
        for (int j = 0; j < N; j++) {
            Eigen::VectorXd V0_perturbed = V0;
            V0_perturbed(j) += eps; // 小扰动

            node_voltages = V0_perturbed;

            //扰动后transient求解
            TRAN_solve_with_initial_value(period_T, tstep);
            J.col(j) = (node_voltages - VT) / eps;
        }
        Eigen::MatrixXd JF = J - Eigen::MatrixXd::Identity(N, N);
        Eigen::VectorXd delta_V = JF.fullPivLu().solve(-F);

        node_voltages = V0 + delta_V;
    }
    if (iter == max_iters){
        std::cout << "Warning: Shooting method did NOT converge within max_iters.\n";
    }
    // //展示节点电压结果
    // std::cout << "PSS Analysis Node Voltages:\n";
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

//==================================================
//shooting method大改专用函数
//==================================================

void solver::init_skeleton(double tstep){
    cap_states.clear();
    cap_skeletons.clear();

    ind_skeletons.clear();
    ind_states.clear();

    ckt.extract_MOS_capacitances();

    // Snapshot original linear devices to avoid mutating the container while iterating.
    const std::vector<device> devices = ckt.linear_devices;

    for (const auto& dev : devices){
        if (dev.type == "C"){
            CapacitorState state;
            state.v_prev = 0.0;
            state.i_prev = 0.0;
            cap_states.push_back(state);

            CapacitorSkeleton sk;
            sk.n1 = dev.nodes[0];
            sk.n2 = dev.nodes[1];
            sk.mid = ckt.allocate_internal_node(); // 永久节点
            sk.Req = tstep / (2.0 * dev.parameters.at("value")); // C 改为 Req

            // n1 -- Req -- mid
            ckt.add_resistor(dev.nodes[0], sk.mid, sk.Req);

            // mid -- Veq -- n2
            sk.veq_index = ckt.add_voltage_source(
                sk.mid, dev.nodes[1], 0.0
            );

            cap_skeletons.push_back(sk);
        }
        else if (dev.type == "L"){
            InductorSkeleton sk;
            sk.n1 = dev.nodes[0];
            sk.n2 = dev.nodes[1];
            sk.L  = dev.parameters.at("value");
            sk.Req = 2.0 * sk.L / tstep;

            // 并联等效电阻
            ckt.add_resistor(sk.n1, sk.n2, sk.Req);

            // 等效电流源（初始为 0）
            sk.ieq_index = ckt.add_current_source(
                sk.n1, sk.n2, 0.0
            );

            ind_skeletons.push_back(sk);

            InductorState st;
            st.i_prev = 0.0;
            st.v_prev = 0.0;
            ind_states.push_back(st);
        }
    }
}

void solver::update_capacitor_rhs(){
    for (size_t i = 0; i < cap_skeletons.size(); i++){
        auto& sk = cap_skeletons[i];
        auto& st = cap_states[i];

        double Veq = st.v_prev + sk.Req * st.i_prev;

        CHECK_IDX(sk.veq_index, (int)ckt.sources.size(), "update_capacitor_rhs: sources[veq_index]");
        ckt.sources[sk.veq_index].parameters["DC"] = Veq;
    }
}

void solver::update_capacitor_state()
{
    for (size_t k = 0; k < cap_skeletons.size(); ++k) {
        auto& sk = cap_skeletons[k];
        auto& st = cap_states[k];
        
        if (sk.n1 > 0) CHECK_IDX(sk.n1 - 1, (int)node_voltages.size(), "update_capacitor_state: node_voltages(n1)");
        if (sk.n2 > 0) CHECK_IDX(sk.n2 - 1, (int)node_voltages.size(), "update_capacitor_state: node_voltages(n2)");
        
        double v1 = (sk.n1 == 0) ? 0.0 : node_voltages[sk.n1 - 1];
        double v2 = (sk.n2 == 0) ? 0.0 : node_voltages[sk.n2 - 1];
        
        CHECK_IDX(sk.mid - 1, (int)node_voltages.size(), "update_capacitor_state: node_voltages(mid)");
        double v_mid = node_voltages[sk.mid - 1];

        double iC = (v1 - v_mid) / sk.Req;
        double vC = v1 - v2;

        st.i_prev = iC;
        st.v_prev = vC;
    }
}

void solver::update_inductor_rhs(){
    for (size_t k = 0; k < ind_skeletons.size(); ++k) {
        auto& sk = ind_skeletons[k];
        auto& st = ind_states[k];

        double Ieq = st.i_prev + st.v_prev / sk.Req;

        CHECK_IDX(sk.ieq_index, (int)ckt.sources.size(), "update_inductor_rhs: sources[ieq_index]");
        ckt.sources[sk.ieq_index].parameters["DC"] = Ieq;
    }
}

void solver::update_inductor_state(){
    for (size_t k = 0; k < ind_skeletons.size(); ++k) {
        auto& sk = ind_skeletons[k];
        auto& st = ind_states[k];

        double v1 = (sk.n1 == 0) ? 0.0 : node_voltages[sk.n1 - 1];
        double v2 = (sk.n2 == 0) ? 0.0 : node_voltages[sk.n2 - 1];

        double vL = v1 - v2;

        double iL = st.i_prev + (vL + st.v_prev) / sk.Req;

        st.i_prev = iL;
        st.v_prev = vL;
    }
}

void solver::reset_dynamic_state(){
    for (auto& st : cap_states){
        st.i_prev = 0.0;
        st.v_prev = 0.0;
    }
    for (auto& sk : cap_skeletons){
        ckt.sources[sk.veq_index].parameters["DC"] = 0.0;
    }
    for (auto& st : ind_states){
        st.i_prev = 0.0;
        st.v_prev = 0.0;
    }
    for (auto& sk : ind_skeletons){
        ckt.sources[sk.ieq_index].parameters["DC"] = 0.0;
    }
}


//瞬态步进
void solver::transient_step(double time){
    // (1) 用上一步状态，更新所有动态器件 RHS
    update_capacitor_rhs();
    update_inductor_rhs();

    // (2) 解一次“静态拓扑 + 当前 RHS”的非线性 MNA
    DC_solve_new(time);

    // (3) 从解中提取新的状态，供下一步使用
    update_capacitor_state();
    update_inductor_state();
}

//瞬态求解
void solver::TRAN_solve_new_new(double tstop, double tstep){
    parse_print_variables();
    init_skeleton(tstep);
    node_voltages = Eigen::VectorXd::Zero(static_cast<int>(ckt.node_list.size()) - 1);
    branch_currents.resize(0);

    update_capacitor_rhs();
    update_inductor_rhs();

    DC_solve_new(0.0);

    init_transient();

    TRAN_solve_for_shooting(tstop, tstep);
}

//瞬态主循环
void solver::TRAN_solve_for_shooting(double tstop, double tstep){
    int steps = static_cast<int>(tstop / tstep);

    tran_plot_data.clear();

    for (int step = 0; step <= steps + 1; ++step) {
        double t = step * tstep;
        transient_step(t);
        std::cout << "Time" << "\t" << t << "\n";
        // ===== 记录波形 =====
        for (int plot_node_id : ckt.plot_node_ids) {
            double v = (plot_node_id == 0) ? 0.0 : node_voltages[plot_node_id - 1];
            tran_plot_data[plot_node_id].push_back({t, v});
            std::cout << "node " << ckt.node_list[plot_node_id] << " " << v << "\n";
        }

        for (int idx : ckt.plot_branch_current_indices) {
            CHECK_IDX(idx, (int)branch_currents.size(), "TRAN_solve_for_shooting: branch_currents");
            double i = branch_currents[idx];
            tran_plot_data[-(idx + 1)].push_back({t, i});
        }
    }
}

//初始化瞬态分析
void solver::init_transient()
{
    // ======================================================
    // 1. 电容：根据当前 node_voltages 初始化 v_prev / i_prev
    // ======================================================
    for (size_t k = 0; k < cap_skeletons.size(); ++k) {
        auto& sk = cap_skeletons[k];
        auto& st = cap_states[k];

        double v1 = (sk.n1 == 0) ? 0.0 : node_voltages[sk.mid - 1];
        double v2 = (sk.n2 == 0) ? 0.0 : node_voltages[sk.n2 - 1];

        st.v_prev = v1 - v2;
        st.i_prev = 0.0;   // t=0⁺ 电容电流未知，设 0（标准做法）
    }

    // ======================================================
    // 2. 电感：根据当前 node_voltages 初始化 v_prev / i_prev
    // ======================================================
    for (size_t k = 0; k < ind_skeletons.size(); ++k) {
        auto& sk = ind_skeletons[k];
        auto& st = ind_states[k];

        double v1 = (sk.n1 == 0) ? 0.0 : node_voltages[sk.n1 - 1];
        double v2 = (sk.n2 == 0) ? 0.0 : node_voltages[sk.n2 - 1];

        st.v_prev = v1 - v2;
        st.i_prev = 0.0;   // t=0⁺ 电感电流未知，设 0
    }

    // ======================================================
    // 3. 根据刚初始化的状态，生成 RHS
    // ======================================================
    update_capacitor_rhs();
    update_inductor_rhs();
}


//计算Jacobian
void solver::compute_shooting_jacobian(
    const Eigen::VectorXd& X0,
    Eigen::MatrixXd& J,
    double period_T,
    double tstep
){
    const int N = X0.size();
    const double eps = 1e-6;

    Eigen::VectorXd F0(N), Fi(N);

    // ===== baseline =====
    node_voltages = X0;
    reset_dynamic_state();
    init_transient();
    TRAN_solve_for_shooting(period_T, tstep);
    F0 = node_voltages - X0;

    // ===== perturb each state =====
    for (int i = 0; i < N; ++i) {
        Eigen::VectorXd Xp = X0;
        Xp[i] += eps;

        node_voltages = Xp;
        reset_dynamic_state();
        init_transient();
        TRAN_solve_with_initial_value(period_T, tstep);

        Fi = node_voltages - Xp;

        J.col(i) = (Fi - F0) / eps;
    }
}


//外层调用函数
void solver::PSS_solve_shooting_new_new(double T, double tstep, int max_it, double tol)
{
    // ===== A) 一次性初始化（只做一次）=====
    parse_print_variables();
    // 重要：init_skeleton 会永久改变电路拓扑（加 mid 节点、加等效源）
    // 所以在 shooting 入口里做一次即可，Newton 迭代里绝对不能重复做
    init_skeleton(tstep);

    // node_count 取“当前”node_list（已包含 mid 节点）
    const int node_count = (int)ckt.node_list.size() - 1;
    node_voltages = Eigen::VectorXd::Zero(node_count);
    branch_currents.setZero();

    // 用当前 node_voltages 初始化 cap/ind 的 v_prev/i_prev，并更新 RHS
    init_transient();

    // ===== B) 组装 shooting 状态向量大小 =====
    const int m = (int)cap_states.size();
    const int n = (int)ind_states.size();
    const int N = m + n;

    // 初猜：可以先用 init_transient() 得到的状态，而不是全零
    // Eigen::VectorXd x0 = Eigen::VectorXd::Zero(N);
    // for (int k = 0; k < m; ++k) x0[k] = cap_states[k].v_prev;
    // for (int k = 0; k < n; ++k) x0[m + k] = ind_states[k].i_prev;

    int N_pre_cycles = 3;
    Eigen::VectorXd x0 = compute_x0_by_prerun_TR(T, tstep, N_pre_cycles);
    Eigen::VectorXd V0 = node_voltages;
    Eigen::VectorXd I0 = branch_currents;

    // ===== C) Newton 外层 =====
    const double eps = 1e-6;

    for (int it = 0; it < max_it; ++it) {


        Eigen::VectorXd xT = propagate_one_period(x0, T, tstep);
        Eigen::VectorXd F  = xT - x0;

        V0 = node_voltages;
        I0 = branch_currents;

        double nF = F.norm();
        std::cerr << "[shooting] it=" << it << " |F|=" << nF << "\n";

        if (nF < tol) {
            std::cerr << "[shooting] converged.\n";

            node_voltages = V0;
            branch_currents = I0;

            run_transient_and_record(T, tstep, x0); // 最终稳态波形
            return;
        }

        // 有限差分 Jacobian
        Eigen::MatrixXd J(N, N);
        for (int i = 0; i < N; ++i) {
            Eigen::VectorXd x0p = x0;
            x0p[i] += eps;

            node_voltages = V0;
            branch_currents = I0;

            Eigen::VectorXd xTp = propagate_one_period(x0p, T, tstep);
            Eigen::VectorXd Fp  = xTp - x0p;
            J.col(i) = (Fp - F) / eps;
        }

        Eigen::VectorXd dx = J.partialPivLu().solve(-F);

        // ===== D) 阻尼（先给你最简单稳定版）=====
        double alpha = 1.0;
        // 如果你想更稳：先用 0.5 起步
        // double alpha = 0.5;

        x0 = x0 + alpha * dx;
    }

    std::cerr << "[shooting] WARNING: did not converge in " << max_it << " iterations\n";

    // 即便不收敛，你也可以选择用当前 x0 输出一次波形，方便观察
    run_transient_and_record(T, tstep, x0);
}


void solver::run_transient_and_record(double T, double tstep, const Eigen::VectorXd& x0_star)
{
    // ===== 1. 重置数值状态 =====
    node_voltages.setZero();
    branch_currents.setZero();
    tran_plot_data.clear();

    // ===== 2. 设置 shooting 收敛得到的初始状态 =====
    set_state_from_x0(x0_star);

    // ===== 3. 时间推进并记录 =====
    int steps = (int)std::round(T / tstep);
    for (int s = 0; s <= steps; ++s) {
        double t = s * tstep;

        transient_step(t);
        //std::cout << "Time" << "\t" << t << "\n";
        // ---- 节点电压 ----
        for (int nid : ckt.plot_node_ids) {
            double v = (nid == 0) ? 0.0 : node_voltages[nid - 1];
            tran_plot_data[nid].push_back({t, v});
            //std::cout << "node " << ckt.node_list[nid] << " " << v << "\n";
        }

        // ---- 支路电流 ----
        for (int idx : ckt.plot_branch_current_indices) {
            if (idx < 0 || idx >= (int)branch_currents.size()) {
                std::cerr << "[WARN] branch idx out of range: " << idx << "\n";
                continue;
            }
            double i = branch_currents[idx];
            tran_plot_data[-(idx + 1)].push_back({t, i});
        }
    }
    std::cout <<"\nresults:\n";
    for (int nid : ckt.plot_node_ids) {
        double v = (nid == 0) ? 0.0 : node_voltages[nid - 1];
        std::cout << "node " << ckt.node_list[nid] << " " << tran_plot_data[nid][tran_plot_data[nid].size()-1].second << "\n";
    }
    for (int idx : ckt.plot_branch_current_indices) {
        if (idx < 0 || idx >= (int)branch_currents.size()) {
            std::cerr << "[WARN] branch idx out of range: " << idx << "\n";
            continue;
        }
        print_branch_current(idx);
    }
}

void solver::DC_solve_new(double time) {
    const int max_iter = 1000;
    const double tol = 1e-9;
    const int node_count = static_cast<int>(ckt.node_list.size()) - 1;

    if (node_count <= 0) {
        return;
    }

    if (node_voltages.size() != node_count) {
        node_voltages = Eigen::VectorXd::Zero(node_count);
    }

    for (int iter = 0; iter < max_iter; ++iter) {
        Eigen::VectorXd prev = node_voltages;

        if ((int)node_voltages.size() != node_count) {
            std::cerr << "[WARN] node_voltages.size=" << node_voltages.size()
                    << " node_count=" << node_count << " (resizing)\n";
            node_voltages.conservativeResize(node_count);
            node_voltages.setZero();
        }
        build_MNA_tran(time);
        if ((int)node_voltages.size() != node_count) {
            std::cerr << "[WARN] node_voltages.size=" << node_voltages.size()
                    << " node_count=" << node_count << " (resizing)\n";
            node_voltages.conservativeResize(node_count);
            node_voltages.setZero();
        }

        // Solve the linearized MNA system for this Newton step.
        Eigen::PartialPivLU<Eigen::MatrixXd> lu(MNA_Y);
        Eigen::VectorXd solution = lu.solve(J);

        // Update node voltages and branch currents from the solution vector.
        node_voltages = solution.head(node_count);
        if (solution.size() > node_count) {
            branch_currents = solution.tail(solution.size() - node_count);
        } else {
            branch_currents.resize(0);
        }

        double max_diff = (node_voltages - prev).cwiseAbs().maxCoeff();
        if (max_diff < tol) {
            return;
        }
    }

    //std::cerr << "Warning: DC_solve_new did not converge\n";
}

void solver::build_MNA_tran(double time){
    const int N = ckt.node_list.size() - 1;

    // ===== 0. 初始化 =====
    MNA_Y.resize(N, N);
    MNA_Y.setZero();

    J.resize(N);
    J.setZero();

    // ===== 1. 线性器件 =====
    stamp_linear_devices();

    // ===== 2. 非线性器件（Newton 线性化）=====
    stamp_nonlinear_devices();


    // ===== 3. 独立源 =====
    stamp_independent_sources(time);
}

void solver::stamp_linear_devices(){
    for (auto& dev : ckt.linear_devices){
        char t = dev.type[0];

        if (t == 'R'){
            int n1 = dev.nodes[0];
            int n2 = dev.nodes[1];
            double g = 1.0 / dev.parameters.at("value");

            if (n1) MNA_Y(n1-1, n1-1) += g;
            if (n2) MNA_Y(n2-1, n2-1) += g;
            if (n1 && n2){
                MNA_Y(n1-1, n2-1) -= g;
                MNA_Y(n2-1, n1-1) -= g;
            }
        }

        // ⚠️ C / L 在 transient 中不直接 stamp
        // 已经被 skeleton 替换成 R + source
    }
}

void solver::stamp_nonlinear_devices(){
    build_nonlinear_MNA();
    // const int node_count = static_cast<int>(ckt.node_list.size()) - 1;

    // if (node_voltages.size() < node_count) {
    //     node_voltages.conservativeResize(node_count);
    //     node_voltages.setZero();
    // }

    // auto node_voltage = [&](int node) -> double {
    //     return (node > 0 && node <= node_count) ? node_voltages[node - 1] : 0.0;
    // };

    // for (const auto& dev : ckt.nonlinear_devices){
    //     if (dev.type != "MOS") continue;

    //     const auto& params = dev.parameters;
    //     auto get_param = [&](const char* key, double default_value) {
    //         auto it = params.find(key);
    //         return it == params.end() ? default_value : it->second;
    //     };

    //     int n1 = dev.nodes.size() > 0 ? dev.nodes[0] : 0;
    //     int ng = dev.nodes.size() > 1 ? dev.nodes[1] : 0;
    //     int n2 = dev.nodes.size() > 2 ? dev.nodes[2] : 0;

    //     double W = get_param("W", 1e-6);
    //     double L = get_param("L", 1e-6);
    //     double type = get_param("TYPE", 1.0);

    //     double MU = 0.0, COX = 0.0, VT = 0.0, LAMBDA = 0.0;
    //     if (const model* pmodel = ckt.findModelConst(dev.model)) {
    //         const auto& mp = pmodel->parameters;
    //         auto grab = [&](const char* key, double& dst) {
    //             auto it = mp.find(key);
    //             if (it != mp.end()) dst = it->second;
    //         };
    //         grab("MU", MU);
    //         grab("COX", COX);
    //         grab("VT", VT);
    //         grab("LAMBDA", LAMBDA);
    //     }

    //     double beta = (L != 0.0) ? (MU * COX * (W / L)) : 0.0;

    //     double V1 = node_voltage(n1);
    //     double Vg = node_voltage(ng);
    //     double V2 = node_voltage(n2);

    //     int nd = n1, ns = n2;
    //     double Vd0 = V1, Vs0 = V2;
    //     if ((type > 0 && V1 <= V2) || (type < 0 && V1 > V2)) {
    //         nd = n2;
    //         ns = n1;
    //         Vd0 = V2;
    //         Vs0 = V1;
    //     }

    //     double sign = type;
    //     double Vgs = sign * (Vg - Vs0);
    //     double Vds = sign * (Vd0 - Vs0);
    //     double Vth = sign * VT;

    //     double gm = 0.0;
    //     double gds = 1e-12;
    //     double Ieq = 0.0;

    //     if (Vgs > Vth) {
    //         double Vov = Vgs - Vth;
    //         if (Vds < Vov) {
    //             double Vds_sq = Vds * Vds;
    //             gm = beta * Vds;
    //             gds = beta * (Vov - Vds);
    //             Ieq = beta * (Vov * Vds - 0.5 * Vds_sq) - gm * Vgs - gds * Vds;
    //         } else {
    //             double Vov_sq = Vov * Vov;
    //             gm = beta * Vov * (1.0 + LAMBDA * Vds);
    //             gds = 1e-12 + 0.5 * LAMBDA * beta * Vov;
    //             Ieq = 0.5 * beta * Vov_sq * (1.0 + LAMBDA * Vds) - gm * Vgs - gds * Vds;
    //         }
    //     }

    //     Ieq *= sign;

    //     auto addY = [&](int r, int c, double v) {
    //         if (r && c) MNA_Y(r - 1, c - 1) += v;
    //     };
    //     auto addJ = [&](int node, double v) {
    //         if (node) J(node - 1) += v;
    //     };

    //     addY(nd, ng, gm);
    //     addY(ns, ng, -gm);
    //     addY(nd, ns, -gm);
    //     addY(ns, ns, gm);

    //     addY(nd, nd, gds);
    //     addY(nd, ns, -gds);
    //     addY(ns, nd, -gds);
    //     addY(ns, ns, gds);

    //     addJ(nd, -Ieq);
    //     addJ(ns, Ieq);
    // }
}

void solver::stamp_independent_sources(double time){
    build_sources_MNA(true, time);
}

void solver::set_state_from_x0(const Eigen::VectorXd& x0) {
    const int m = (int)cap_states.size();
    const int n = (int)ind_states.size();
    if (x0.size() != m + n) {
        throw std::runtime_error("x0 size mismatch");
    }

    for (int k = 0; k < m; ++k) {
        cap_states[k].v_prev = x0[k];
        cap_states[k].i_prev = 0.0;     // 先固定为0（可选：也纳入状态）
    }
    for (int k = 0; k < n; ++k) {
        ind_states[k].i_prev = x0[m + k];
        ind_states[k].v_prev = 0.0;
    }

    update_capacitor_rhs();
    update_inductor_rhs();
}

Eigen::VectorXd solver::propagate_one_period(const Eigen::VectorXd& x0, double T, double tstep) {
    // 关键：每次 propagation 只重置“数值状态”，不要再 init_skeleton()
    // node_voltages/branch_currents 这些都回到初始即可
    //node_voltages.setZero();
    //branch_currents.setZero();

    set_state_from_x0(x0);

    int steps = (int)std::round(T / tstep);
    for (int s = 0; s <= steps; ++s) {
        double t = s * tstep;
        transient_step(t);  // 你现成：update_rhs -> DC_solve_new -> update_state
    }

    const int m = (int)cap_states.size();
    const int n = (int)ind_states.size();
    Eigen::VectorXd xT(m + n);
    for (int k = 0; k < m; ++k) xT[k] = cap_states[k].v_prev;
    for (int k = 0; k < n; ++k) xT[m + k] = ind_states[k].i_prev;
    return xT;
}

Eigen::VectorXd solver::compute_x0_by_prerun_TR(double T, double tstep, int N_pre_cycles)
{
    const int m = (int)cap_states.size();
    const int n = (int)ind_states.size();
    const int N = m + n;

    Eigen::VectorXd x0 = Eigen::VectorXd::Zero(N);

    // ---- 安全检查 ----
    if (N == 0) {
        std::cerr << "[pre-run] No dynamic states (no C/L). x0 is empty.\n";
        return x0;
    }
    if (tstep <= 0 || T <= 0) {
        std::cerr << "[pre-run] invalid T/tstep.\n";
        return x0;
    }
    if (N_pre_cycles <= 0) N_pre_cycles = 1;

    // 每周期步数（尽量整数步）
    int steps_per_cycle = (int)std::round(T / tstep);
    if (steps_per_cycle < 1) steps_per_cycle = 1;

    // ======================================================
    // A) 复位数值状态（不要重复 init_skeleton！）
    // ======================================================
    // 复位动态器件内部状态 & RHS 源
    reset_dynamic_state();

    // 复位节点/支路未知量
    const int node_count = (int)ckt.node_list.size() - 1;
    node_voltages = Eigen::VectorXd::Zero(node_count);
    branch_currents.resize(0);

    // ======================================================
    // B) 做一次 t=0 的 DC（用当前 RHS：此时 RHS=0）
    //    然后 init_transient() 用 DC 解初始化状态
    // ======================================================
    update_capacitor_rhs();
    update_inductor_rhs();

    DC_solve_new(0.0);
    init_transient();

    // ======================================================
    // C) 预跑 N 个周期（梯形法 transient_step）
    // ======================================================
    const int total_steps = N_pre_cycles * steps_per_cycle;
    double t = 0.0;
    for (int s = 1; s <= total_steps; ++s) {
        t = s * tstep;
        transient_step(t);

        // 可选：每跑完一个周期打印一下状态范数，方便观察是否趋于周期稳态
        if (s % steps_per_cycle == 0) {
            double cap_norm = 0.0, ind_norm = 0.0;
            for (auto& cs : cap_states) cap_norm += cs.v_prev * cs.v_prev;
            for (auto& is : ind_states) ind_norm += is.i_prev * is.i_prev;
            std::cerr << "[pre-run] cycle=" << (s / steps_per_cycle)
                      << " cap_rms=" << std::sqrt(cap_norm / std::max(1, m))
                      << " ind_rms=" << std::sqrt(ind_norm / std::max(1, n))
                      << "\n";
        }
    }

    // ======================================================
    // D) 打包末端状态为 x0
    // ======================================================
    for (int k = 0; k < m; ++k) x0[k] = cap_states[k].v_prev;
    for (int k = 0; k < n; ++k) x0[m + k] = ind_states[k].i_prev;

    return x0;
}

// ======================================================
// Backward Euler: fast skeleton-based shooting method
// - Capacitor: equivalent voltage source (series with Req)
// - Inductor : equivalent current source (parallel with Req)
// ======================================================

void solver::init_skeleton_BE(double tstep)
{
    if (be_skeleton_initialized) {
        if (be_skeleton_tstep == tstep) return;
        std::cerr << "[BE] WARNING: init_skeleton_BE called again with different tstep; ignoring.\n";
        return;
    }

    be_skeleton_initialized = true;
    be_skeleton_tstep = tstep;

    cap_states_BE.clear();
    cap_skeletons_BE.clear();
    ind_states_BE.clear();
    ind_skeletons_BE.clear();

    ckt.extract_MOS_capacitances();

    const std::vector<device> devices = ckt.linear_devices;
    cap_states_BE.reserve(devices.size());
    cap_skeletons_BE.reserve(devices.size());
    ind_states_BE.reserve(devices.size());
    ind_skeletons_BE.reserve(devices.size());

    for (const auto& dev : devices) {
        if (dev.type == "C") {
            double C = dev.parameters.at("value");

            int st_idx = (int)cap_states_BE.size();
            cap_states_BE.push_back({0.0, 0.0});

            CapacitorSkeleton sk;
            sk.n1 = dev.nodes[0];
            sk.n2 = dev.nodes[1];
            sk.mid = ckt.allocate_internal_node();
            sk.Req = tstep / C;
            sk.state_index = st_idx;

            ckt.add_resistor(sk.n1, sk.mid, sk.Req);
            sk.veq_index = ckt.add_voltage_source(sk.mid, sk.n2, 0.0);

            cap_skeletons_BE.push_back(sk);
        } else if (dev.type == "L") {
            double L = dev.parameters.at("value");

            int st_idx = (int)ind_states_BE.size();
            ind_states_BE.push_back({0.0, 0.0});

            InductorSkeleton sk;
            sk.n1 = dev.nodes[0];
            sk.n2 = dev.nodes[1];
            sk.L  = L;
            sk.Req = L / tstep;
            sk.state_index = st_idx;

            ckt.add_resistor(sk.n1, sk.n2, sk.Req);
            sk.ieq_index = ckt.add_current_source(sk.n1, sk.n2, 0.0);

            ind_skeletons_BE.push_back(sk);
        }
    }
}

void solver::update_capacitor_rhs_BE()
{
    const size_t K = cap_skeletons_BE.size();
    for (size_t k = 0; k < K; ++k) {
        const auto& sk = cap_skeletons_BE[k];
        const auto& st = cap_states_BE[k];
        ckt.sources[sk.veq_index].parameters["DC"] = st.v_prev;
    }
}

void solver::update_capacitor_state_BE()
{
    const int nv = (int)node_voltages.size();
    const size_t K = cap_skeletons_BE.size();

    for (size_t k = 0; k < K; ++k) {
        auto& sk = cap_skeletons_BE[k];
        auto& st = cap_states_BE[k];

        double v1 = (sk.n1 == 0) ? 0.0 : ((sk.n1 - 1 < nv) ? node_voltages[sk.n1 - 1] : 0.0);
        double v2 = (sk.n2 == 0) ? 0.0 : ((sk.n2 - 1 < nv) ? node_voltages[sk.n2 - 1] : 0.0);
        double vm = (sk.mid == 0) ? 0.0 : ((sk.mid - 1 < nv) ? node_voltages[sk.mid - 1] : 0.0);

        double iC = (sk.Req != 0.0) ? (v1 - vm) / sk.Req : 0.0;
        double vC = v1 - v2;

        st.i_prev = iC;
        st.v_prev = vC;
    }
}

void solver::update_inductor_rhs_BE()
{
    const size_t K = ind_skeletons_BE.size();
    for (size_t k = 0; k < K; ++k) {
        const auto& sk = ind_skeletons_BE[k];
        const auto& st = ind_states_BE[k];
        ckt.sources[sk.ieq_index].parameters["DC"] = st.i_prev;
    }
}

void solver::update_inductor_state_BE()
{
    const int nv = (int)node_voltages.size();
    const size_t K = ind_skeletons_BE.size();

    for (size_t k = 0; k < K; ++k) {
        auto& sk = ind_skeletons_BE[k];
        auto& st = ind_states_BE[k];

        double v1 = (sk.n1 == 0) ? 0.0 : ((sk.n1 - 1 < nv) ? node_voltages[sk.n1 - 1] : 0.0);
        double v2 = (sk.n2 == 0) ? 0.0 : ((sk.n2 - 1 < nv) ? node_voltages[sk.n2 - 1] : 0.0);
        double vL = v1 - v2;

        double iL = st.i_prev + ((sk.Req != 0.0) ? (vL / sk.Req) : 0.0);
        st.i_prev = iL;
        st.v_prev = vL;
    }
}

void solver::reset_dynamic_state_BE()
{
    for (auto& st : cap_states_BE) {
        st.v_prev = 0.0;
        st.i_prev = 0.0;
    }
    for (const auto& sk : cap_skeletons_BE) {
        ckt.sources[sk.veq_index].parameters["DC"] = 0.0;
    }

    for (auto& st : ind_states_BE) {
        st.i_prev = 0.0;
        st.v_prev = 0.0;
    }
    for (const auto& sk : ind_skeletons_BE) {
        ckt.sources[sk.ieq_index].parameters["DC"] = 0.0;
    }
}

void solver::init_transient_BE()
{
    const int nv = (int)node_voltages.size();
    auto V = [&](int node) -> double {
        if (node == 0) return 0.0;
        int idx = node - 1;
        return (idx >= 0 && idx < nv) ? node_voltages[idx] : 0.0;
    };

    for (size_t k = 0; k < cap_skeletons_BE.size(); ++k) {
        auto& sk = cap_skeletons_BE[k];
        auto& st = cap_states_BE[k];
        st.v_prev = V(sk.n1) - V(sk.n2);
        st.i_prev = 0.0;
    }

    for (size_t k = 0; k < ind_skeletons_BE.size(); ++k) {
        auto& sk = ind_skeletons_BE[k];
        auto& st = ind_states_BE[k];
        st.v_prev = V(sk.n1) - V(sk.n2);
        st.i_prev = 0.0;
    }

    update_capacitor_rhs_BE();
    update_inductor_rhs_BE();
}

void solver::transient_step_BE(double time)
{
    update_capacitor_rhs_BE();
    update_inductor_rhs_BE();

    DC_solve_new(time);

    update_capacitor_state_BE();
    update_inductor_state_BE();
}

void solver::set_state_from_x0_BE(const Eigen::VectorXd& x0)
{
    const int m = (int)cap_states_BE.size();
    const int n = (int)ind_states_BE.size();
    if (x0.size() != m + n) {
        throw std::runtime_error("x0 size mismatch (BE)");
    }

    for (int k = 0; k < m; ++k) {
        cap_states_BE[k].v_prev = x0[k];
        cap_states_BE[k].i_prev = 0.0;
    }
    for (int k = 0; k < n; ++k) {
        ind_states_BE[k].i_prev = x0[m + k];
        ind_states_BE[k].v_prev = 0.0;
    }

    update_capacitor_rhs_BE();
    update_inductor_rhs_BE();
}

Eigen::VectorXd solver::propagate_one_period_BE(const Eigen::VectorXd& x0, double T, double tstep)
{
    node_voltages.setZero();
    branch_currents.setZero();

    set_state_from_x0_BE(x0);

    int steps = (int)std::round(T / tstep);
    for (int s = 0; s <= steps; ++s) {
        double t = s * tstep;
        transient_step_BE(t);
    }

    const int m = (int)cap_states_BE.size();
    const int n = (int)ind_states_BE.size();
    Eigen::VectorXd xT(m + n);
    for (int k = 0; k < m; ++k) xT[k] = cap_states_BE[k].v_prev;
    for (int k = 0; k < n; ++k) xT[m + k] = ind_states_BE[k].i_prev;
    return xT;
}

Eigen::VectorXd solver::compute_x0_by_prerun_BE(double T, double tstep, int N_pre_cycles)
{
    const int m = (int)cap_states_BE.size();
    const int n = (int)ind_states_BE.size();
    const int N = m + n;

    Eigen::VectorXd x0 = Eigen::VectorXd::Zero(N);

    // ---- 安全检查 ----
    if (N == 0) {
        std::cerr << "[pre-run] No dynamic states (no C/L). x0 is empty.\n";
        return x0;
    }
    if (tstep <= 0 || T <= 0) {
        std::cerr << "[pre-run] invalid T/tstep.\n";
        return x0;
    }
    if (N_pre_cycles <= 0) N_pre_cycles = 1;

    int steps_per_cycle = (int)std::round(T / tstep);
    if (steps_per_cycle < 1) steps_per_cycle = 1;

    reset_dynamic_state_BE();

    const int node_count = (int)ckt.node_list.size() - 1;
    node_voltages = Eigen::VectorXd::Zero(node_count);
    branch_currents.resize(0);

    update_capacitor_rhs_BE();
    update_inductor_rhs_BE();
    DC_solve_new(0.0);
    init_transient_BE();

    const int total_steps = N_pre_cycles * steps_per_cycle;
    for (int s = 1; s <= total_steps; ++s) {
        double t = s * tstep;
        transient_step_BE(t);

        // 可选：每跑完一个周期打印一下状态范数，方便观察是否趋于周期稳态
        if (s % steps_per_cycle == 0) {
            double cap_norm = 0.0, ind_norm = 0.0;
            for (auto& cs : cap_states_BE) cap_norm += cs.v_prev * cs.v_prev;
            for (auto& is : ind_states_BE) ind_norm += is.i_prev * is.i_prev;
            std::cerr << "[pre-run] cycle=" << (s / steps_per_cycle)
                      << " cap_rms=" << std::sqrt(cap_norm / std::max(1, m))
                      << " ind_rms=" << std::sqrt(ind_norm / std::max(1, n))
                      << "\n";
        }
    }

    for (int k = 0; k < m; ++k) x0[k] = cap_states_BE[k].v_prev;
    for (int k = 0; k < n; ++k) x0[m + k] = ind_states_BE[k].i_prev;
    return x0;
}

void solver::run_transient_and_record_BE(double T, double tstep, const Eigen::VectorXd& x0_star)
{
    node_voltages.setZero();
    branch_currents.setZero();
    tran_plot_data.clear();

    set_state_from_x0_BE(x0_star);

    int steps = (int)std::round(T / tstep);
    for (int s = 0; s <= steps; ++s) {
        double t = s * tstep;
        transient_step_BE(t);

        //std::cout << "Time" << "\t" << t << "\n";
        for (int nid : ckt.plot_node_ids) {
            double v = (nid == 0) ? 0.0 : node_voltages[nid - 1];
            tran_plot_data[nid].push_back({t, v});
            //std::cout << "node " << ckt.node_list[nid] << " " << v << "\n";
        }
        for (int idx : ckt.plot_branch_current_indices) {
            if (idx < 0 || idx >= (int)branch_currents.size()) continue;
            double i = branch_currents[idx];
            tran_plot_data[-(idx + 1)].push_back({t, i});
        }
    }
    std::cout <<"\nresults:\n";
    for (int nid : ckt.plot_node_ids) {
        double v = (nid == 0) ? 0.0 : node_voltages[nid - 1];
        std::cout << "node " << ckt.node_list[nid] << " " << tran_plot_data[nid][tran_plot_data[nid].size()-1].second << "\n";
    }
    for (int idx : ckt.plot_branch_current_indices) {
        if (idx < 0 || idx >= (int)branch_currents.size()) {
            std::cerr << "[WARN] branch idx out of range: " << idx << "\n";
            continue;
        }
        print_branch_current(idx);
    }
}

void solver::PSS_solve_shooting_backward_euler(double T, double tstep, int max_it, double tol)
{
    parse_print_variables();

    init_skeleton_BE(tstep);

    const int node_count = (int)ckt.node_list.size() - 1;
    node_voltages = Eigen::VectorXd::Zero(node_count);
    branch_currents.setZero();

    init_transient_BE();

    const int m = (int)cap_states_BE.size();
    const int n = (int)ind_states_BE.size();
    const int N = m + n;

    int N_pre_cycles = 0;
    Eigen::VectorXd x0 = compute_x0_by_prerun_BE(T, tstep, N_pre_cycles);
    Eigen::VectorXd V0 = node_voltages;
    Eigen::VectorXd I0 = branch_currents;

    const double eps = 1e-6;
    for (int it = 0; it < max_it; ++it) {
        V0 = node_voltages;
        I0 = branch_currents;

        Eigen::VectorXd xT = propagate_one_period_BE(x0, T, tstep);
        Eigen::VectorXd F = xT - x0;
        double nF = F.norm();
        std::cerr << "[shooting-BE] it=" << it << " |F|=" << nF << "\n";

        if (nF < tol) {
            std::cerr << "[shooting-BE] converged.\n";
            node_voltages = V0;
            branch_currents = I0;
            run_transient_and_record_BE(T, tstep, x0);
            return;
        }

        Eigen::MatrixXd J(N, N);
        for (int i = 0; i < N; ++i) {
            node_voltages = V0;
            branch_currents = I0;
            Eigen::VectorXd x0p = x0;
            x0p[i] += eps;

            Eigen::VectorXd xTp = propagate_one_period_BE(x0p, T, tstep);
            Eigen::VectorXd Fp = xTp - x0p;
            J.col(i) = (Fp - F) / eps;
        }

        Eigen::VectorXd dx = J.partialPivLu().solve(-F);
        x0 = x0 + dx;
    }

    std::cerr << "[shooting-BE] WARNING: did not converge in " << max_it << " iterations\n";
    run_transient_and_record_BE(T, tstep, x0);
}

void solver::PSS_solve_shooting_backward_euler_sensitivity(double T, double tstep, int max_it, double tol)
{
    parse_print_variables();

    init_skeleton_BE(tstep);

    const int node_count = (int)ckt.node_list.size() - 1;
    node_voltages = Eigen::VectorXd::Zero(std::max(0, node_count));
    branch_currents.resize(0);

    init_transient_BE();

    const int m = (int)cap_states_BE.size();
    const int n = (int)ind_states_BE.size();
    const int N = m + n;

    int N_pre_cycles = 0;
    Eigen::VectorXd x0 = compute_x0_by_prerun_BE(T, tstep, N_pre_cycles);

    if (N == 0) {
        run_transient_and_record_BE(T, tstep, x0);
        return;
    }

    int steps = (int)std::round(T / tstep);
    if (steps < 1) steps = 1;

    using RowMajorMat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

    Eigen::MatrixXd dJ_dx0;
    Eigen::MatrixXd dsol_dx0;
    const Eigen::MatrixXd I = Eigen::MatrixXd::Identity(N, N);

    const int max_iter_dc = 100;
    const double tol_dc = 1e-9;

    auto dc_solve_with_sens = [&](double time, const RowMajorMat& S, Eigen::MatrixXd& dsol_out) -> bool {
        if (node_count <= 0) return true;

        if (node_voltages.size() != node_count) {
            node_voltages = Eigen::VectorXd::Zero(node_count);
        }

        for (int iter = 0; iter < max_iter_dc; ++iter) {
            Eigen::VectorXd prev = node_voltages;

            build_MNA_tran(time);

            Eigen::PartialPivLU<Eigen::MatrixXd> lu(MNA_Y);
            Eigen::VectorXd solution = lu.solve(J);

            node_voltages = solution.head(node_count);
            if (solution.size() > node_count) {
                branch_currents = solution.tail(solution.size() - node_count);
            } else {
                branch_currents.resize(0);
            }

            double max_diff = (node_voltages - prev).cwiseAbs().maxCoeff();
            if (max_diff < tol_dc) {
                const int M = (int)solution.size();
                if (dJ_dx0.rows() != M || dJ_dx0.cols() != N) {
                    dJ_dx0.resize(M, N);
                }
                dJ_dx0.setZero();

                for (int k = 0; k < m; ++k) {
                    const int src_idx = cap_skeletons_BE[k].veq_index;
                    const int bc_idx = ckt.sources[src_idx].branch_current_index;
                    const int row = node_count + bc_idx;
                    if (row >= 0 && row < M) {
                        dJ_dx0.row(row) = S.row(k);
                    }
                }

                for (int k = 0; k < n; ++k) {
                    const int idx = m + k;
                    const auto& sk = ind_skeletons_BE[k];
                    if (sk.n1 != 0) dJ_dx0.row(sk.n1 - 1) -= S.row(idx);
                    if (sk.n2 != 0) dJ_dx0.row(sk.n2 - 1) += S.row(idx);
                }

                dsol_out = lu.solve(dJ_dx0);
                return true;
            }
        }

        std::cerr << "[shooting-BE-sens] WARNING: DC solve did not converge at t=" << time << "\n";
        return false;
    };

    for (int it = 0; it < max_it; ++it) {
        node_voltages.setZero();
        branch_currents.setZero();
        set_state_from_x0_BE(x0);

        RowMajorMat S = RowMajorMat::Identity(N, N);

        for (int s = 0; s <= steps; ++s) {
            double t = s * tstep;

            update_capacitor_rhs_BE();
            update_inductor_rhs_BE();

            if (!dc_solve_with_sens(t, S, dsol_dx0)) {
                break;
            }

            for (int k = 0; k < m; ++k) {
                const auto& sk = cap_skeletons_BE[k];
                S.row(k).setZero();
                if (sk.n1 != 0) S.row(k) += dsol_dx0.row(sk.n1 - 1);
                if (sk.n2 != 0) S.row(k) -= dsol_dx0.row(sk.n2 - 1);
            }

            for (int k = 0; k < n; ++k) {
                const auto& sk = ind_skeletons_BE[k];
                if (sk.Req != 0.0) {
                    const double invReq = 1.0 / sk.Req;
                    if (sk.n1 != 0) S.row(m + k) += dsol_dx0.row(sk.n1 - 1) * invReq;
                    if (sk.n2 != 0) S.row(m + k) -= dsol_dx0.row(sk.n2 - 1) * invReq;
                }
            }

            update_capacitor_state_BE();
            update_inductor_state_BE();
        }

        Eigen::VectorXd xT(N);
        for (int k = 0; k < m; ++k) xT[k] = cap_states_BE[k].v_prev;
        for (int k = 0; k < n; ++k) xT[m + k] = ind_states_BE[k].i_prev;

        Eigen::VectorXd F = xT - x0;
        double nF = F.norm();
        std::cerr << "[shooting-BE-sens] it=" << it << " |F|=" << nF << "\n";

        if (nF < tol) {
            std::cerr << "[shooting-BE-sens] converged.\n";
            run_transient_and_record_BE(T, tstep, x0);
            return;
        }

        Eigen::MatrixXd JF = S - I;
        Eigen::VectorXd dx = JF.partialPivLu().solve(-F);
        x0 += dx;
    }

    std::cerr << "[shooting-BE-sens] WARNING: did not converge in " << max_it << " iterations\n";
    run_transient_and_record_BE(T, tstep, x0);
}
