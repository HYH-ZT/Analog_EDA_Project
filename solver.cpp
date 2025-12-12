#include "solver.hpp"
#include <Eigen/Dense>
#include <iostream>
#include "circuit.hpp"
#include <map>
#include <string>
#include <fstream>
#include <chrono>


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
    for (int k = 0; k < n; ++k){
        if (std::abs(MNA_Y(k, k)) < 1e-6){
            //寻找所有行中绝对值最大的元素进行交换
            int max_row = k;
            for (int i = 0; i < n; ++i){
                if (std::abs(MNA_Y(i, k)) > std::abs(MNA_Y(max_row, k))){
                    max_row = i;
                }
            }
            if (max_row != k){
                MNA_Y.row(k).swap(MNA_Y.row(max_row));
                std::swap(J(k), J(max_row));
                // std::cout << "Swapping rows " << k << " and " << max_row << " to avoid zero diagonal\n";
            }
        }
    }

    //使用对角占优化处理，增强收敛性
    Eigen::MatrixXd inversePerm = diagMax();


    //Debug: 输出调整后的MNA矩阵和J向量
    std::cout << "Adjusted MNA_Y Matrix:\n" << MNA_Y << "\n";
    std::cout << "Adjusted J Vector:\n" << J << "\n";

    // std::cout << "Solving MNA equations using Gauss-Jacobi Iteration:\n";
    // //展示特征值
    Eigen::EigenSolver<Eigen::MatrixXd> es(MNA_Y);
    // std::cout << "Eigenvalues of MNA_Y:\n" << es.eigenvalues() << "\n";
    // //如果存在特征值绝对值大于1，则可能不收敛
    for (int i = 0; i < es.eigenvalues().size(); ++i) {
        if (std::abs(es.eigenvalues()(i)) > 1.0) {
            std::cout << "Warning: Eigenvalue " << es.eigenvalues()(i) << " may indicate non-convergence\n";
            break;
        }
    }

    Eigen::VectorXd x_old = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd x_new = Eigen::VectorXd::Zero(n);
    const int max_iterations = 100;
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
        
        //Debug: 输出每次迭代的结果
        std::cout << "Iteration " << iter + 1 << ": x = " << x_new << "\n";

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
    //交换回来
    x_new = inversePerm * x_new;
    
    //存储节点电压结果
    node_voltages.resize(ckt.node_map.size() - 1);
    for (int i = 0; i < ckt.node_map.size() - 1; ++i){
        node_voltages[i] = x_new(i);
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


void solver::build_sources_MNA(bool in_tran,double time){
    //清空打印支路电流索引
    ckt.print_branch_current_indices.clear();
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

        // Debug: 输出当前迭代的 MNA 矩阵和 J 向量
        std::cout << "Iteration " << iter + 1 << ":\n";
        //输出节点序号与名称对照
        for (const auto& pair : ckt.node_map) {
            std::cout << "Node ID " << pair.second << ": " << pair.first << "\n";
        }
        std::cout << "MNA_Y:\n" << MNA_Y << "\n";
        std::cout << "J:\n" << J << "\n";

        // 5. 求解线性方程
        //保存旧节点电压用于收敛性检查
        Eigen::VectorXd old_node_voltages = node_voltages;
        solve_linear_MNA();

        // Debug: 输出当前迭代的节点电压
        std::cout << "Node Voltages:\n" << node_voltages << std::endl << std::endl; 

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
void solver::DC_solve(const Eigen::VectorXd& initial_node_voltages, bool in_tran,double time) {
    const int maxNewtonIter = 500;
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
        // std::cout << "MNA_Y:\n" << MNA_Y << "\n";
        // std::cout << "J:\n" << J << "\n";


        // 5. 求解线性方程
        //保存旧节点电压用于收敛性检查
        Eigen::VectorXd old_node_voltages = node_voltages;
        solve_linear_MNA();

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
            device equiv_current;
            equiv_current.name = "I" + dev.name + "_transient_I";
            equiv_current.type = "I";
            equiv_current.nodes = {n1, mid_node_id};
            equiv_current.node_names = {dev.node_names[0], mid_node_name};
            equiv_current.parameters["DC"] = I_eq;
            equiv_current.original_device_name = dev.name; // 电感等效器件标识
            ckt.sources.push_back(equiv_current);
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
void solver::TRAN_solve(double tstop, double tstep){
    //提取电容信息
    ckt.extract_MOS_capacitances();
    //确定需要打印的变量
    parse_print_variables();

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

    //进行前向欧拉瞬态分析
    // double tstop = analysis_type.parameters["tstop"];
    // double tstep = analysis_type.parameters["tstep"];
    int steps = static_cast<int>(tstop / tstep);
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
        DC_solve(node_voltages, true);

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
            //对Ids_time进行DFT变换，得到频域分量
            // Ids_time 是单一变量在 N 个时间点上的序列（长度 N），
            // 此处直接用 hb_DFT_matrix（N x N）乘以序列，避免将长度为 N 的向量传入
            // 期望长度为 base_size * N 的 hb_DFT 函数导致的越界。
            // Eigen::VectorXcd Ids_freq = hb_DFT_matrix * Ids_time.cast<std::complex<double> >();
            // //将Ids_freq加入到hb_MNA_Y和hb_J中
            // for(int k = 0; k < (2 * hb_params.num_harmonics + 1); ++k){
            //     int row_index;
            //     //对漏极节点加入Ids_freq
            //     row_index = (k) * base_size + nd - 1;
            //     hb_J(row_index) -= Ids_freq(k);
            //     //对源极节点加入-Ids_freq
            //     row_index = (k) * base_size + ns - 1;
            //     hb_J(row_index) += Ids_freq(k);
            // }
           
            
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
    //线性元件对于频域雅可比矩阵的贡献
    hb_jacobian = hb_T2F_matrix * t_jacobian * hb_F2T_matrix;

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

    //确定需要打印的节点
    parse_print_variables();
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
                // for (const auto &d : ckt.sources){
                //     if (d.printI) hdr << "\tI(" << d.name << ")";
                // }
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
                //提取支路电流
                // Eigen::VectorXd branch_currents(ckt.sources.size());
                // for(int i = 0; i < ckt.sources.size(); ++i){
                //     int branch_index = ckt.sources[i].branch_current_index + (ckt.node_list.size() -1);
                //     branch_currents(i) = hb_xt(branch_index + n * base_size).real();
                // }
                out << time;
                for (int node_id : ckt.print_node_ids) {
                    double v = 0.0;
                    if (node_id == 0) v = 0.0;
                    else if (node_id - 1 >= 0 && node_id - 1 < hb_node_voltages.size()) v = hb_node_voltages[node_id - 1];
                    out << "\t" << v;
                }
                out << "\n";
            // for (int current_dev_index : ckt.print_branch_current_indices) {
            //     if(current_dev_index >=0 && current_dev_index < ckt.sources.size()){
            //         out << "\t" << branch_currents[current_dev_index];
            //     }
            // }
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
Eigen::VectorXd solver::run_transient_once(double T, double tstep, const Eigen::VectorXd &init_x)
{
    // 设置初始节点电压
    node_voltages = init_x;

    // 计算步数
    int steps = static_cast<int>(T / tstep);

    for (int step = 0; step <= steps; step++) {
        double time = step * tstep;

        // 构建瞬态分析电路
        build_transient_ckt(tstep);

        // 以当前 node_voltages 为初值求解非线性 MNA
        DC_solve(node_voltages, true);

        // node_voltages 会在 DC_solve 中被更新，无需额外处理
        std::cout << node_voltages.size() << std::endl;
    }

    return node_voltages;   // 即 v(T)
}

void solver::PSS_solve_shooting(double period_T, double tstep, int max_iters, double tol){
    // 节点个数
    int N = ckt.node_map.size() - 1;

    // ---- Step 0：初始化初始条件 X0 ----
    // 你可以用 DC 解，或者直接用 0
    Eigen::VectorXd X0 = Eigen::VectorXd::Zero(N);

    std::cout << "Shooting Method Start: N = " << N << "\n";

    for (int iter = 0; iter < max_iters; iter++)
    {
        std::cout << "=== Shooting Iteration " << iter << " ===\n";

        // ---- Step1：从 X0 出发运行瞬态，得到周期末 v(T) ----
        Eigen::VectorXd XT = run_transient_once(period_T, tstep, X0);

        // ---- Step2：误差 F = XT - X0 ----
        if (XT.size() != X0.size()) {
            X0 = Eigen::VectorXd::Zero(XT.size());
        }
        Eigen::VectorXd F = XT - X0;
        double err = F.norm();

        std::cout << "Error norm = " << err << "\n";

        // ---- Step3：检查收敛 ----
        if (err < tol) {
            std::cout << "Shooting method converged.\n";
            node_voltages = XT;    // 最终稳态
            break;
        }

        // ---- Step4：更新初始条件 ----
        // 松弛法：X0 ← X0 + α*(XT - X0)
        double alpha = 0.5;   // 可调，0.3~0.8 之间效果较好
        X0 = X0 + alpha * F;
    }

    std::cout << "Warning: Shooting method did NOT converge within max_iters.\n";
    node_voltages = X0;

    //展示节点电压结果
    std::cout << "PSS Analysis Node Voltages:\n";
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