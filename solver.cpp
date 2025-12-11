#include "solver.hpp"
#include <Eigen/Dense>
#include <iostream>
#include "circuit.hpp"
#include <map>
#include <string>
#include <fstream>


solver::solver(circuit& ckt_, analysis& analysis_, 
               LinearSolverMethod lsm, TransientMethod tm,
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
        //消元过程
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
    // //Debug: 输出调整后的MNA矩阵和J向量
    // std::cout << "Adjusted MNA_Y Matrix:\n" << MNA_Y << "\n";
    // std::cout << "Adjusted J Vector:\n" << J << "\n";

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
            source.name = "L_short_" + dev.name;
            source.type = "V_DC";
            source.nodes = dev.nodes;
            source.parameters["DC"] = 0.0;
            source.original_device_name = dev.name; // 标记为电感的等效器件
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
            continue;
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
                Id0 = 0.5 * beta * Vov * Vov * (1.0 + LAMBDA * Vds);
                // gm = ∂Id/∂Vg = beta * Vov * (1 + lambda*Vds)
                gm = beta * Vov * (1.0 + LAMBDA * Vds);
                // gds = ∂Id/∂Vd = 0.5 * lambda * beta * Vov^2
                gds = 0.5 * LAMBDA * beta * Vov;
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
        addToY(ns, ns, -gds);
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
            //记录改电压源对应的支路电流变量索引
            dev.branch_current_index = new_var_index- ((int)ckt.node_list.size() - 1);
            //如果需要打印支路电流，存在ckt.print_branch_current_indices,需要减去节点数偏移
            if(dev.printI){
                ckt.print_branch_current_indices.push_back(new_var_index - ((int)ckt.node_list.size() - 1));
            }
            // 如果是瞬态等效电压源，更新动态器件映射
            if (!dev.original_device_name.empty()) {
                dynamic_device_current_map[dev.original_device_name] = new_var_index- ((int)ckt.node_list.size() - 1);
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
    const int maxNewtonIter = 50;
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


    for (int iter = 0; iter < maxNewtonIter; iter++) {

        // 2. 每次迭代重新构造 MNA
        MNA_Y = liner_Y;
        J = Eigen::VectorXd::Zero(MNA_Y.rows());

        // 3. MOS stamp
        build_nonlinear_MNA();

        // 4. 电源 stamp
        build_sources_MNA();

        // Debug: 输出当前迭代的 MNA 矩阵和 J 向量
        std::cout << "Iteration " << iter + 1 << ":\n";
        std::cout << "MNA_Y:\n" << MNA_Y << "\n";
        std::cout << "J:\n" << J << "\n";

        // 5. 求解线性方程
        //保存旧节点电压用于收敛性检查
        Eigen::VectorXd old_node_voltages = node_voltages;
        solve_linear_MNA();

        // // Debug: 输出当前迭代的节点电压
        std::cout << "Node Voltages:\n" << node_voltages << std::endl << std::endl; 

        // 6. 检查收敛性
        double max_diff = (node_voltages - old_node_voltages).cwiseAbs().maxCoeff();
        // std::cout << "Max voltage change: " << max_diff << "\n";
        if (max_diff < tol) {
            std::cout << "Converged after " << iter + 1 << " iterations.\n";
            break;
        }
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

// 使用节点名和电压值的映射来设置初值
void solver::DC_solve(const std::map<std::string, double>& node_voltage_map, bool in_tran) {
    const int maxNewtonIter = 50;
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

    for (int iter = 0; iter < maxNewtonIter; iter++) {

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
void solver::DC_solve(const Eigen::VectorXd& initial_node_voltages, bool in_tran) {
    const int maxNewtonIter = 50;
    const double tol = 1e-9;

    // 1. 只构建一次线性矩阵
    build_linear_MNA(in_tran);

        // // Debug: 输出当前迭代的 MNA 矩阵和 J 向量

        // std::cout << "liner_Y:\n" << liner_Y << "\n";
        // std::cout << "J:\n" << J << "\n";

    int nodeCount = (int)ckt.node_list.size() - 1;
    node_voltages = initial_node_voltages;

    for (int iter = 0; iter < maxNewtonIter; iter++) {

        // 2. 每次迭代重新构造 MNA
        MNA_Y = liner_Y;
        J = Eigen::VectorXd::Zero(MNA_Y.rows());

        // 3. MOS stamp
        build_nonlinear_MNA();

        // // Debug: 输出当前迭代的 MNA 矩阵和 J 向量
        // std::cout << "Iteration " << iter + 1 << ":\n";
        // std::cout << "MNA_Y:\n" << MNA_Y << "\n";
        // std::cout << "J:\n" << J << "\n";

        // 4. 电源 stamp
        build_sources_MNA();

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
            std::cout << "Converged after " << iter + 1 << " iterations.\n";
            break;
        }
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
            double R_eq = tstep / (2.0 * C);
            

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
            double V_eq = (V1_prev - V2_prev) + R_eq * I_prev;
            
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
            double R_eq = 2 * L / tstep;
            
            // 获取上一时刻的电感电流
            // 直接从映射中获取该电感对应的支路电流
            double I_prev = 0.0;
            auto it = dynamic_device_current_map.find(dev.name);
            if (it != dynamic_device_current_map.end() && 
                it->second >= 0 && it->second < branch_currents.size()) {
                I_prev = branch_currents[it->second];
            }
            // 如果找不到，使用初始值0

            double V1_prev = V_prev(n1);
            double V2_prev = V_prev(n2);
            
            double V_eq = -2 * L * I_prev / tstep - (V1_prev - V2_prev);
            
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
            equiv_resistor.original_device_name = dev.name; // 电感等效器件标识
            ckt.linear_devices.push_back(equiv_resistor);
            
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
        std::cout << "Transient Analysis Time: " << time << " s\n";
        //构建瞬态分析电路
        build_transient_ckt(tstep);
        
        //Debug: 输出当前电路的线性器件和源列表
        std::cout << "Linear Devices:\n";
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
        std::cout << "Sources:\n";
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
        std::cout << "\n";

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
    
//瞬态分析，参数可调节
void solver::TRAN_solve(double tstop, double tstep){
    //提取电容信息
    ckt.extract_MOS_capacitances();
    //确定需要打印的变量
    parse_print_variables();

    //先进行直流分析，获得初始条件
    //DC_solve();
    //零初值条件
    //node_voltages = Eigen::VectorXd::Zero(ckt.node_map.size() - 1);

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

        //Debug:展示节点电压结果
        std::cout << "Node Voltages at time " << time << " s:\n";
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
        std::cout << "\n\n";


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


//     ////////////////11.27用于验证线性直流求解结果正确/////////////////////
//     //不同方法求解MNA方程
//     //solve_linear_MNA_Gauss();
//     solve_linear_MNA_LU();
//     // get_linear_MNA_LU_manual();
//     // solve_with_LU_matrices();
//     //展示节点电压结果
//     std::cout << "DC Analysis Node Voltages:\n";
//     for (const auto& pair : ckt.node_map){
//         const std::string& node_name = pair.first;
//         int node_id = pair.second;
//         if (node_id == 0){
//             std::cout << "Node " << node_name << " (ID " << node_id << "): 0 V (Ground)\n";
//         }
//         else{
//             std::cout << "Node " << node_name << " (ID " << node_id << "): " << node_voltages[node_id - 1] << " V\n";
//         }
//     }
// }