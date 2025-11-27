#include "solver.hpp"
#include <Eigen/Dense>
#include <iostream>
#include "circuit.hpp"

solver::solver(circuit& ckt_, analysis& analysis_):ckt(ckt_), analysis_type(analysis_){
    liner_Y.resize(ckt.node_map.size() - 1, ckt.node_map.size() - 1); //不包含地节点
    liner_Y.setZero();
    J.resize(ckt.node_map.size() - 1);
    J.setZero();
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
    
    //检查解的精度
    double relative_error = (MNA_Y * x - J).norm() / J.norm();
    std::cout << "LU Decomposition Relative Error: " << relative_error << "\n";
    
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
        std::cout << "Branch Currents:\n";
        for (int i = 0; i < branch_currents.size(); ++i) {
            std::cout << "I" << i << ": " << branch_currents[i] << " A\n";
        }
    }
}

//完整LU分解法（手动实现）
void solver::get_linear_MNA_LU_manual(){
    std::cout << "Manual LU Decomposition:\n";
    std::cout << "Original MNA Matrix:\n" << MNA_Y << "\n";
    std::cout << "Original J Vector:\n" << J << "\n";
    
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
            A.row(k).swap(A.row(max_row));
            L.row(k).swap(L.row(max_row));
            //同时交换右端向量J的对应元素
            std::swap(b(k), b(max_row));
            std::swap(pivot_order[k], pivot_order[max_row]);
            std::cout << "Swapping rows " << k << " and " << max_row << " (including J vector)\n";
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
    
    //验证LU分解的正确性
    Eigen::MatrixXd product = L * U;
    Eigen::MatrixXd error_matrix = product - MNA_Y;
    double max_error = error_matrix.cwiseAbs().maxCoeff();
    
    std::cout << "LU Decomposition Results:\n";
    std::cout << "L Matrix:\n" << L << "\n";
    std::cout << "U Matrix:\n" << U << "\n";
    std::cout << "L*U Product:\n" << product << "\n";
    std::cout << "Permuted J Vector:\n" << J_permuted << "\n";
    std::cout << "Maximum decomposition error: " << max_error << "\n";
    
    if (max_error < 1e-10) {
        std::cout << "LU decomposition successful!\n";
    } else {
        std::cout << "Warning: LU decomposition may have numerical issues\n";
    }
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
    std::cout << "Solving using pre-computed L and U matrices:\n";
    std::cout << "Permuted J Vector:\n" << J_permuted << "\n";
    
    //第一步：前向替代求解 L*y = J_permuted
    Eigen::VectorXd y(n);
    for (int i = 0; i < n; ++i) {
        double sum = J_permuted(i);  //使用经过置换的J向量
        for (int j = 0; j < i; ++j) {
            sum -= L(i, j) * y(j);
        }
        y(i) = sum / L(i, i);
    }
    
    std::cout << "Intermediate vector y (L*y = J_permuted):\n" << y << "\n";
    
    //第二步：后向替代求解 U*x = y
    Eigen::VectorXd x(n);
    for (int i = n - 1; i >= 0; --i) {
        double sum = y(i);
        for (int j = i + 1; j < n; ++j) {
            sum -= U(i, j) * x(j);
        }
        x(i) = sum / U(i, i);
    }
    
    std::cout << "Solution vector x (U*x = y):\n" << x << "\n";
    
    //验证解的正确性（使用原始的MNA矩阵和J向量）
    Eigen::VectorXd residual = MNA_Y * x - J;
    double residual_norm = residual.norm();
    std::cout << "Residual norm ||MNA_Y*x - J||: " << residual_norm << "\n";
    
    //存储节点电压结果
    node_voltages.resize(std::min((int)(ckt.node_map.size() - 1), (int)x.size()));
    for (int i = 0; i < node_voltages.size(); ++i){
        node_voltages[i] = x(i);
    }
    
    //如果有额外的变量（如支路电流），也存储起来
    if (x.size() > ckt.node_map.size() - 1) {
        branch_currents.resize(x.size() - (ckt.node_map.size() - 1));
        for (int i = ckt.node_map.size() - 1; i < x.size(); ++i) {
            branch_currents[i - (ckt.node_map.size() - 1)] = x(i);
        }
        std::cout << "Branch Currents:\n";
        for (int i = 0; i < branch_currents.size(); ++i) {
            std::cout << "I" << i << ": " << branch_currents[i] << " A\n";
        }
    }
}

//Gauss-Jacobi迭代法求解线性MNA方程
void solver::solve_linear_MNA_Gauss_Jacobi(){
    std::cout << "Solving MNA equations using Gauss-Jacobi Iteration:\n";
    //展示特征值
    Eigen::EigenSolver<Eigen::MatrixXd> es(MNA_Y);
    std::cout << "Eigenvalues of MNA_Y:\n" << es.eigenvalues() << "\n";
    //如果存在特征值绝对值大于1，则可能不收敛
    for (int i = 0; i < es.eigenvalues().size(); ++i) {
        if (std::abs(es.eigenvalues()(i)) > 1.0) {
            std::cout << "Warning: Eigenvalue " << es.eigenvalues()(i) << " may indicate non-convergence\n";
            break;
        }
    }
    int n = MNA_Y.rows();
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
                std::cout << "Warning: Near-zero diagonal element at row " << i << "\n";
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
    node_voltages.resize(ckt.node_map.size() - 1);
    for (int i = 0; i < ckt.node_map.size() - 1; ++i){
        node_voltages[i] = x_new(i);
    }
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
        // if (c == 'M'){
        //     int d = dev.nodes[0];
        //     int g = dev.nodes[1];
        //     int s = dev.nodes[2];
        //     //简单的DC工作点假设：Vgs=2V, Vds=2V, kn=0.5mA/V^2, Vth=1V
        //     double Vg = (g == 0) ? 0.0 : node_voltages[g - 1];
        //     double Vs = (s == 0) ? 0.0 : node_voltages[s - 1];
        //     double Vd = (d == 0) ? 0.0 : node_voltages[d - 1];
        //     double Vgs = Vg - Vs;
        //     if ()
        //     double kn = dev.parameters["KP"]; //单位A/V^2
        //     double Vth = dev.parameters["VTO"]; //单位V
        //     double Ids = 0.0;
        //     if (Vgs <= Vth){
        //         Ids = 0.0; //截止区
        //     }
        //     else if (Vds < (Vgs - Vth)){
        //         Ids = kn * ((Vgs - Vth) * Vds - 0.5 * Vds * Vds); //线性区
        //     }
        //     else{
        //         Ids = 0.5 * kn * (Vgs - Vth) * (Vgs - Vth); //饱和区
        //     }
        //     if (d != 0){
        //         J(d - 1) -= Ids; //流出漏极为负
        //     }
        //     if (s != 0){
        //         J(s - 1) += Ids; //流入源极为正
        //     }
        // }
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


    ////////////////11.27用于验证线性直流求解结果正确/////////////////////
    //不同方法求解MNA方程
    //solve_linear_MNA_Gauss();
    solve_linear_MNA_LU();
    // get_linear_MNA_LU_manual();
    // solve_with_LU_matrices();
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