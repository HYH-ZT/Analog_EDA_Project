#include "solver.hpp"
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
//存储解变�?
void solver::store_linear_solution(const Eigen::VectorXd& x) {
    const int node_count = (int)ckt.node_map.size() - 1;
    node_voltages.resize(std::max(0, node_count));
    for (int i = 0; i < std::min(node_count, (int)x.size()); ++i) {
        node_voltages[i] = x(i);
    }
    if (x.size() > node_count) {
        branch_currents.resize(x.size() - node_count);
        for (int i = node_count; i < x.size(); ++i) {
            branch_currents[i - node_count] = x(i);
        }
    } else {
        branch_currents.resize(0);
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
        //如果需要，交换�?
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
    store_linear_solution(x);
    //如果有额外的变量（如支路电流），也存储起�?
    // branch_currents are handled inside store_linear_solution
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
    //使用Eigen的LU分解求解�?
    Eigen::PartialPivLU<Eigen::MatrixXd> lu(MNA_Y);
    //检查矩阵是否可�?
    if (lu.determinant() == 0.0) {
        std::cout << "Warning: Matrix is singular or near-singular\n";
    }
    //求解线性方程组 MNA_Y * x = J
    Eigen::VectorXd x = lu.solve(J);
    // //Debug检查解的精�?
    // double relative_error = (MNA_Y * x - J).norm() / J.norm();
    // std::cout << "LU Decomposition Relative Error: " << relative_error << "\n";
    //存储节点电压结果
    store_linear_solution(x);
    //如果有额外的变量（如支路电流），也存储起�?
    // branch_currents are handled inside store_linear_solution
}
//完整LU分解法（手动实现�?
void solver::get_linear_MNA_LU_manual(){ 
    int n = MNA_Y.rows();
    if (n == 0) {
        std::cout << "Error: Empty MNA matrix\n";
        return;
    }
    //初始化L和U矩阵
    L = Eigen::MatrixXd::Identity(n, n);  //L矩阵初始为单位矩�?
    U = Eigen::MatrixXd::Zero(n, n);      //U矩阵初始为零矩阵
    //复制原矩阵到临时矩阵A进行操作
    Eigen::MatrixXd A = MNA_Y;
    //同样需要复制J向量，因为主元交换时需要同步交�?
    Eigen::VectorXd b = J;
    //部分主元选择的置换向�?
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
        //如果需要，交换�?
        if (max_row != k) {
            // 只交换A矩阵和b向量的行，不交换L矩阵
            A.row(k).swap(A.row(max_row));
            std::swap(b(k), b(max_row));
            std::swap(pivot_order[k], pivot_order[max_row]);
            // 交换L矩阵中已计算的部分（第k列之前的元素�?
            for (int j = 0; j < k; ++j) {
                std::swap(L(k, j), L(max_row, j));
            }
            // std::cout << "Swapping rows " << k << " and " << max_row << " (including J vector)\n";
        }
        //检查主元是否接近零
        if (std::abs(A(k, k)) < 1e-12) {
            std::cout << "Warning: Near-zero pivot at position (" << k << "," << k << "): " << A(k, k) << "\n";
        }
        //计算U矩阵的第k�?
        for (int j = k; j < n; ++j) {
            U(k, j) = A(k, j);
        }
        //计算L矩阵的第k列并进行消元
        for (int i = k + 1; i < n; ++i) {
            if (std::abs(U(k, k)) > 1e-12) {
                L(i, k) = A(i, k) / U(k, k);
                //消元：第i行减去L(i,k)*第k�?
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
    //验证解的正确性（使用原始的MNA矩阵和J向量�?
    Eigen::VectorXd residual = MNA_Y * x - J;
    double residual_norm = residual.norm();
    // std::cout << "Residual norm ||MNA_Y*x - J||: " << residual_norm << "\n";
    //存储节点电压结果
    store_linear_solution(x);
}
//Gauss-Jacobi迭代法求解线性MNA方程
void solver::solve_linear_MNA_Gauss_Jacobi(){
    // //Debug：打印MNA矩阵和J向量
    // std::cout << "MNA_Y:\n" << MNA_Y << "\n";
    // std::cout << "J:\n" << J << "\n";
    //先进行行交换，确保对角线元素不为�?
    int n = MNA_Y.rows();
    //对于MNA(i,i)为零的情况，必然是电压源，寻找上方为1或�?1的行，哪个行对应的对角元小，与那行交�?
    for (int i = 0; i < n; ++i){
        if (std::abs(MNA_Y(i, i)) < 1e-12){
            //寻找上方�?或�?1的行
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
            //选择对角元较小的行进行交�?
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
    //构造迭代矩阵，判断收敛�?
    Eigen::MatrixXd D = MNA_Y.diagonal().asDiagonal();
    Eigen::MatrixXd D_inv = D.inverse();
    Eigen::MatrixXd L_plus_U = MNA_Y - D;
    Eigen::MatrixXd iteration_matrix = -D_inv * L_plus_U;
    //计算特征值判断收敛�?
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
        // //Debug: 输出每次迭代的结�?
        // std::cout << "Iteration " << iter + 1 << ": x = " << x_new << "\n";
        //检查收敛�?
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
    //如果有额外的变量（如支路电流），也存储起�?
    if (x_new.size() > ckt.node_map.size() - 1) {
        branch_currents.resize(x_new.size() - (ckt.node_map.size() - 1));
        for (int i = ckt.node_map.size() - 1; i < x_new.size(); ++i) {
            branch_currents[i - (ckt.node_map.size() - 1)] = x_new(i);
        }
    }
}
//Gauss-Seidel迭代法求解线性MNA方程
void solver::solve_linear_MNA_Gauss_Seidel(){
    //先进行行交换，确保对角线元素不为�?
    int n = MNA_Y.rows();
    //对于MNA(i,i)为零的情况，必然是电压源，寻找上方为1或�?1的行，哪个行对应的对角元小，与那行交�?
    for (int i = 0; i < n; ++i){
        if (std::abs(MNA_Y(i, i)) < 1e-12){
            //寻找上方�?或�?1的行
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
            //选择对角元较小的行进行交�?
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
//    Eigen::VectorXd x_old = Eigen::VectorXd::Zero(n);
    //从上次的解开始迭�?
    Eigen::VectorXd x_old = Eigen::VectorXd::Zero(n);
    //把电压和电流结果合并到x_old中，注意可能是空�?
    for (int i = 0; i < std::min((int)node_voltages.size(), n); ++i){
        x_old(i) = node_voltages[i];
    }
    for (int i = 0; i < std::min((int)branch_currents.size(), n - (int)node_voltages.size()); ++i){
        x_old(ckt.node_map.size() - 1 + i) = branch_currents[i];
    }
    Eigen::VectorXd x_new = Eigen::VectorXd::Zero(n);
    const int max_iterations = 5000;
    const double tolerance = 1e-9;
    for (int iter = 0; iter < max_iterations; ++iter) {
        x_new = x_old; //每次迭代开始时，先复制旧�?
        for (int i = 0; i < n; ++i) {
            double sum = J(i);
            for (int j = 0; j < n; ++j) {
                if (j != i) {
                    sum -= MNA_Y(i, j) * x_new(j); //使用最新的x_new�?
                }
            }
            if (std::abs(MNA_Y(i, i)) < 1e-12) {
                // std::cout << "Warning: Near-zero diagonal element at row " << i << "\n";
                x_new(i) = 0.0;
            } else {
                x_new(i) = sum / MNA_Y(i, i);
            }
        }
        //检查收敛�?
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
    //如果有额外的变量（如支路电流），也存储起�?
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
//根据设置的方法选择线性方程求解算�?
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
//method: 0-高斯消去法，1-LU分解法，2-手动LU分解法，3-Gauss-Jacobi迭代�?
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
