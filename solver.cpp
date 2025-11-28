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

void solver::build_linear_MNA(){
    //构建线性MNA矩阵liner_Y和电流源向量J
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
}

void solver::build_nonlinear_MNA() {
    // 假设 MNA_Y 已经被初始化为 liner_Y 的副本（见你 DC_solve 片段）
    // 假设 J 已预先被清零并且大小等于节点数（node_count）
    // node_voltages 存的是上一次迭代的节点电压，索引对应节点编号-1，节点0(ground)不在该向量中

    const double gmin = 1e-12; // cutoff 时的小导纳容错
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
        double KP = 1e-4;    // 默认值μCox（你可按需要改）
        double VTO = 1.0;    // 默认阈值
        double LAMBDA = 0.0; // channel-length modulation
        if (pmodel) {
            if (pmodel->parameters.count("KP")) KP = pmodel->parameters.at("KP");
            if (pmodel->parameters.count("VTO")) VTO = pmodel->parameters.at("VTO");
            if (pmodel->parameters.count("LAMBDA")) LAMBDA = pmodel->parameters.at("LAMBDA");
        }

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
        //根据电压高低确定 Drain 和 Source
        int nd, ns;
        double Vd0, Vs0;

        if (V1 > V2) {
            nd = n1;  Vd0 = V1;
            ns = n2; Vs0 = V2;
        } else {
            nd = n2;  Vd0 = V2;
            ns = n1; Vs0 = V1;
        }

        // 对于 PMOS，把电压和阈值变换为针对模型的常用形式
        // 这里我们通过 TYPE 字段处理：如果 TYPE==-1（p），把电压极性翻转
        double sign = type; // n: +1, p: -1
        double Vgs = sign * (Vg0 - Vs0);
        double Vds = sign * (Vd0 - Vs0);
        double Vth = VTO; // 参考极性已由 sign 处理

        // 计算 Id0, gm, gds
        double Id0 = 0.0;
        double gm = 0.0;
        double gds = 0.0;

        if (Vgs <= Vth) {
            // cutoff
            Id0 = 0.0;
            gm = 0.0;
            gds = gmin;
        } else {
            // 依据 Vds 与 Vgs-Vth 判定工作区
            double Vov = Vgs - Vth; // overdrive
            if (Vds < Vov) {
                // 线性区 (triode)
                // Id = beta * ( (Vov)*Vds - 0.5*Vds^2 ) * (1 + lambda*Vds)
                double Id_no_lambda = beta * (Vov * Vds - 0.5 * Vds * Vds);
                Id0 = Id_no_lambda * (1.0 + LAMBDA * Vds);
                // 导数计算
                // gm = ∂Id/∂Vg = beta * Vds * (1 + lambda*Vds)
                gm = beta * Vds * (1.0 + LAMBDA * Vds);
                // gds = ∂Id/∂Vd = beta * (Vov - Vds) * (1 + lambda*Vds) + Id_no_lambda * lambda
                gds = beta * (Vov - Vds) * (1.0 + LAMBDA * Vds) + Id_no_lambda * LAMBDA;
                if (gds < gmin) gds = gmin;
            } else {
                // 饱和区 (saturation)
                // Id = 0.5 * beta * Vov^2 * (1 + lambda*Vds)
                double Id_no_lambda = 0.5 * beta * Vov * Vov;
                Id0 = Id_no_lambda * (1.0 + LAMBDA * Vds);
                // gm = ∂Id/∂Vg = beta * Vov * (1 + lambda*Vds)
                gm = beta * Vov * (1.0 + LAMBDA * Vds);
                // gds = ∂Id/∂Vd = Id_no_lambda * lambda
                gds = Id_no_lambda * LAMBDA;
                if (gds < gmin) gds = gmin;
            }
        }

        // 对 PMOS，我们计算的 Id0/gm/gds 是基于 sign * voltages 的正向定义。
        // 物理上上述 Id0 表示从 drain -> source（按 sign 方向）。 为简化，我们保持 Id0 表示流出 drain 到 source（即 drain→source）。
        // 若 TYPE==-1（p），sign=-1 已经在 Vgs/Vds 中处理好，Id0 的符号应是正确的。

        // 现在构造线性方程：小信号 i ≈ gds*(vd - vs) + gm*(vg - vs) + Iconst
        // 求常量 Iconst 使在工作点成立：
        // Iconst = Id0 - gds*Vd0 - gm*Vg0 + (gds+gm)*Vs0
        double Iconst = Id0 - gds * Vd0 - gm * Vg0 + (gds + gm) * Vs0;

        // --- 把 gm, gds, Iconst stamp 进 MNA_Y 和 J ---
        // 注意：节点编号为 0 表示地。MNA_Y 行列对应节点编号 1..N -> 索引 0..N-1
        auto addToY = [&](int rowNode, int colNode, double val) {
            if (rowNode == 0 || colNode == 0) return;
            int r = rowNode - 1;
            int c = colNode - 1;
            MNA_Y(r, c) += val;
        };
        auto addToJ = [&](int node, double val) {
            if (node == 0) return;
            int idx = node - 1;
            J(idx) += val;
        };

        // Drain row/col
        // row d: +gds at (d,d); +gm at (d,g); -(gds+gm) at (d,s)
        if (nd != 0) addToY(nd, nd, gds);
        if (ng != 0) addToY(nd, ng, gm);
        if (ns != 0) addToY(nd, ns, -(gds + gm));

        // Source row/col
        // row s: -gds at (s,d); -gm at (s,g); +(gds+gm) at (s,s)
        if (ns != 0) addToY(ns, ns, (gds + gm));
        if (nd != 0) addToY(ns, nd, -gds);
        if (ng != 0) addToY(ns, ng, -gm);

        // Gate 列：我们已经在 (d,g) 和 (s,g) 加了非对角项。Gate 自身对角不增加导纳（理想 MOS gate 无直流电流）
        // 如果 nd/ns/ng 中有 ground(0)，上面对地的 addToY 会被忽略 —— 但要把等效电流接回地

        // 把等效常数电流 Iconst 当作从 drain -> source 的独立电流源来处理：
        // 我们采用约定：J 向量表示注入到节点（positive into node）。
        // 对于一个从 drain 流向 source 的电流 Iconst：
        //   在 drain 节点的注入量应减少 Iconst（因为电流离开 drain），所以 J[d] -= Iconst
        //   在 source 节点的注入量应增加 Iconst（因为电流流入 source），所以 J[s] += Iconst
        if (nd != 0) addToJ(nd, -Iconst);
        if (ns != 0) addToJ(ns, +Iconst);

        // 备注：如果需要考虑 body/bulk (第四个节点)，可以在此加入 body effect （影响 Vth）或把 body 当成源处理。
    } // end for each MOS
}

//直流分析
void solver::DC_solve(){
    //构建直流MNA矩阵
    build_linear_MNA();
    MNA_Y = liner_Y;
    
    for (auto &dev : ckt.nonlinear_devices){
        char c = toupper(dev.name[0]);

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