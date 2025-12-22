#include "solver.hpp"
#include "solver_internal.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <utility>
#include <vector>

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

