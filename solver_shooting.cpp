#include "solver.hpp"
#include "solver_internal.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>
static bool almost_integer(double x, double tol = 1e-9) {
    return std::fabs(x - std::llround(x)) <= tol;
}

// 将一组 double 频率近似转成整数并求 gcd，返回基频（Hz）。
// 思路：找一个 scale，使得 fi * scale 都“足够接近整数”；再对整数 gcd。
bool gcd_fundamental_freq(const std::vector<double>& freqs,
                          double& f0_out,
                          double tol)
{
    if (freqs.empty()) return false;

    // 过滤非正数，避免异常
    std::vector<double> f;
    f.reserve(freqs.size());
    for (double x : freqs) {
        if (x > 0.0 && std::isfinite(x)) f.push_back(x);
    }
    if (f.empty()) return false;

    // 尝试 1, 10, 100, ... 1e12 的 scale，使 fi*scale 都接近整数
    // 1e12 对应皮赫兹分辨率；一般电路作业足够用
    const long long scale_candidates[] = {
        1LL, 10LL, 100LL, 1000LL, 10000LL, 100000LL,
        1000000LL, 10000000LL, 100000000LL, 1000000000LL,
        10000000000LL, 100000000000LL, 1000000000000LL
    };

    long long scale = 0;
    for (long long s : scale_candidates) {
        bool ok = true;
        for (double fi : f) {
            double x = fi * static_cast<double>(s);
            if (!almost_integer(x, tol * static_cast<double>(s))) { // 随 scale 放宽一点
                ok = false;
                break;
            }
        }
        if (ok) {
            scale = s;
            break;
        }
    }

    // 如果找不到合适 scale，就退化为“微赫兹量化”（仍然比随便选一个源强）
    if (scale == 0) {
        scale = 1000000LL; // 1e-6 Hz 分辨率
    }

    // 转整数并 gcd
    long long g = 0;
    for (double fi : f) {
        // 注意溢出风险：若频率特别大，可改用 long double 或更小 scale
        long long vi = static_cast<long long>(std::llround(fi * static_cast<double>(scale)));
        if (vi <= 0) continue;
        g = (g == 0) ? vi : std::gcd(g, vi);
        if (g == 1) {
            // gcd 已经为 1（在该 scale 下），可以提前结束
            // 但仍可继续；这里提前结束可省时间
        }
    }

    if (g <= 0) return false;

    f0_out = static_cast<double>(g) / static_cast<double>(scale);
    return (f0_out > 0.0 && std::isfinite(f0_out));
}


//==================================================
//shooting method大改专用函数
//==================================================


void solver::init_skeleton_tr(double tstep){
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


void solver::update_capacitor_rhs_tr(){
    for (size_t i = 0; i < cap_skeletons.size(); i++){
        auto& sk = cap_skeletons[i];
        auto& st = cap_states[i];

        double Veq = st.v_prev + sk.Req * st.i_prev;

        CHECK_IDX(sk.veq_index, (int)ckt.sources.size(), "update_capacitor_rhs_tr: sources[veq_index]");
        ckt.sources[sk.veq_index].parameters["DC"] = Veq;
    }
}


void solver::update_capacitor_state_tr()
{
    for (size_t k = 0; k < cap_skeletons.size(); ++k) {
        auto& sk = cap_skeletons[k];
        auto& st = cap_states[k];
        
        if (sk.n1 > 0) CHECK_IDX(sk.n1 - 1, (int)node_voltages.size(), "update_capacitor_state_tr: node_voltages(n1)");
        if (sk.n2 > 0) CHECK_IDX(sk.n2 - 1, (int)node_voltages.size(), "update_capacitor_state_tr: node_voltages(n2)");
        
        double v1 = (sk.n1 == 0) ? 0.0 : node_voltages[sk.n1 - 1];
        double v2 = (sk.n2 == 0) ? 0.0 : node_voltages[sk.n2 - 1];
        
        CHECK_IDX(sk.mid - 1, (int)node_voltages.size(), "update_capacitor_state_tr: node_voltages(mid)");
        double v_mid = node_voltages[sk.mid - 1];

        double iC = (v1 - v_mid) / sk.Req;
        double vC = v1 - v2;

        st.i_prev = iC;
        st.v_prev = vC;
    }
}


void solver::update_inductor_rhs_tr(){
    for (size_t k = 0; k < ind_skeletons.size(); ++k) {
        auto& sk = ind_skeletons[k];
        auto& st = ind_states[k];

        double Ieq = st.i_prev + st.v_prev / sk.Req;

        CHECK_IDX(sk.ieq_index, (int)ckt.sources.size(), "update_inductor_rhs_tr: sources[ieq_index]");
        ckt.sources[sk.ieq_index].parameters["DC"] = Ieq;
    }
}


void solver::update_inductor_state_tr(){
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


void solver::reset_dynamic_state_tr(){
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

void solver::transient_step_tr(double time){
    // (1) 用上一步状态，更新所有动态器件 RHS
    update_capacitor_rhs_tr();
    update_inductor_rhs_tr();

    // (2) 解一次“静态拓扑 + 当前 RHS”的非线性 MNA
    DC_solve(node_voltages, true, time);

    // (3) 从解中提取新的状态，供下一步使用
    update_capacitor_state_tr();
    update_inductor_state_tr();
}

//瞬态求解

void solver::TRAN_solve_new_new(double tstop, double tstep){
    parse_print_variables();
    init_skeleton_tr(tstep);
    node_voltages = Eigen::VectorXd::Zero(static_cast<int>(ckt.node_list.size()) - 1);
    branch_currents.resize(0);

    update_capacitor_rhs_tr();
    update_inductor_rhs_tr();

    DC_solve(node_voltages, true, 0.0);

    init_transient_tr();

    TRAN_solve_for_shooting_tr(tstop, tstep);
}

//瞬态主循环

void solver::TRAN_solve_for_shooting_tr(double tstop, double tstep){
    int steps = static_cast<int>(tstop / tstep);

    tran_plot_data.clear();

    for (int step = 0; step <= steps + 1; ++step) {
        double t = step * tstep;
        transient_step_tr(t);
        std::cout << "Time" << "\t" << t << "\n";
        // ===== 记录波形 =====
        for (int plot_node_id : ckt.plot_node_ids) {
            double v = (plot_node_id == 0) ? 0.0 : node_voltages[plot_node_id - 1];
            tran_plot_data[plot_node_id].push_back({t, v});
            std::cout << "node " << ckt.node_list[plot_node_id] << " " << v << "\n";
        }

        for (int idx : ckt.plot_branch_current_indices) {
            CHECK_IDX(idx, (int)branch_currents.size(), "TRAN_solve_for_shooting_tr: branch_currents");
            double i = branch_currents[idx];
            tran_plot_data[-(idx + 1)].push_back({t, i});
        }
    }
}

//初始化瞬态分析

void solver::init_transient_tr()
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
    update_capacitor_rhs_tr();
    update_inductor_rhs_tr();
}

//外层调用函数

void solver::PSS_solve_shooting_trapezoidal(double T, double tstep, int max_it, double tol, int pre_run_cycles)
{
    // ===== A) 一次性初始化（只做一次）=====
    parse_print_variables();
    // 重要：init_skeleton_tr 会永久改变电路拓扑（加 mid 节点、加等效源）
    // 所以在 shooting 入口里做一次即可，Newton 迭代里绝对不能重复做
    init_skeleton_tr(tstep);

    // node_count 取“当前”node_list（已包含 mid 节点）
    const int node_count = (int)ckt.node_list.size() - 1;
    node_voltages = Eigen::VectorXd::Zero(node_count);
    branch_currents.setZero();

    // 用当前 node_voltages 初始化 cap/ind 的 v_prev/i_prev，并更新 RHS
    init_transient_tr();

    // ===== B) 组装 shooting 状态向量大小 =====
    const int m = (int)cap_states.size();
    const int n = (int)ind_states.size();
    const int N = m + n;

    // 初猜：可以先用 init_transient_tr() 得到的状态，而不是全零
    // Eigen::VectorXd x0 = Eigen::VectorXd::Zero(N);
    // for (int k = 0; k < m; ++k) x0[k] = cap_states[k].v_prev;
    // for (int k = 0; k < n; ++k) x0[m + k] = ind_states[k].i_prev;

    int N_pre_cycles = pre_run_cycles;
    if (N_pre_cycles < 0) {
        N_pre_cycles = 3;
    }
    Eigen::VectorXd x0 = compute_x0_by_prerun_tr(T, tstep, N_pre_cycles);
    Eigen::VectorXd V0 = node_voltages;
    Eigen::VectorXd I0 = branch_currents;

    // ===== C) Newton 外层 =====
    const double eps = 1e-6;

    for (int it = 0; it < max_it; ++it) {


        Eigen::VectorXd xT = propagate_one_period_tr(x0, T, tstep);
        Eigen::VectorXd F  = xT - x0;

        V0 = node_voltages;
        I0 = branch_currents;

        double nF = F.norm();
        std::cerr << "[shooting] it=" << it << " |F|=" << nF << "\n";

        if (nF < tol) {
            std::cerr << "[shooting] converged.\n";

            node_voltages = V0;
            branch_currents = I0;

            run_transient_and_record_tr(T, tstep, x0); // 最终稳态波形
            return;
        }

        // 有限差分 Jacobian
        Eigen::MatrixXd J(N, N);
        for (int i = 0; i < N; ++i) {
            Eigen::VectorXd x0p = x0;
            x0p[i] += eps;

            node_voltages = V0;
            branch_currents = I0;

            Eigen::VectorXd xTp = propagate_one_period_tr(x0p, T, tstep);
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
    run_transient_and_record_tr(T, tstep, x0);
}


void solver::PSS_solve_shooting_trapezoidal_sensitivity(double T, double tstep, int max_it, double tol, int pre_run_cycles)
{
    parse_print_variables();

    init_skeleton_tr(tstep);

    const int node_count = (int)ckt.node_list.size() - 1;
    node_voltages = Eigen::VectorXd::Zero(std::max(0, node_count));
    branch_currents.resize(0);

    init_transient_tr();

    const int m = (int)cap_states.size();
    const int n = (int)ind_states.size();
    const int N = m + n;

    int N_pre_cycles = pre_run_cycles;
    if (N_pre_cycles < 0) {
        N_pre_cycles = 3;
    }
    Eigen::VectorXd x0 = compute_x0_by_prerun_tr(T, tstep, N_pre_cycles);

    if (N == 0) {
        run_transient_and_record_tr(T, tstep, x0);
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

    auto dc_solve_with_sens = [&](double time,
                                  const RowMajorMat& S,
                                  const RowMajorMat& SiC,
                                  const RowMajorMat& SvL,
                                  Eigen::MatrixXd& dsol_out) -> bool {
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
                    const auto& sk = cap_skeletons[k];
                    const int src_idx = sk.veq_index;
                    const int bc_idx = ckt.sources[src_idx].branch_current_index;
                    const int row = node_count + bc_idx;
                    if (row >= 0 && row < M) {
                        dJ_dx0.row(row) = S.row(k) + sk.Req * SiC.row(k);
                    }
                }

                for (int k = 0; k < n; ++k) {
                    const auto& sk = ind_skeletons[k];
                    const int idx = m + k;
                    Eigen::RowVectorXd dIeq = S.row(idx);
                    if (sk.Req != 0.0) {
                        dIeq += SvL.row(k) * (1.0 / sk.Req);
                    }
                    if (sk.n1 != 0) dJ_dx0.row(sk.n1 - 1) -= dIeq;
                    if (sk.n2 != 0) dJ_dx0.row(sk.n2 - 1) += dIeq;
                }

                dsol_out = lu.solve(dJ_dx0);
                return true;
            }
        }

        std::cerr << "[shooting-TR-sens] WARNING: DC solve did not converge at t=" << time << "\n";
        return false;
    };

    for (int it = 0; it < max_it; ++it) {
        node_voltages.setZero();
        branch_currents.setZero();
        set_state_from_x0_tr(x0);

        RowMajorMat S = RowMajorMat::Identity(N, N);
        RowMajorMat SiC = RowMajorMat::Zero(m, N);
        RowMajorMat SvL = RowMajorMat::Zero(n, N);

        for (int s = 0; s <= steps; ++s) {
            double t = s * tstep;

            update_capacitor_rhs_tr();
            update_inductor_rhs_tr();

            if (!dc_solve_with_sens(t, S, SiC, SvL, dsol_dx0)) {
                break;
            }

            for (int k = 0; k < m; ++k) {
                const auto& sk = cap_skeletons[k];
                S.row(k).setZero();
                if (sk.n1 != 0) S.row(k) += dsol_dx0.row(sk.n1 - 1);
                if (sk.n2 != 0) S.row(k) -= dsol_dx0.row(sk.n2 - 1);

                SiC.row(k).setZero();
                if (sk.Req != 0.0) {
                    const double invReq = 1.0 / sk.Req;
                    if (sk.n1 != 0) SiC.row(k) += dsol_dx0.row(sk.n1 - 1) * invReq;
                    if (sk.mid != 0) SiC.row(k) -= dsol_dx0.row(sk.mid - 1) * invReq;
                }
            }

            for (int k = 0; k < n; ++k) {
                const auto& sk = ind_skeletons[k];

                Eigen::RowVectorXd SvL_prev = SvL.row(k);
                SvL.row(k).setZero();
                if (sk.n1 != 0) SvL.row(k) += dsol_dx0.row(sk.n1 - 1);
                if (sk.n2 != 0) SvL.row(k) -= dsol_dx0.row(sk.n2 - 1);

                if (sk.Req != 0.0) {
                    const double invReq = 1.0 / sk.Req;
                    S.row(m + k) += (SvL_prev + SvL.row(k)) * invReq;
                }
            }

            update_capacitor_state_tr();
            update_inductor_state_tr();
        }

        Eigen::VectorXd xT(N);
        for (int k = 0; k < m; ++k) xT[k] = cap_states[k].v_prev;
        for (int k = 0; k < n; ++k) xT[m + k] = ind_states[k].i_prev;

        Eigen::VectorXd F = xT - x0;
        double nF = F.norm();
        std::cerr << "[shooting-TR-sens] it=" << it << " |F|=" << nF << "\n";

        if (nF < tol) {
            std::cerr << "[shooting-TR-sens] converged.\n";
            run_transient_and_record_tr(T, tstep, x0);
            return;
        }

        Eigen::MatrixXd JF = S - I;
        Eigen::VectorXd dx = JF.partialPivLu().solve(-F);
        x0 += dx;
    }

    std::cerr << "[shooting-TR-sens] WARNING: did not converge in " << max_it << " iterations\n";
    run_transient_and_record_tr(T, tstep, x0);
}


void solver::run_transient_and_record_tr(double T, double tstep, const Eigen::VectorXd& x0_star)
{
    run_transient_and_record_common(T, tstep, x0_star, false);
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
}


void solver::stamp_independent_sources(double time){
    build_sources_MNA(true, time);
}


void solver::set_state_from_x0_tr(const Eigen::VectorXd& x0) {
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

    update_capacitor_rhs_tr();
    update_inductor_rhs_tr();
}


Eigen::VectorXd solver::propagate_one_period_tr(const Eigen::VectorXd& x0, double T, double tstep) {
    return propagate_one_period_common(x0, T, tstep, false);
}


Eigen::VectorXd solver::compute_x0_by_prerun_tr(double T, double tstep, int N_pre_cycles)
{
    return compute_x0_by_prerun_common(T, tstep, N_pre_cycles, false);
}

Eigen::VectorXd solver::propagate_one_period_common(const Eigen::VectorXd& x0, double T, double tstep, bool use_be)
{
    node_voltages.setZero();
    branch_currents.setZero();

    if (use_be) {
        set_state_from_x0_BE(x0);
    } else {
        set_state_from_x0_tr(x0);
    }

    int steps = (int)std::round(T / tstep);
    for (int s = 0; s <= steps; ++s) {
        double t = s * tstep;
        if (use_be) {
            transient_step_BE(t);
        } else {
            transient_step_tr(t);
        }
    }

    auto& caps = use_be ? cap_states_BE : cap_states;
    auto& inds = use_be ? ind_states_BE : ind_states;
    const int m = (int)caps.size();
    const int n = (int)inds.size();
    Eigen::VectorXd xT(m + n);
    for (int k = 0; k < m; ++k) xT[k] = caps[k].v_prev;
    for (int k = 0; k < n; ++k) xT[m + k] = inds[k].i_prev;
    return xT;
}

Eigen::VectorXd solver::compute_x0_by_prerun_common(double T, double tstep, int N_pre_cycles, bool use_be)
{
    auto& caps = use_be ? cap_states_BE : cap_states;
    auto& inds = use_be ? ind_states_BE : ind_states;
    const int m = (int)caps.size();
    const int n = (int)inds.size();
    const int N = m + n;

    Eigen::VectorXd x0 = Eigen::VectorXd::Zero(N);

    if (N == 0) {
        std::cerr << "[pre-run] No dynamic states (no C/L). x0 is empty.\n";
        return x0;
    }
    if (tstep <= 0 || T <= 0) {
        std::cerr << "[pre-run] invalid T/tstep.\n";
        return x0;
    }
    if (N_pre_cycles < 0) N_pre_cycles = 1;
    if (N_pre_cycles == 0) return x0;

    int steps_per_cycle = (int)std::round(T / tstep);
    if (steps_per_cycle < 1) steps_per_cycle = 1;

    if (use_be) {
        reset_dynamic_state_BE();
    } else {
        reset_dynamic_state_tr();
    }

    const int node_count = (int)ckt.node_list.size() - 1;
    node_voltages = Eigen::VectorXd::Zero(node_count);
    branch_currents.resize(0);

    if (use_be) {
        update_capacitor_rhs_BE();
        update_inductor_rhs_BE();
    } else {
        update_capacitor_rhs_tr();
        update_inductor_rhs_tr();
    }

    DC_solve(node_voltages, true, 0.0);

    if (use_be) {
        init_transient_BE();
    } else {
        init_transient_tr();
    }

    const int total_steps = N_pre_cycles * steps_per_cycle;
    for (int s = 1; s <= total_steps; ++s) {
        double t = s * tstep;
        if (use_be) {
            transient_step_BE(t);
        } else {
            transient_step_tr(t);
        }

        if (s % steps_per_cycle == 0) {
            double cap_norm = 0.0, ind_norm = 0.0;
            for (auto& cs : caps) cap_norm += cs.v_prev * cs.v_prev;
            for (auto& is : inds) ind_norm += is.i_prev * is.i_prev;
            std::cerr << "[pre-run] cycle=" << (s / steps_per_cycle)
                      << " cap_rms=" << std::sqrt(cap_norm / std::max(1, m))
                      << " ind_rms=" << std::sqrt(ind_norm / std::max(1, n))
                      << "\n";
        }
    }

    for (int k = 0; k < m; ++k) x0[k] = caps[k].v_prev;
    for (int k = 0; k < n; ++k) x0[m + k] = inds[k].i_prev;
    return x0;
}

void solver::run_transient_and_record_common(double T, double tstep, const Eigen::VectorXd& x0_star, bool use_be)
{
    node_voltages.setZero();
    branch_currents.setZero();
    tran_plot_data.clear();

    if (use_be) {
        set_state_from_x0_BE(x0_star);
    } else {
        set_state_from_x0_tr(x0_star);
    }

    int steps = (int)std::round(T / tstep);
    for (int s = 0; s <= steps; ++s) {
        double t = s * tstep;
        if (use_be) {
            transient_step_BE(t);
        } else {
            transient_step_tr(t);
        }

        for (int nid : ckt.plot_node_ids) {
            double v = (nid == 0) ? 0.0 : node_voltages[nid - 1];
            tran_plot_data[nid].push_back({t, v});
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

    DC_solve(node_voltages, true, time);

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
    return propagate_one_period_common(x0, T, tstep, true);
}


Eigen::VectorXd solver::compute_x0_by_prerun_BE(double T, double tstep, int N_pre_cycles)
{
    return compute_x0_by_prerun_common(T, tstep, N_pre_cycles, true);
}


void solver::run_transient_and_record_BE(double T, double tstep, const Eigen::VectorXd& x0_star)
{
    run_transient_and_record_common(T, tstep, x0_star, true);
}


void solver::PSS_solve_shooting_backward_euler(double T, double tstep, int max_it, double tol, int pre_run_cycles)
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

    int N_pre_cycles = pre_run_cycles;
    if (N_pre_cycles < 0) {
        N_pre_cycles = 0;
    }
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


void solver::PSS_solve_shooting_backward_euler_sensitivity(double T, double tstep, int max_it, double tol, int pre_run_cycles)
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

    int N_pre_cycles = pre_run_cycles;
    if (N_pre_cycles < 0) {
        N_pre_cycles = 0;
    }
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
