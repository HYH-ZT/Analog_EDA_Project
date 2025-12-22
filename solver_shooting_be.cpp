#include "solver.hpp"
#include "solver_internal.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

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
