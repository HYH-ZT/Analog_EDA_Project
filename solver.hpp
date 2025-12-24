#pragma once
#include "circuit.hpp"
#include "parse_netlist.hpp"
#include <Eigen/Dense>
#include <vector>

bool gcd_fundamental_freq(const std::vector<double>& freqs,
                          double& f0_out,
                          double tol = 1e-9);

// 线性方程求解方法枚举
enum class LinearSolverMethod {
    GAUSS_ELIMINATION = 0,  // 高斯消去法
    LU_DECOMPOSITION = 1,   // LU分解法（库函数）
    MANUAL_LU = 2,          // 手动LU分解法
    GAUSS_JACOBI = 3,        // Gauss-Jacobi迭代法
    GAUSS_SEIDEL = 4
};

// 瞬态分析方法枚举
enum class TransientMethod {
    FORWARD_EULER = 0,      // 前向欧拉法
    BACKWARD_EULER = 1,     // 后向欧拉法（隐式）
    TRAPEZOIDAL = 2         // 梯形积分法
};

// 稳态分析方法枚举
enum class SteadyStateMethod {
    SHOOTING = 0,     // 这玩意中文不会真叫射击法吧
    HARMONIC_BALANCE = 1,      // 谐波平衡法
    //CONTINUATION = 2        // 延拓法 这是啥
};

// HB求解器选择
enum class HBLinearSolverMethod {
    EIGEN_LU = 0,
    MANUAL_LU = 1
};

//HB相关参数
class HB_params {
    public:
        double fundamental_omega; // 基频角频率
        int num_harmonics; // 谐波数量
        int max_iterations; // 最大迭代次数
        double tolerance; // 收敛容限
        HBLinearSolverMethod hb_solver_method; // HB线性求解方法
        //初始频域解
        Eigen::VectorXcd initial_xw;
};

class solver {
    private:
        circuit ckt;
        analysis analysis_type;
        
        // 求解方法选择
        LinearSolverMethod linear_solver_method;
        TransientMethod transient_method;
        SteadyStateMethod steady_state_method;
        
        Eigen::VectorXd J; //电流源向量

        //经过主元置换后的J向量
        Eigen::VectorXd J_permuted;
        //线性MNA矩阵
        Eigen::MatrixXd liner_Y;

        //LU分解矩阵
        Eigen::MatrixXd L;
        Eigen::MatrixXd U;
        //最终MNA矩阵
        Eigen::MatrixXd MNA_Y;
        //各节点电压
        Eigen::VectorXd node_voltages;
        //支路电流（用于电压源等）
        Eigen::VectorXd branch_currents;
        //动态器件名称到支路电流索引的映射
        std::map<std::string, int> dynamic_device_current_map;
        //瞬态分析需要绘图的节点电压时间序列
        std::map<int, std::vector<std::pair<double, double>>> tran_plot_data; //节点ID -> (时间, 电压)对列表

        //HB相关变量
        Eigen::VectorXcd hb_J; //HB电流源向量
        Eigen::MatrixXcd hb_liner_Y;
        Eigen::MatrixXcd hb_MNA_Y; //频率点MNA矩阵
        Eigen::MatrixXcd hb_jacobian; //HB雅可比矩阵
        Eigen::MatrixXd t_jacobian;  //时域雅可比矩阵
        //DFT变换矩阵
        Eigen::MatrixXcd hb_DFT_matrix;
        Eigen::MatrixXcd hb_iDFT_matrix;
        Eigen::MatrixXcd hb_F2T_matrix; //把频域解变换为时域解的矩阵
        Eigen::MatrixXcd hb_T2F_matrix; //把时域解变换为频域解的矩阵
        Eigen::VectorXcd hb_xw; //HB频域解
        //Eigen::VectorXcd hb_xt; //HB时域解

        //构建线性MNA矩阵和J向量
        void build_linear_MNA(bool in_tran = false);
        //非线性器件的处理
        void build_nonlinear_MNA();
        //电源的处理
        void build_sources_MNA();
        void build_sources_MNA_ramp(double alpha);
        void build_sources_MNA(bool in_tran,double time);

        //高斯消去法线性MNA方程求解
        void solve_linear_MNA_Gauss();
        //库函数LU分解法求解MNA方程
        void solve_linear_MNA_LU();
        //手动得到LU分解矩阵
        void get_linear_MNA_LU_manual();
        //使用已计算的L和U矩阵求解方程
        void solve_with_LU_matrices();
        //Gauss-Jacobi迭代法求解线性MNA方程
        void solve_linear_MNA_Gauss_Jacobi();
        //Gauss_Seidel迭代法求解线性MNA方程
        void solve_linear_MNA_Gauss_Seidel();
        //根据设置的方法选择线性方程求解算法
        void solve_linear_MNA();
        //根据输入参数，选择矩阵方程求解方法（兼容旧接口）
        //method: 0-高斯消去法，1-LU分解法，2-手动LU分解法，3-Gauss-Jacobi迭代法
        void solve_linear_MNA(int method);
        
        //MNA矩阵操作辅助函数
        void addToY(int rowNode, int colNode, double val);
        void addToJ(int node, double val);
        void store_linear_solution(const Eigen::VectorXd& x);
        //解析打印变量，填充ckt.print_node_ids
        void parse_print_variables();

        //构建瞬态分析电路
        void build_transient_ckt(double tstep);

        //===========================================
        //Harmonic Balance专用变量与函数
        //===========================================
        // int base_size;  //在build_linear_MNA中初始化，表示单频点矩阵大小
        void build_linear_MNA_frequency(Eigen::MatrixXcd& liner_Y_freq, double omega);
        void hb_build_linear_MNA();
        void hb_initialize_DFT_matrices();
        Eigen::VectorXcd hb_DFT(Eigen::VectorXcd xt);
        Eigen::VectorXcd hb_iDFT(Eigen::VectorXcd xw);
        void hb_solve_linear_MNA();
        void hb_build_nonlinear_MNA();
        void hb_build_sources_MNA();
        void hb_compute_jacobian();
        void hb_build_TF_matrix();


        //===========================================
        //shooting method专用变量与函数
        //===========================================
        std::vector<CapacitorState> cap_states;
        std::vector<CapacitorSkeleton> cap_skeletons;

        std::vector<InductorState> ind_states;
        std::vector<InductorSkeleton> ind_skeletons;

        // Backward-Euler shooting dedicated storage (kept separate from TR skeleton/state)
        bool be_skeleton_initialized = false;
        double be_skeleton_tstep = 0.0;
        std::vector<CapacitorState> cap_states_BE;
        std::vector<CapacitorSkeleton> cap_skeletons_BE;
        std::vector<InductorState> ind_states_BE;
        std::vector<InductorSkeleton> ind_skeletons_BE;

        // Common helpers shared by TR/BE shooting workflows
        Eigen::VectorXd propagate_one_period_common(const Eigen::VectorXd& x0, double T, double tstep, bool use_be);
        Eigen::VectorXd compute_x0_by_prerun_common(double T, double tstep, int N_pre_cycles, bool use_be);
        void run_transient_and_record_common(double T, double tstep, const Eigen::VectorXd& x0_star, bool use_be);

    public:
        Eigen::VectorXcd hb_xt;
        int base_size; //基频点矩阵大小
        solver(circuit& ckt_, analysis& analysis_, 
               LinearSolverMethod lsm = LinearSolverMethod::LU_DECOMPOSITION,
               TransientMethod tm = TransientMethod::TRAPEZOIDAL,
               SteadyStateMethod ssm = SteadyStateMethod::SHOOTING);
        // 给节点电压向量设置初值
        void set_initial_node_voltages(std::string node_name, double voltage);
        // 获取ckt中要观察的节点列表
        const std::vector<int>& get_plot_node_ids() const { return ckt.plot_node_ids; }
        const std::vector<int>& get_plot_current_ids() const { return ckt.plot_branch_current_indices; }
        void print_branch_current(int branch_index);
        const int get_voltage_node_size() const { return ckt.node_list.size() - 1; } //不含地节点
        // 输出所有节点电压
        void print_node_voltages();
        // 输出指定节点电压
        void print_node_voltage(const std::string& node_name);
        void print_node_voltage(int node_id);
        // 设置求解方法
        void setLinearSolverMethod(LinearSolverMethod method) { linear_solver_method = method; }
        void setTransientMethod(TransientMethod method) { transient_method = method; }
        void setSteadyStateMethod(SteadyStateMethod method) { steady_state_method = method; }
        void setHBLinearSolverMethod(HBLinearSolverMethod method) { hb_params.hb_solver_method = method; }
        
        // 获取当前求解方法
        LinearSolverMethod getLinearSolverMethod() const { return linear_solver_method; }
        TransientMethod getTransientMethod() const { return transient_method; }
        SteadyStateMethod getSteadyStateMethod() const { return steady_state_method; }
        HBLinearSolverMethod getHBLinearSolverMethod() const { return hb_params.hb_solver_method; }
        
        //直流分析
        void DC_solve();    //默认初值为0
        void DC_solve_ramp();
        void DC_solve(const Eigen::VectorXd& initial_voltages,bool in_tran = false,double time = 0.0);    //根据输入的节点电压向量设置初值，给瞬态用的
        void DC_solve(const std::map<std::string, double>& node_voltage_map,bool in_tran = false);  //根据节点名和电压值的映射设置初值，给人用的
        //瞬态分析
        void TRAN_solve();
        //
        void TRAN_solve(double tstop, double tstep,int use_initial_dc); 
        //
        void TRAN_solve_with_initial_value(double tstop, double tstep);   
        Eigen::MatrixXd TRAN_solve_return(double tstop, double tstep,int use_initial_dc);
        //获取瞬态绘图数据
        const std::map<int, std::vector<std::pair<double, double>>>& get_tran_plot_data() const {
            return tran_plot_data;
        }
        //稳态分析
        //HB相关变量
        HB_params hb_params;
        //设置HB初始频域解
        void HB_set_initial_xw();
        void HB_set_initial_xw_from_transient();
        void PSS_solve();
        void PSS_solve_harmonic_balance();
        void PSS_solve_harmonic_balance(analysis& analysis, int ic_choice, int max_iters = 50, double tol = 1e-6);
        void print_hb_time_domain_results();
        void print_hb_frequency_domain_results();
        void plot_hb_time_domain_results();


    //============================================
    //shooting method专用函数
    //============================================
    void PSS_solve_shooting_trapezoidal(double T, double tstep, int max_it = 100, double tol = 1e-8, int pre_run_cycles = -1);
    void PSS_solve_shooting_trapezoidal_sensitivity(double T, double tstep, int max_it = 100, double tol = 1e-8, int pre_run_cycles = -1);
    void PSS_solve_shooting_backward_euler(double T, double tstep, int max_it = 100, double tol = 1e-8, int pre_run_cycles = -1);
    void PSS_solve_shooting_backward_euler_sensitivity(double T, double tstep, int max_it = 100, double tol = 1e-8, int pre_run_cycles = -1);
    void TRAN_solve_for_shooting_tr(double tstop, double tstep);
    void transient_step_tr(double time);

    void init_skeleton_tr(double tstep);
    void init_transient_tr();

    void update_capacitor_rhs_tr();
    void update_capacitor_state_tr();

    void update_inductor_rhs_tr();
    void update_inductor_state_tr();

    void reset_dynamic_state_tr();
    void set_state_from_x0_tr(const Eigen::VectorXd& x0);

    void TRAN_solve_new_new(double tstop, double tstep);    //测试瞬态
    void build_MNA_tran(double time);

    void stamp_linear_devices();
    void stamp_nonlinear_devices();
    void stamp_independent_sources(double time);

    Eigen::VectorXd propagate_one_period_tr(const Eigen::VectorXd& x0, double T, double tstep);
    void run_transient_and_record_tr(double T, double tstep, const Eigen::VectorXd& x0_star);

    Eigen::VectorXd compute_x0_by_prerun_tr(double T, double tstep, int N_pre_cycles);

    // Backward-Euler shooting (fast skeleton version)
    void init_skeleton_BE(double tstep);
    void init_transient_BE();
    void transient_step_BE(double time);

    void update_capacitor_rhs_BE();
    void update_capacitor_state_BE();

    void update_inductor_rhs_BE();
    void update_inductor_state_BE();

    void reset_dynamic_state_BE();
    void set_state_from_x0_BE(const Eigen::VectorXd& x0);

    Eigen::VectorXd propagate_one_period_BE(const Eigen::VectorXd& x0, double T, double tstep);
    void run_transient_and_record_BE(double T, double tstep, const Eigen::VectorXd& x0_star);
    Eigen::VectorXd compute_x0_by_prerun_BE(double T, double tstep, int N_pre_cycles);
};
