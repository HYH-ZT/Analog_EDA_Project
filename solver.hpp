#pragma once
#include "circuit.hpp"
#include "parse_netlist.hpp"
#include <Eigen/Dense>

// 线性方程求解方法枚举
enum class LinearSolverMethod {
    GAUSS_ELIMINATION = 0,  // 高斯消去法
    LU_DECOMPOSITION = 1,   // LU分解法（库函数）
    MANUAL_LU = 2,          // 手动LU分解法
    GAUSS_JACOBI = 3        // Gauss-Jacobi迭代法
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

//HB相关参数
class HB_params {
    public:
        double fundamental_omega; // 基频角频率
        int num_harmonics; // 谐波数量
        int max_iterations; // 最大迭代次数
        double tolerance; // 收敛容限
        double relaxation_factor; // 松弛因子
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
        Eigen::VectorXcd hb_xt; //HB时域解

        //构建线性MNA矩阵和J向量
        void build_linear_MNA(bool in_tran = false);
        //非线性器件的处理
        void build_nonlinear_MNA();
        //电源的处理
        void build_sources_MNA();
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
        //根据设置的方法选择线性方程求解算法
        void solve_linear_MNA();
        //根据输入参数，选择矩阵方程求解方法（兼容旧接口）
        //method: 0-高斯消去法，1-LU分解法，2-手动LU分解法，3-Gauss-Jacobi迭代法
        void solve_linear_MNA(int method);
        
        //MNA矩阵操作辅助函数
        void addToY(int rowNode, int colNode, double val);
        void addToJ(int node, double val);
        //解析打印变量，填充ckt.print_node_ids
        void parse_print_variables();

        //构建瞬态分析电路
        void build_transient_ckt(double tstep);

        //HB相关函数
        int base_size;  //在build_linear_MNA中初始化，表示单频点矩阵大小
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


    public:
        solver(circuit& ckt_, analysis& analysis_, 
               LinearSolverMethod lsm = LinearSolverMethod::LU_DECOMPOSITION,
               TransientMethod tm = TransientMethod::TRAPEZOIDAL,
               SteadyStateMethod ssm = SteadyStateMethod::SHOOTING);
        // 给节点电压向量设置初值
        void set_initial_node_voltages(std::string node_name, double voltage);
        // 获取ckt中要观察的节点列表
        const std::vector<int>& get_plot_node_ids() const { return ckt.plot_node_ids; }
        // 输出所有节点电压
        void print_node_voltages();
        // 输出指定节点电压
        void print_node_voltage(const std::string& node_name);
        void print_node_voltage(int node_id);
        // 设置求解方法
        void setLinearSolverMethod(LinearSolverMethod method) { linear_solver_method = method; }
        void setTransientMethod(TransientMethod method) { transient_method = method; }
        void setSteadyStateMethod(SteadyStateMethod method) { steady_state_method = method; }
        
        // 获取当前求解方法
        LinearSolverMethod getLinearSolverMethod() const { return linear_solver_method; }
        TransientMethod getTransientMethod() const { return transient_method; }
        SteadyStateMethod getSteadyStateMethod() const { return steady_state_method; }
        
        //直流分析
        void DC_solve();    //默认初值为0
        void DC_solve(const Eigen::VectorXd& initial_voltages,bool in_tran = false,double time = 0.0);    //根据输入的节点电压向量设置初值，给瞬态用的
        void DC_solve(const std::map<std::string, double>& node_voltage_map,bool in_tran = false);  //根据节点名和电压值的映射设置初值，给人用的
        //瞬态分析
        void TRAN_solve();
        //
        void TRAN_solve(double tstop, double tstep);    //测试用
        //稳态分析
        //HB相关变量
        HB_params hb_params;
        //设置HB初始频域解
        void HB_set_initial_xw(const std::map<std::string, double>& node_voltage_map);
        void PSS_solve();
        void PSS_solve_shooting();
        void PSS_solve_harmonic_balance();



        //神秘实验内容：主元最优置换
    bool bfs_match(const std::vector<std::vector<int>>& graph, 
                   std::vector<int>& pairU, 
                   std::vector<int>& pairV, 
                   std::vector<int>& dist, 
                   int N);

    bool dfs_match(int u, const std::vector<std::vector<int>>& graph,
                   std::vector<int>& pairU, std::vector<int>& pairV,
                   std::vector<int>& dist);

    std::vector<int> hungarianMatch(const Eigen::MatrixXd& matrix);
    Eigen::MatrixXd diagMax();
    static Eigen::VectorXd restoreOrder(const Eigen::VectorXd& solution, 
                                        const Eigen::MatrixXd& inversePerm);
    int checkDiagonalDominance(double threshold = 1e-6) const ;

    Eigen::VectorXd run_transient_once(double T, double tstep, const Eigen::VectorXd &init_x);
    void PSS_solve_shooting(double period_T, double tstep, int max_iters, double tol);

};

