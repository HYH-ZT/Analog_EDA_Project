#include "solver.hpp"
#include "solver_internal.hpp"
#include <chrono>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>

namespace {
Eigen::VectorXcd solve_lu_partial_pivot(const Eigen::MatrixXcd& A, const Eigen::VectorXcd& b){
    const int n = A.rows();
    Eigen::MatrixXcd lu = A;
    Eigen::VectorXi piv(n);
    for(int i = 0; i < n; ++i){
        piv(i) = i;
    }

    for(int k = 0; k < n; ++k){
        int pivot_row = k;
        double max_abs = 0.0;
        for(int i = k; i < n; ++i){
            double val = std::abs(lu(i, k));
            if(val > max_abs){
                max_abs = val;
                pivot_row = i;
            }
        }
        if(max_abs == 0.0){
            return Eigen::VectorXcd::Zero(n);
        }
        if(pivot_row != k){
            lu.row(k).swap(lu.row(pivot_row));
            std::swap(piv(k), piv(pivot_row));
        }
        const std::complex<double> pivot = lu(k, k);
        for(int i = k + 1; i < n; ++i){
            lu(i, k) /= pivot;
            for(int j = k + 1; j < n; ++j){
                lu(i, j) -= lu(i, k) * lu(k, j);
            }
        }
    }

    Eigen::VectorXcd pb(n);
    for(int i = 0; i < n; ++i){
        pb(i) = b(piv(i));
    }

    Eigen::VectorXcd y(n);
    for(int i = 0; i < n; ++i){
        std::complex<double> sum = pb(i);
        for(int j = 0; j < i; ++j){
            sum -= lu(i, j) * y(j);
        }
        y(i) = sum;
    }

    Eigen::VectorXcd x(n);
    for(int i = n - 1; i >= 0; --i){
        std::complex<double> sum = y(i);
        for(int j = i + 1; j < n; ++j){
            sum -= lu(i, j) * x(j);
        }
        x(i) = sum / lu(i, i);
    }
    return x;
}
}

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

    // //线性元件对于频域雅可比矩阵的贡献
    // hb_jacobian = hb_T2F_matrix * t_jacobian * hb_F2T_matrix;

    //分块来计算矩阵乘法
    //分块法
    // 变换与矩阵乘法（按块实现）用时开始
    auto start_transform = std::chrono::high_resolution_clock::now();

    // 说明：hb_T2F_matrix 和 hb_F2T_matrix 是由每个变量重复的 DFT/IDFT 小块构成。
    // 可以对 t_jacobian 按变量对 (i,j) 提取 NxN 小块 B，计算 hb_DFT_matrix * B * hb_iDFT_matrix，
    // 再把结果写回 hb_jacobian 的相应位置。这样复杂度从 O((base_size*N)^3) 降为 O(base_size^2 * N^3)。

    int N = 2 * hb_params.num_harmonics + 1;
    int total_size = base_size * N;
    hb_jacobian = Eigen::MatrixXcd::Zero(total_size, total_size);

    // 预分配临时矩阵以减少重复分配开销
    Eigen::MatrixXcd block(N, N);
    Eigen::MatrixXcd trans(N, N);

    // 如果希望启用并行化，请在编译时加入 -fopenmp，并取消下面的注释（确保 Eigen 支持线程安全写入）
    // #pragma omp parallel for collapse(2) private(block, trans)
    for(int i = 0; i < base_size; ++i){
        for(int j = 0; j < base_size; ++j){
            // 提取 N x N 小块：block(k,n) = t_jacobian(i + k*base_size, j + n*base_size)
            for(int k = 0; k < N; ++k){
                for(int n = 0; n < N; ++n){
                    block(k, n) = t_jacobian(i + k * base_size, j + n * base_size);
                }
            }

            // 进行小块的频域变换
            trans.noalias() = hb_DFT_matrix * block * hb_iDFT_matrix;

            // 写回结果到 hb_jacobian
            for(int k = 0; k < N; ++k){
                for(int n = 0; n < N; ++n){
                    hb_jacobian(i + k * base_size, j + n * base_size) = trans(k, n);
                }
            }
        }
    }

    // 变换与矩阵乘法用时结束
    auto end_transform = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> transform_duration = end_transform - start_transform;
    std::cout << "Time-Frequency Transform Time: " << transform_duration.count() << " seconds.\n";


    // //Debug: 输出时域雅可比矩阵
    // std::cout << "Time-Domain Jacobian Matrix (t_jacobian):\n" << t_jacobian << "\n";
    // //Debug: 输出频域雅可比矩阵
    // std::cout << "MOS Frequency-Domain Jacobian Matrix (hb_jacobian):\n" << hb_jacobian << "\n";

    hb_jacobian += hb_liner_Y;
}



void solver::hb_solve_linear_MNA(){
    //求解多频率点的线性MNA方程
    Eigen::VectorXcd hb_x;
    if (hb_params.hb_solver_method == HBLinearSolverMethod::MANUAL_LU) {
        hb_x = solve_lu_partial_pivot(hb_MNA_Y, hb_J);
    } else {
        Eigen::PartialPivLU<Eigen::MatrixXcd> lu(hb_MNA_Y);
        hb_x = lu.solve(hb_J);
    }
    // //Debug: 输出多频率点的解向量
    // std::cout << "Harmonic Balance Linear MNA Solution (x):\n" << hb_x << "\n";
    //进行IDFT变换，得到时域解
    // Eigen::VectorXcd hb_xt = hb_iDFT(hb_x);
    hb_xt = hb_F2T_matrix * hb_x;
    // //Debug: 输出时域解向量
    // std::cout << "Harmonic Balance Time-Domain Solution (xt):\n" << hb_xt << "\n";
}

//设置初始频域解

void solver::HB_set_initial_xw(){
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

    // //Debug: 输出初始频域解
    // std::cout << "Initial Harmonic Balance Frequency-Domain Solution (xw):\n" << hb_xw << "\n";
}

//先进行一次瞬态仿真，得到初始时域解，然后进行DFT变换得到频域解，来设定初始频域解
void solver::HB_set_initial_xw_from_transient(){

    //复制一份电路，因为瞬态仿真会修改器件参数
    circuit original_ckt = ckt;

    build_linear_MNA(false);
    MNA_Y = Eigen::MatrixXd::Zero(liner_Y.rows(), liner_Y.cols());
    MNA_Y = liner_Y;
    build_sources_MNA();
    //这样就得到的MNA_Y大小就是每个频率点的矩阵大小
    base_size = MNA_Y.rows();

    hb_initialize_DFT_matrices();

    //根据谐波数量和基频，计算总仿真时间
    double fundamental_period = 1.0 / (hb_params.fundamental_omega / (2.0 * M_PI));
    double time_step = fundamental_period/(2 * hb_params.num_harmonics + 1); //和谐波分析采样点数相同
    double total_time = fundamental_period-time_step; //仿真一个周期
    //进行瞬态仿真
    Eigen::MatrixXd transient_results = TRAN_solve_return(total_time, time_step, false);
    //提取节点电压，进行DFT变换，得到频域解
    int N = 2 * hb_params.num_harmonics + 1;
    //初始化时域解向量
    hb_xt = Eigen::VectorXcd::Zero(base_size * N);
    hb_xw = Eigen::VectorXcd::Zero(base_size * N);
    //提取返回的瞬态结果行数
    int transient_rows = transient_results.rows();
    int transient_cols = transient_results.cols();
    // std::cout << "Transient simulation returned " << transient_rows << " nodes and " << transient_cols << " time points.\n";
    // system("pause");
    //遍历所有节点，提取时域电压
    for(int n = 0; n < transient_rows; ++n){
        // std::cout << n << "\n";
        Eigen::VectorXcd V_time = Eigen::VectorXcd::Zero(N);
        for(int t = 0; t < N; ++t){
            V_time(t) = transient_results(n, t);
        }
        //将时域电压进行DFT变换，得到频域分量
        Eigen::VectorXcd V_freq = hb_DFT_matrix*V_time;

        //将频域分量加入到hb_xw中
        for(int k = 0; k < N; ++k){
            hb_xw(n + k * base_size) = V_freq(k);
        }
    }

    //还原电路
    ckt = original_ckt;
    // HB_set_initial_xw({}); //先设置为空
    build_linear_MNA(false);
    MNA_Y = Eigen::MatrixXd::Zero(liner_Y.rows(), liner_Y.cols());
    MNA_Y = liner_Y;
    build_sources_MNA();
    //这样就得到的MNA_Y大小就是每个频率点的矩阵大小
    base_size = MNA_Y.rows();

}

void solver::PSS_solve_harmonic_balance(){

    //构建多频率点的线性MNA矩阵
    hb_build_linear_MNA();
    // //Debug: 输出多频率点的线性MNA矩阵
    // std::cout << "Harmonic Balance Linear MNA Matrix (Y):\n" << hb_liner_Y << "\n";

    //初始化DFT和IDFT矩阵
    hb_initialize_DFT_matrices();
    hb_build_TF_matrix();

    // //确定需要打印的节点
    // parse_print_variables();
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
    //根据HB求解器设置选择LU实现
    Eigen::VectorXcd delta_xw;
    if (hb_params.hb_solver_method == HBLinearSolverMethod::MANUAL_LU) {
        delta_xw = solve_lu_partial_pivot(hb_jacobian, delta_F);
    } else {
        Eigen::PartialPivLU<Eigen::MatrixXcd> lu(hb_jacobian);
        delta_xw = lu.solve(delta_F);
    }

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
    //给出稳态条件
    node_voltages = Eigen::VectorXd(ckt.node_list.size() -1);
    for(int i = 0; i < (ckt.node_list.size() -1); ++i){
        node_voltages(i) = hb_xt(i).real();
        //std::cout << "Steady-State Voltage at Node " << ckt.node_list[i +1] << ": " << node_voltages(i) << " V\n";
    }
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
        // {
        //     // 输出文件: hb_print.txt
        //         std::ofstream hdr("hb_print.txt", std::ios::out);
        //         hdr << "Time(s)";
        //         for (int node_id : ckt.print_node_ids) {
        //             std::string name = "NODE";
        //             if (node_id >= 0 && node_id < (int)ckt.node_list.size()) name = ckt.node_list[node_id];
        //             hdr << "\tV(" << name << ")";
        //         }
        //         // //只需要遍历所有sources，按顺序输出支路电流表头
        //         for (const auto &d : ckt.sources){
        //             if (d.printI) hdr << "\tI(" << d.name << ")";
        //         }
        //         //关闭
        //         hdr << "\n";
        //         hdr.close();
            

        //     std::ofstream out("hb_print.txt", std::ios::app);
        //     //遍历所有时域点，输出需要打印的节点电压和支路电流
        //     int N = 2 * hb_params.num_harmonics + 1;
        //     double T = 1.0 / (hb_params.fundamental_omega / (2.0 * M_PI)); //周期
        //     double time = 0.0;
        //     for(int n = 0; n < N; ++n){
        //         time = n * T / N;
        //         //提取节点电压
        //         Eigen::VectorXd hb_node_voltages(ckt.node_list.size() -1);
        //         for(int i = 0; i < (ckt.node_list.size() -1); ++i){
        //             hb_node_voltages(i) = hb_xt(i + n * base_size).real();
        //         }
        //         //1215提取支路电流
        //         Eigen::VectorXd branch_currents(ckt.sources.size());
        //         for(int i = 0; i < ckt.sources.size(); ++i){
        //             int branch_index = ckt.sources[i].branch_current_index + (ckt.node_list.size() -1);
        //             branch_currents(i) = hb_xt(branch_index + n * base_size).real();
        //         }
        //         out << time;
        //         for (int node_id : ckt.print_node_ids) {
        //             double v = 0.0;
        //             if (node_id == 0) v = 0.0;
        //             else if (node_id - 1 >= 0 && node_id - 1 < hb_node_voltages.size()) v = hb_node_voltages[node_id - 1];
        //             out << "\t" << v;
        //         }
        //         out << "\n";

        //     //1215 电流
        //     for (int current_dev_index : ckt.print_branch_current_indices) {
        //         if(current_dev_index >=0 && current_dev_index < ckt.sources.size()){
        //             out << "\t" << branch_currents[current_dev_index];
        //         }
        //     }
        //     }

        //     //打印频域结果
        //     out << "\nFrequency Domain Results:\n";
        //     out << "Harmonic\tFrequency(Hz)";
        //     for (int node_id : ckt.print_node_ids) {
        //         std::string name = "NODE";
        //         if (node_id >= 0 && node_id < (int)ckt.node_list.size()) name = ckt.node_list[node_id];
        //         out << "\tV(" << name << ")";
        //     }
        //     out << "\n";
        //     for (int h = 0; h < N; ++h) {
        //         Eigen::VectorXcd hb_node_vw(ckt.node_list.size() -1);
        //         for(int i = 0; i < (ckt.node_list.size() -1); ++i){
        //             hb_node_vw(i) = hb_xw(i + h * base_size);
        //         }
        //         out << h - hb_params.num_harmonics << "\t" << ((h - hb_params.num_harmonics) * (hb_params.fundamental_omega / (2.0 * M_PI)));
        //         for (int node_id : ckt.print_node_ids) {
        //             std::complex<double> v = 0.0;
        //             if (node_id == 0) v = 0.0;
        //             else if (node_id - 1 >= 0 && node_id - 1 < hb_xw.size()) v = hb_node_vw[node_id - 1];
        //             out << "\t" << v;
        //         }
        //         out << "\n";
        //     }

        //     out << "\n";
        //     out.close();
        // }
}


void solver::PSS_solve_harmonic_balance(analysis& analysis, int ic_choice, int max_iters, double tol){
    //开始计时
    auto start_time_pss = std::chrono::high_resolution_clock::now();
    
    //根据网表设定参数
    hb_params.fundamental_omega = 2 * 3.14159265358979323846 * analysis.parameters["freq"]; // 基频角频率
    hb_params.num_harmonics = static_cast<int>(analysis.parameters["harm"]); // 谐波数量
    hb_params.max_iterations = max_iters; // 最大迭代次数
    hb_params.tolerance = tol; // 收敛容限
    ckt.extract_MOS_capacitances(); //提取MOS管寄生电容
        //确定需要打印的节点
    parse_print_variables();
    if(ic_choice == 2){
        //使用瞬态仿真结果作为初始解
        HB_set_initial_xw_from_transient();
    }
    else{
        //使用零初始解
        HB_set_initial_xw(); //空参数表示使用默认初始解
    }
    PSS_solve_harmonic_balance();
    std::cout << "Harmonic Balance Analysis Completed.\n";

    //结束计时
    auto end_time_pss = std::chrono::high_resolution_clock::now();
    auto duration_pss = std::chrono::duration<double>(end_time_pss - start_time_pss).count();
    std::cout << "Total Harmonic Balance Analysis Time: " << duration_pss << " seconds.\n";
}


void solver::print_hb_time_domain_results(){
    //输出到文件中查看,把std::cout改为文件输出
    std::ofstream out("hb_time_domain_results.txt");
    out << "Time(s)";
    for (int node_id : ckt.print_node_ids) {
        std::string name = "NODE";
        if (node_id >= 0 && node_id < (int)ckt.node_list.size()) name = ckt.node_list[node_id];
        out << "\tV(" << name << ")";
    }
    for(const auto &d : ckt.sources){
        if (d.printI) out << "\tI(" << d.name << ")";
    }
    out << "\n";

    int N = 2 * hb_params.num_harmonics + 1;
    double T = 1.0 / (hb_params.fundamental_omega / (2.0 * M_PI)); //周期
    double time = 0.0;
    for(int n = 0; n < N; ++n){
        time = n * T / N;
        //提取节点电压

        Eigen::VectorXd hb_node_voltages(ckt.node_list.size() -1);
        for(int i = 0; i < (ckt.node_list.size() - 1); ++i){
            hb_node_voltages(i) = hb_xt(i + n * base_size).real();
        }
        //提取支路电流
        Eigen::VectorXd branch_currents(ckt.sources.size());
        for(int i = 0; i < ckt.sources.size(); ++i){
            int branch_index = ckt.sources[i].branch_current_index + (ckt.node_list.size() -1);
            branch_currents(i) = hb_xt(branch_index + n * base_size).real();
        }
        out << time;
        for (int node_id : ckt.print_node_ids) {
            double v = 0.0;
            if (node_id == 0) v = 0.0;
            else if (node_id - 1 >= 0 && node_id - 1 < hb_node_voltages.size()) v = hb_node_voltages[node_id - 1];
            out << "\t" << v;
        }
        //打印支路电流
        for (int current_dev_index : ckt.print_branch_current_indices) {
            // std::cout << current_dev_index << "\n";
            if(current_dev_index >=0 && current_dev_index < ckt.sources.size()){
                out << "\t" << branch_currents[current_dev_index];
            }
        }
        out << "\n";
    }
    out.close();
}


void solver::print_hb_frequency_domain_results(){
    //打印频域结果到文件中
    int N = 2 * hb_params.num_harmonics + 1;
    std::ofstream out("hb_frequency_domain_results.txt");
    out << "\nFrequency Domain Results:\n";
    out << "Harmonic\tFrequency(Hz)";
    for (int node_id : ckt.print_node_ids) {
        std::string name = "NODE";
        if (node_id >= 0 && node_id < (int)ckt.node_list.size()) name = ckt.node_list[node_id];
        out << "\tV(" << name << ")";
    }
    for(const auto &d : ckt.sources){
        if (d.printI) out << "\tI(" << d.name << ")";
    }
    out << "\n";
    for (int h = 0; h < N; ++h) {
        Eigen::VectorXcd hb_node_vw(base_size);
        for(int i = 0; i < base_size; ++i){
            hb_node_vw(i) = hb_xw(i + h * base_size);
        }
        //提取支路电流
        Eigen::VectorXcd branch_currents_w(ckt.sources.size());
        for(int i = 0; i < ckt.sources.size(); ++i){
            int branch_index = ckt.sources[i].branch_current_index + (ckt.node_list.size() -1);
            branch_currents_w(i) = hb_xw(branch_index + h * base_size);
        }
        out << h - hb_params.num_harmonics << "\t" << ((h - hb_params.num_harmonics) * (hb_params.fundamental_omega / (2.0 * M_PI)));
        for (int node_id : ckt.print_node_ids) {
            std::complex<double> v = 0.0;
            if (node_id == 0) v = 0.0;
            else if (node_id - 1 >= 0 && node_id - 1 < hb_node_vw.size()) v = hb_node_vw[node_id - 1];
            out << "\t" << v;
        }
        // //打印支路电流
        // for (int current_dev_index : ckt.print_branch_current_indices) {
        //     if(current_dev_index >=0 && current_dev_index < ckt.sources.size()){
        //         std::complex<double> i = hb_node_vw[current_dev_index + (ckt.node_list.size() -1)];
        //         out << "\t" << i;
        //     }
        // }
        //打印支路电流
        for (int current_dev_index : ckt.print_branch_current_indices) {
            // //Debug:
            // std::cout << current_dev_index << "\n";
            if(current_dev_index >=0 && current_dev_index < ckt.sources.size()){
                out << "\t" << branch_currents_w[current_dev_index];
            }
        }
        out << "\n";
    }
    out.close();
}

// void solver::plot_hb_time_domain_results() {
//     // 时间采样点数
//     int N = 2 * hb_params.num_harmonics + 1;
//     double T = 1.0 / (hb_params.fundamental_omega / (2.0 * M_PI)); //周期

//     // 构造时间轴
//     std::vector<double> time_points;
//     for (int n = 0; n < N; ++n) {
//         time_points.push_back(n * T / N);
//     }

//     // 创建一个图
//     plt::figure();

//     for (int node_id : get_plot_node_ids()) {
//         std::string name = "NODE";
//         if (node_id >= 0 && node_id < (int)ckt.node_list.size())
//             name = ckt.node_list[node_id];
//         // 提取节点电压
//         std::vector<double> voltages;
//         for (int n = 0; n < N; ++n) {
//             voltages.push_back(hb_xt((node_id - 1) + n * base_size).real());
//             std::cout << time_points[n] << "\t" << voltages[n] << "\n";
//         }
//         // 绘制节点电压波形
//         plt::plot(time_points, voltages, {{"label", "Node " + name}});
//         std::cout << "Plotted Node " << name << " voltage.\n";
//     }
//     plt::legend();
//     plt::title("HB Time-Domain Node Voltages");
//     plt::xlabel("Time (s)");
//     plt::ylabel("Voltage (V)");
//     plt::grid(true);
//     plt::show();

//     std::cout << "Plotted HB time-domain voltages for all nodes.\n";
// }



