static inline std::string trim(const std::string& s);
bool isComment(const std::string& s);
bool isModel(const std::string &line);
bool isAnalysis(const std::string &line);
static void parseModelLine(const std::string& line, circuit& ckt);
static void parseDeviceLine(const std::string& line, circuit& ckt);
static void parseAnalysisLine(const std::string& line, std::vector<analysis>& analysis_list);
void parseNetlistFile(const std::string& filename,
                      circuit &ckt,
                      std::vector<analysis>& analyses);
