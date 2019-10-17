#include <sys/resource.h>
#include <chrono>
#include <iostream>
#include <clipp.h>

#include "seedlib.hpp"

struct Timer
{
    Timer()
      : _start(std::chrono::high_resolution_clock::now())
    {}

    std::chrono::milliseconds::duration::rep get_ms() const
    {
        auto stop = std::chrono::high_resolution_clock::now();
        auto dur  = stop - _start;
        return std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
    }

    std::chrono::high_resolution_clock::time_point _start;
};

void
build_index(seedlib::index_config& config)
{
    /// FIXME: rework paritioning such that this is not needed...
    rlimit limits{};
    getrlimit(RLIMIT_NOFILE, &limits);
    limits.rlim_cur = limits.rlim_max;
    setrlimit(RLIMIT_NOFILE, &limits);

    Timer          t;
    seedlib::index index(config);
    std::cerr << "index construction time: " << t.get_ms() << "ms \n";
}

void
stat_index(const std::string& index_name)
{
    seedlib::index index(index_name);
    index.print_memreport();
}

void
query_index(const std::vector<std::string>& fq_in, const std::string& index_name)
{
    using namespace seedlib;

    seedlib::index index(index_name);

    auto on_read
      = [&](query::read_id_t query_id, query::read_pos_t query_length, const std::vector<query::mapping_t>& mappings) {
            //            doNotOptimize(mappings);
            //            return true;
            std::cout << query_id << "(" << query_length << "): \n";
            for (const auto& mapping : mappings) {
                if (mapping.target_id == query_id || mapping.size() < 5) continue;

                std::cout << "\t" << mapping.target_id << "(" << mapping.target_length << "): ";
                for (auto& seed : mapping)
                    std::cout << seed.query_pos << "/" << seed.target_pos << " ";
                std::cout.put('\n');
            }

            return true;
        };

    Timer t;
    query q{index, on_read};
    for (auto& fin : fq_in)
        q.read_fastx(fin);
    std::cerr << "query time: " << t.get_ms() << "ms \n";
}

int
main(int argc, char* argv[])
{
    using namespace clipp;
    // variables storing the parsing result; initialized with their default values
    enum class mode { index, stat, query, help };
    seedlib::index_config config;
    double                downsampling = 1.0;
    mode                  selected     = mode::help;
    bool                  split = false, progr = false;

    auto indexArg   = required("-i", "--index") & value("index", config.output);
    auto inputsArgs = values("fast(a|q) files", config.inputs);

    auto indexCmd = command("index").set(selected, mode::index);
    auto indexMode
      = (indexCmd,
         indexArg,
         option("--b1_len", "-1") & value("b1_len", config.b1) % "b1 length (int [1-8])",
         option("--b2_len", "-2") & value("b2_len", config.b2) % "b2 length (int [0-8], =0 use b1)",
         option("--b3_len", "-3") & value("b3_len", config.b3) % "b3 length (int [0-8], =0 use b1)",
         option("--complexity", "-c") & value("2nuc entropy", config.min_entropy) % "min entropy in bits (float [0-4])",
         option("--downsample", "-d") & value("downsample", downsampling) % "inverse of fraction of positions retained (float [0-1])",
         option("--b1vsb1b2", "-r")
           & value("b1vsb1b2", config.b1_vs_b1b2_drop) % "0 to downsample by b1, 1 to downsample by b1b2 (float [0-1])",
         inputsArgs);

    auto statCmd = command("stat").set(selected, mode::stat);
    auto statMode = (statCmd, indexArg);

    auto findCmd = command("query").set(selected, mode::query);
    auto findMode = (findCmd, indexArg, values("infile", config.inputs));

    auto cli = (indexMode | statMode | findMode | command("help").set(selected, mode::help));

    if (parse(argc, argv, cli)) {
        config.sampling = 1.0 / downsampling;
        switch (selected) {
            case mode::index: build_index(config); break;
            case mode::stat: stat_index(config.output); break;
            case mode::query: query_index(config.inputs, config.output); break;
            case mode::help: std::cout << make_man_page(cli, argv[0]); break;
        }
    } else {
        std::cout << usage_lines(cli, argv[0]) << '\n';
    }
}
