#include <sys/resource.h>
#include <chrono>
#include <iostream>
#include <unordered_map>
#include <unordered_set>

#include <clipp.h>

#include "seedlib/seedlib.hpp"

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
query_index_clique(const std::vector<std::string>& fq_in, const seedlib::index& index)
{
    using namespace seedlib;

    auto on_read
      = [&](query::read_id_t query_id, query::read_pos_t query_length, const std::vector<query::mapping_t>& mappings) {
            //            doNotOptimize(mappings);
            //            return true;
            std::cout << query_id << '(' << query_length << "): \n";
            for (const auto& mapping : mappings) {
                // if (mapping.target_id() == query_id || mapping.size() < 5) continue;

                std::cout.put('\t');
                if (mapping.is_rev()) std::cout.put('-');
                std::cout << mapping.target_id() << '(' << mapping.target_length() << "): ";
                for (auto& seed : mapping)
                    std::cout << seed.query_pos << '/' << seed.target_pos << ' ';
                std::cout.put('\n');
            }

            return true;
        };

    query q{index, on_read};
    for (auto& fin : fq_in)
        q.read_fastx(fin);
}

struct read_info_t
{

    uint16_t in_edges  = 0;
    uint16_t out_edges = 0;
    uint32_t in_seeds  = 0;
    uint32_t out_seeds = 0;

    size_t seeds() const { return in_seeds + out_seeds; }
    size_t edges() const { return in_edges + out_edges; }

    void add_in_edge(const seedlib::query::mapping_t& mapping)
    {
        in_seeds += mapping.size();
        in_edges += 1;
    }

    void add_out_edges(uint32_t seeds, uint16_t edges)
    {
        out_seeds += seeds;
        out_edges += edges;
    }

    friend std::ostream& operator<<(std::ostream& out, const read_info_t& that)
    {
        out << "in:(" << that.in_edges << ", " << that.in_seeds << ") out:(" << that.out_edges << ", " << that.out_seeds
            << ")";
        return out;
    }
};

void
query_index_star(const std::vector<std::string>& fq_in, const seedlib::index& index)
{
    using namespace seedlib;
    using read_id_t  = query::read_id_t;
    using read_pos_t = query::read_pos_t;

    std::vector<read_info_t> read_info;
    auto                     get_read_info = [&](read_id_t id) -> read_info_t& {
        if (unlikely(id >= read_info.size())) {
            size_t new_size = id + 1;
            new_size        = new_size + new_size / 2; // 1.5 growth factor
            read_info.resize(new_size);
        }
        assume(id < read_info.size(), "wtf");
        return read_info[id];
    };

    auto on_read = [&](read_id_t query_id, read_pos_t query_length, const std::vector<query::mapping_t>& mappings) {
        if (!mappings.empty()) {

            auto& query_info    = get_read_info(query_id);
            auto  prev_in_edges = query_info.in_edges;
            std::cout << query_id << "(sz:" << query_length << ", " << query_info << ") -> (";

            size_t total_seeds    = 0;
            size_t total_mappings = 0;
            for (const auto& mapping : mappings) {
                auto target_id = mapping.target_id();

                if (target_id == query_id) continue;
                total_mappings++;
                total_seeds += mapping.size();
                //                assert(mapping.size() >= 10, "wtf");

                if (prev_in_edges < 1) {
                    get_read_info(target_id).add_in_edge(mapping);
                    //                    std::cout << query_id << '\t' << target_id << '\t' << mapping.size() << '\n';
                } else {
                    //                    std::cout << query_id << '\t' << target_id << '\t' << mapping.size() << "c\n";
                }

                //                continue;

                std::cout.put('\t');
                if (mapping.is_rev()) std::cout.put('-');
                std::cout << mapping.target_id() << '(' << mapping.target_length() << "): ";
                for (auto& seed : mapping)
                    std::cout << seed.query_pos << '/' << seed.target_pos << ' ';
                std::cout.put('\n');
            }

            {
                auto& query_info = get_read_info(query_id);
                query_info.add_out_edges(total_seeds, total_mappings);

                std::cout << query_info << ")\n";
            }
        }

        auto& next_read_info = get_read_info(query_id + 1);

        if (next_read_info.in_edges < 1)
            return true;
        else {
            return true;
            std::cerr << query_id + 1 << "(skip: " << next_read_info << ")\n";
            return false;
        }
    };

    query q{index, on_read};
    for (auto& fin : fq_in)
        q.read_fastx(fin);
}

void
query_index(const std::vector<std::string>& fq_in, const std::string& index_name, bool star = true)
{

    seedlib::index index(index_name);

    Timer t;
    if (star)
        query_index_star(fq_in, index);
    else
        query_index_clique(fq_in, index);
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

    auto indexArg   = required("-i", "--index") & value("index", config.output);
    auto inputsArgs = values("fast(a|q) files", config.inputs);

    auto indexCmd = command("index").set(selected, mode::index);
    auto indexMode
      = (indexCmd,
         indexArg,
         option("--b1_len", "-1") & value("b1_len", config.b1) % "b1 length (int [1-8])",
         option("--b2_len", "-2") & value("b2_len", config.b2) % "b2 length (int [0-8], =0 use b1)",
         option("--b3_len", "-3") & value("b3_len", config.b3) % "b3 length (int [0-8], =0 use b1)",
         option("--complexity", "-c")
           & value("entropy (bits)", config.min_entropy)
               % "min entropy (in bits, float [0-4]), of the 2mers distribution in b1+b2+b3-mers",
         option("--downsample", "-d")
           & value("downsample", downsampling) % "inverse of fraction of positions retained (float [0-1])",
         option("--b1vsb1b2", "-r")
           & value("b1vsb1b2", config.b1_vs_b1b2_drop) % "0 to downsample by b1, 1 to downsample by b1b2 (float [0-1])",
         inputsArgs);

    auto statCmd  = command("stat").set(selected, mode::stat);
    auto statMode = (statCmd, indexArg);

    auto findCmd  = command("query").set(selected, mode::query);
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
