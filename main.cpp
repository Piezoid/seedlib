#include <iostream>
#include "clipp.h"

#include "common.hpp"
#include "indexer.hpp"

using namespace std;
using namespace seedlib;

int
main(int argc, char* argv[])
{
    using namespace clipp;
    // variables storing the parsing result; initialized with their default values
    enum class mode { index, find, help };
    mode      selected = mode::help;
    file_list input;
    string    index, out;
    bool      split = false, progr = false;

    auto indexArg = required("-d", "--index") & value("index", index);

    auto indexMode = (command("index").set(selected, mode::index),
                      indexArg,
                      option("--progress", "-p").set(progr) % "show progress",
                      values("fastq files", input));

    auto findMode
      = (command("find").set(selected, mode::find),
         values("infile", input),
         indexArg,
         (option("-o", "--output") & value("outfile", out)) % "write to file instead of stdout",
         (option("-split").set(split, true) | option("-nosplit").set(split, false)) % "(do not) split output");

    auto cli = ((indexMode | findMode | command("help").set(selected, mode::help)),
                option("-v", "--version").call([] { cout << "version 0.0.1\n\n"; }).doc("show version"));

    if (parse(argc, argv, cli)) {
        switch (selected) {
            case mode::index: build_index(input, index); break;
            case mode::find: break;
            case mode::help: cout << make_man_page(cli, argv[0]); break;
        }
    } else {
        cout << usage_lines(cli, argv[0]) << '\n';
    }
}
