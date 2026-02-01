#include <sstream>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/search.hpp>

int main(int argc, char const* const* argv) {
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    seqan3::argument_parser parser{"fmindex_pigeon_search", argc, argv, seqan3::update_notifications::off};

    parser.info.author = "SeqAn-Team";
    parser.info.version = "1.0.0";

    auto index_path = std::filesystem::path{};
    parser.add_option(index_path, '\0', "index", "path to the query file");

    auto reference_file = std::filesystem::path{};
    parser.add_option(reference_file, '\0', "reference", "path to the reference file");

    auto query_file = std::filesystem::path{};
    parser.add_option(query_file, '\0', "query", "path to the query file");

    auto number_of_queries = size_t{100};
    parser.add_option(number_of_queries, '\0', "query_ct", "number of query, if not enough queries, these will be duplicated");

    auto number_of_errors = uint8_t{0};
    parser.add_option(number_of_errors, '\0', "errors", "number of allowed hamming distance errors");

    try {
         parser.parse();
    } catch (seqan3::argument_parser_error const& ext) {
        seqan3::debug_stream << "Parsing error. " << ext.what() << "\n";
        return EXIT_FAILURE;
    }

    // loading our files
    auto reference_stream = seqan3::sequence_file_input{reference_file};
    auto query_stream     = seqan3::sequence_file_input{query_file};

    // read reference into memory
    std::vector<std::vector<seqan3::dna5>> reference;
    for (auto& record : reference_stream) {
        reference.push_back(record.sequence());
    }

    // read query into memory
    std::vector<std::vector<seqan3::dna5>> queries;
    for (auto& record : query_stream) {
        queries.push_back(record.sequence());
    }

    // loading fm-index into memory
    using Index = decltype(seqan3::fm_index{std::vector<std::vector<seqan3::dna5>>{}}); // Some hack
    Index index; // construct fm-index
    {
        seqan3::debug_stream << "Loading 2FM-Index ... " << std::flush;
        std::ifstream is{index_path, std::ios::binary};
        cereal::BinaryInputArchive iarchive{is};
        iarchive(index);
        seqan3::debug_stream << "done\n";
    }

    // duplicate input until its large enough
    while (queries.size() < number_of_queries) {
        auto old_count = queries.size();
        queries.resize(2 * old_count);
        std::copy_n(queries.begin(), old_count, queries.begin() + old_count);
    }
    queries.resize(number_of_queries); // will reduce the amount of searches
    
    std::chrono::steady_clock::time_point end_load = std::chrono::steady_clock::now();
    seqan3::debug_stream << "Loading time " << std::chrono::duration_cast<std::chrono::milliseconds>(end_load - begin).count() << " ms\n";
    std::chrono::steady_clock::time_point begin_search std::chrono::steady_clock::now();

    seqan3::configuration const cfg = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{0}};
    //!TODO !ImplementMe use the seqan3::search to find a partial error free hit, verify the rest inside the text
    // Pseudo code (might be wrong):
    // for query in queries:
    //      parts[3] = cut_query(3, query);
    //      for p in {0, 1, 2}:
    //          for (pos in search(index, part[p]):
    //              if (verify(ref, query, pos +- ....something)):
    //                  std::cout << "found something\n"

    size_t parts = number_of_errors + 1;
    for (auto const & query : queries)
    {
        size_t part_length = query.size() / parts;
        for (size_t p = 0; p < parts; ++p)
        {
            // Teilstring bestimmen
            auto part_start = query.begin() + p * part_length;
            auto part_end = (p == parts - 1) ? query.end() : part_start + part_length;
            std::vector<seqan3::dna5> query_part(part_start, part_end);

            // Suche Teilstück mit 0 Fehlern
            for (auto && result : seqan3::search(query_part, index, cfg))
            {
                size_t ref_id = result.reference_id();
                size_t pos = result.reference_begin_position();

                // Startposition der gesamten Query im Referenztext berechnen
                size_t query_start = (pos > p * part_length) ? pos - p * part_length : 0;

                // Überprüfen, ob Query an dieser Position mit <= number_of_errors passt
                if (query_start + query.size() <= reference[ref_id].size())
                {
                    size_t errors = 0;
                    for (size_t i = 0; i < query.size(); ++i){
                        if (query[i] != reference[ref_id][query_start + i]){
                            ++errors;
                        } 
                    }
                    // if (errors <= number_of_errors){
                    //     std::cout << "Treffer in Referenz " << ref_id << " an Position " << query_start << " mit " << errors << " Fehlern\n";
                    //     }
                    
                }
            }
        }
    }
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    seqan3::debug_stream << "Search time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin_search).count() << " ms\n";

    return 0;
}
