#include <divsufsort.h>
#include <sstream>
#include <iostream>
#include <chrono>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/search.hpp>

int compare_suffix(std::vector<seqan3::dna5> const & text,
                   size_t sa_pos,
                   std::vector<seqan3::dna5> const & pattern){
    size_t i = 0;
    while (i < pattern.size() && sa_pos + i < text.size())
    {
        // if text smaller than pattern -> return -1
        if (text[sa_pos + i] < pattern[i]) return -1; 
        // if text larger than pattern -> return -1
        if (text[sa_pos + i] > pattern[i]) return 1;
        // if text equal to pattern -> next character
        ++i;
    }

    if (i == pattern.size())
        return 0; // full match

    return -1; // if suffix shorter than pattern
}

size_t left_border(std::vector<seqan3::dna5> const & text,
               std::vector<saidx_t> const & sa,
               std::vector<seqan3::dna5> const & pattern){
    size_t L = 0, R = sa.size();
    while (L < R)
    {
        size_t M = (L + R) / 2;
        if (compare_suffix(text, sa[M], pattern) >= 0)
            R = M;
        else
            L = M + 1;
    }
    return L;
}

size_t right_border(std::vector<seqan3::dna5> const & text,
               std::vector<saidx_t> const & sa,
               std::vector<seqan3::dna5> const & pattern){
    size_t L = 0, R = sa.size();
    while (L < R)
    {
        size_t M = (L + R) / 2;
        if (compare_suffix(text, sa[M], pattern) > 0)
            R = M;
        else
            L = M + 1;
    }
    return L;
}

int main(int argc, char const* const* argv) {
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    seqan3::argument_parser parser{"suffixarray_search", argc, argv, seqan3::update_notifications::off};

    parser.info.author = "SeqAn-Team";
    parser.info.version = "1.0.0";

    auto reference_file = std::filesystem::path{};
    parser.add_option(reference_file, '\0', "reference", "path to the reference file");

    auto query_file = std::filesystem::path{};
    parser.add_option(query_file, '\0', "query", "path to the query file");

    auto number_of_queries = size_t{100};
    parser.add_option(number_of_queries, '\0', "query_ct", "number of query, if not enough queries, these will be duplicated");

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
    // Attention: we are concatenating all sequences into one big combined sequence
    //            this is done to simplify the implementation of suffix_arrays
    std::vector<seqan3::dna5> reference;
    for (auto& record : reference_stream) {
        auto r = record.sequence();
        reference.insert(reference.end(), r.begin(), r.end());
    }

    // read query into memory
    std::vector<std::vector<seqan3::dna5>> queries;
    for (auto& record : query_stream) {
        queries.push_back(record.sequence());
    }

    // duplicate input until its large enough
    while (queries.size() < number_of_queries) {
        auto old_count = queries.size();
        queries.resize(2 * old_count);
        std::copy_n(queries.begin(), old_count, queries.begin() + old_count);
    }
    queries.resize(number_of_queries); // will reduce the amount of searches

    // Array that should hold the future suffix array
    std::vector<saidx_t> suffixarray;
    suffixarray.resize(reference.size()); // resizing the array, so it can hold the complete SA

    //!TODO !ImplementMe implement suffix array sort
    //Hint, if can use libdivsufsort (already integrated in this repo)
    //      https://github.com/y-256/libdivsufsort
    //      To make the `reference` compatible with libdivsufsort you can simply
    //      cast it by calling:
    //      `sauchar_t const* str = reinterpret_cast<sauchar_t const*>(reference.data());`

    sauchar_t const* str = reinterpret_cast<sauchar_t const*>(reference.data());
    divsufsort(str, suffixarray.data(), reference.size());
    // seqan3::debug_stream << "test" << suffixarray << "/n";
    // Calculates left boarder
    std::chrono::steady_clock::time_point end_SA = std::chrono::steady_clock::now();
    seqan3::debug_stream << "Elapsed time for SA build-up: " << std::chrono::duration_cast<std::chrono::seconds>(end_SA - begin).count() << " s\n";

    std::chrono::steady_clock::time_point begin2 = std::chrono::steady_clock::now();
    for (auto& q : queries) {
        //!TODO !ImplementMe apply binary search and find q  in reference using binary search on `suffixarray`
        // You can choose if you want to use binary search based on left_border

        size_t l = left_border(reference, suffixarray, q);
        size_t r = right_border(reference, suffixarray, q);
        std::vector<size_t> occurences;
        for (size_t i = l; i < r; ++i) {
            occurences.push_back(suffixarray[i]);
        }
        // seqan3::debug_stream << "Query" << q << "found at positions: " << occurences << "\n";
    }
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    seqan3::debug_stream << "Elapsed time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin2).count() << " ms\n";

    return 0;
}
