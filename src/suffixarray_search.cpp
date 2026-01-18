#include <divsufsort.h>
#include <sstream>
#include <iostream>
#include <chrono>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/search.hpp>

// Calculates left boarder
// auto left_border = [&](std::vector<seqan3::dna5> const & reference, const std::vector<seqan3::dna5>& query, std::vector<saidx_t> const&  suffixarray) -> size_t {
//     size_t left = 0, right = reference.size();
//     while (left < right) {
//         size_t mid = left + (right - left) / 2;
//         auto suf_it = reference.begin() + suffixarray[mid];
//         if (std::lexicographical_compare(suf_it, reference.end(), query.begin(), query.end())) {
//             left = mid + 1;
//         } else {
//             right = mid;
//         }
//     }
//     return left;
// };

// // Calculates right boarder
// auto right_border = [&](std::vector<seqan3::dna5> const & reference, const std::vector<seqan3::dna5>& query, std::vector<saidx_t> const & suffixarray) -> size_t {
//     size_t left = 0, right = reference.size();
//     while (left < right) {
//         size_t mid = left + (right - left) / 2;
//         auto suf_it = reference.begin() + suffixarray[mid];
//         if (std::lexicographical_compare(query.begin(), query.end(), suf_it, reference.end())) {
//             right = mid;
//         } else {
//             left = mid + 1;
//         }
//     }
//     return left;
// };

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

    // auto compare_suffix(std::vector<seqan3::dna5> const & reference, size_t const & mid, std::vector<seqan3::dna5> const & query){
    //     // seq_ref = reference[]
    //     // if ref from pos mid not lexicographically smaller than query, then query is smaller than suffix
    //     if(!std::lexicographical_compare(reference.begin()+sufix_array[mid], reference.begin()+sufix_array[mid]+query.size())){
    //         return true;
    //     }

    //     if(std::equal(ref.begin() + i, ref.begin() + i + query.size(), query.begin())){
    //         occurences.push_back(i);
    //     }
    // }

    auto left_border = [&](std::vector<seqan3::dna5> const & reference, const std::vector<seqan3::dna5>& query) -> size_t {
        if(!std::lexicographical_compare(reference.begin()+suffixarray[0], reference.begin()+suffixarray[0]+query.size()-1, query.begin(), query.end())){
            return 0;
        }
        else if(std::lexicographical_compare(reference.begin()+suffixarray[suffixarray.size()-1], reference.begin()+suffixarray[suffixarray.size()-1]+query.size()-1, query.begin(), query.end())){
            return suffixarray.size();
        }
        else{
            size_t left = 0, right = suffixarray.size()-1;
            while((right - left) > 1){
                size_t mid = std::ceil((left + right) / 2);
                auto suf_it = reference.begin() + suffixarray[mid];
                if(!std::lexicographical_compare(suf_it, reference.begin()+suffixarray[mid]+query.size()-1, query.begin(), query.end())){
                    right = mid;
                }
                else{
                    left = mid;
                }
            }
            return right;
        }
    };
    
    auto right_border = [&](std::vector<seqan3::dna5> const & reference, const std::vector<seqan3::dna5>& query) -> size_t {
        if(std::lexicographical_compare(query.begin(), query.end(), reference.begin()+suffixarray[0], reference.begin()+suffixarray[0]+query.size()-1)){
            return 0;
        }
        // else if(!std::lexicographical_compare(reference.begin()+suffixarray[suffixarray.size()-1], reference.begin()+suffixarray[suffixarray.size()-1]+query.size()-1, query.begin(), query.end())){
        //     return suffixarray.size();
        // }
        else{
            size_t left = 0, right = suffixarray.size()-1;
            while((right - left) > 1){
                size_t mid = std::ceil((left + right) / 2);
                auto suf_it = reference.begin() + suffixarray[mid];
                if(!std::lexicographical_compare(query.begin(), query.end(), suf_it, reference.begin()+suffixarray[mid]+query.size()-1)){
                    left = mid;
                }
                else{
                    right = mid;
                }
            }
            return left;
        }
    };

    // auto left_border = [&](std::vector<seqan3::dna5> const & reference, const std::vector<seqan3::dna5>& query) -> size_t {
    //     size_t left = 0, right = reference.size();
    //     while (left < right) {
    //         size_t mid = left + (right - left) / 2;
    //         auto suf_it = reference.begin() + suffixarray[mid];
    //         if (std::lexicographical_compare(suf_it, reference.end(), query.begin(), query.end())) {
    //             left = mid + 1;
    //         } else {
    //             right = mid;
    //         }
    //     }
    //     return left;
    // };

    // Calculates right boarder
    // auto right_border(std::vector<seqan3::dna5> const & reference, const std::vector<seqan3::dna5>& query) -> size_t {
    //     size_t left = 0, right = reference.size();
    //     while (left < right) {
    //         size_t mid = left + (right - left) / 2;
    //         auto suf_it = reference.begin() + suffixarray[mid];
    //         if (std::lexicographical_compare(query.begin(), query.end(), suf_it, reference.end())) {
    //             right = mid;
    //         } else {
    //             left = mid + 1;
    //         }
    //     }
    //     return left;
    // };


    for (auto& q : queries) {
        //!TODO !ImplementMe apply binary search and find q  in reference using binary search on `suffixarray`
        // You can choose if you want to use binary search based on left_border
        size_t l = left_border(reference, q);
        size_t r = right_border(reference, q);
        seqan3::debug_stream << "l: " << l << "\n";
        seqan3::debug_stream << "r: " << r << "\n";
        std::vector<size_t> occurences;
        for (size_t i = l; i < r; ++i) {
            occurences.push_back(suffixarray[i]);
            // Prüfe, ob das Suffix wirklich mit der Query übereinstimmt (optional, falls exakte Matches gewünscht)
            // auto suf_it = reference.begin() + suffixarray[i];
            // if (std::equal(q.begin(), q.end(), suf_it, reference.end())) {
            //     occurences.push_back(suffixarray[i]);
            // }
        }
        // Optional: seqan3::debug_stream << "Query found at positions: " << occurences << "\n";
        seqan3::debug_stream << "Query" << q << "found at positions: " << occurences << "\n";
    }
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    seqan3::debug_stream << "Elapsed time: " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " s\n";

    return 0;
}
