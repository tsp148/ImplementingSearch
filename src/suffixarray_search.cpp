#include <divsufsort.h>
#include <sstream>
#include <iostream>
#include <chrono>
#include <vector>
#include <algorithm>
#include <string>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/search.hpp>

std::string text1;
bool rb_check;

struct functor
{   
    bool operator()(const uint32_t a, const uint32_t b)const{
    //Abfrage welches Suffix größer ist --> Anzahl Vergleiche = Anzahl Buchstaben des kürzeren Suffixes
    uint32_t j = text1.size()-std::max(a,b); // a ist größerer Index also ist das dazugehörige Suffix weiter hinten im text und damit kürzer

    j = text1.size()-std::max(a,b);

    for (uint32_t i = 0; ((text1[a+i] <= text1[b+i]) && (j>0) && a>b) || text1[a+i] < text1[b+i]; ++i){
        if (text1[a+i] < text1[b+i])
            return true;
        --j;
    }
    if (j==0) return true;
    return false;
    }
}sort_suffix_indices;

void construct(std::vector<saidx_t>& sa, const std::string& t){
    sa.clear();
    if(!t.empty()){
    text1 = t;
    for (uint32_t i = 0; i < t.size(); ++i){
        sa.push_back(i);
    }
    sort(sa.begin(), sa.end(), sort_suffix_indices);
    }
}

bool is_smaller (const std::string& pat, const uint32_t sa_value, const bool b, const uint32_t mlr){
    uint32_t j = std::min(pat.size(), text1.size()-sa_value) - mlr;
    //std::cout << "Startindex: " << start_index << "\n";
    for (uint32_t i = mlr; (pat[i] <= text1[sa_value+i]) && (j > 0) ; ++i){
        if (pat[i] < text1[sa_value+i]){
            return true;
        }
        else 
            j--;
    }
    if (b == false && pat.size() <= text1.size()-sa_value){
        if (j==0) 
        return true;
    }
    return false;
}

uint32_t lcp_value(const std::string& pat,const uint32_t n){
    uint32_t lcp_value = 0; // counter, der erhöht wird wenn sich die Suffixe um einen zusätzlichen Buchstaben gleichen
    uint32_t j = std::min(pat.size(),text1.size()-n);

    for (unsigned i = 0; pat[i] == text1[n+i] && j > 0; ++i){
        lcp_value += 1;
        --j;
    }
    return lcp_value;
}

uint32_t left_border (const std::string& pat, const std::vector<saidx_t>& sa){
    

    uint32_t lb, M, L, R, l, r;

    if (is_smaller(pat, sa[0], false, 0) == true)
        lb = 0;
    else if (is_smaller(pat, sa[sa.size()-1], false, 0)!= true)
        lb = sa.size();
    else{
        L = 0;
        R = sa.size()-1;

        while (R-L>1){
            M = (L+R)/2;
            l = lcp_value(pat, sa[L]); //lcp-Wert von Pattern und linker Grenze
            r = lcp_value(pat, sa[R]); //lcp-Wert von Pattern und rechter Grenze
            if (is_smaller(pat, sa[M], false, std::min(l,r)) == true)
                R = M;
            else 
                L = M;
        }
        lb = R;
    }
    return lb;
}

int right_border (const std::string& pat, const std::vector<saidx_t>& sa){
    uint32_t rb;
    uint32_t M, L, R, l, r;

    if (is_smaller(pat, sa[0], true,0) == true){
        rb_check = true;
        rb = 0;
    }
    else if (is_smaller(pat, sa[sa.size()-1], true,0) != true)
        rb = sa.size()-1;
    else{
        L = 0;
        R = sa.size()-1;
        while (R-L>1){
            M = (L+R)/2;
            
            l = lcp_value(pat, sa[L]);
            r = lcp_value(pat, sa[R]);
            
            if (is_smaller(pat, sa[M], true, std::min(l,r)) == true)
                R = M;
            else
                L = M;
        }
        rb = L;
    }
    return rb;
}

void find(const std::string& query, const std::vector<saidx_t>& sa, const std::string& text, std::vector<uint32_t>& hits){
    text1 = text;
    hits.clear();
    if (query != " " && !query.empty() && !text.empty()){
    uint32_t lb = left_border(query, sa);
    uint32_t rb = right_border(query, sa);
    if (rb_check != true)
        if (lb <= rb){
            for(uint32_t i = lb; i <= rb; ++i){
            	hits.push_back(sa[i]);
	    }    
	    sort(hits.begin(), hits.end());
    	}
    }
}

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


    // from here
    /*
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
    */
    // to here

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
        hits = std::vector<uint32_t>{};
        find(q, suffixarray, reference, hits);
        seqan3::debug_stream << "Query" << q << "found at positions: \n" << hits << "\n";

        // size_t l = left_border(reference, q);
        // size_t r = right_border(reference, q);
        // seqan3::debug_stream << "l: " << l << "\n";
        // seqan3::debug_stream << "r: " << r << "\n";
        // std::vector<size_t> occurences;
        // for (size_t i = l; i < r; ++i) {
        //     occurences.push_back(suffixarray[i]);
            // Prüfe, ob das Suffix wirklich mit der Query übereinstimmt (optional, falls exakte Matches gewünscht)
            // auto suf_it = reference.begin() + suffixarray[i];
            // if (std::equal(q.begin(), q.end(), suf_it, reference.end())) {
            //     occurences.push_back(suffixarray[i]);
            // }
        }
        // Optional: seqan3::debug_stream << "Query found at positions: " << occurences << "\n";
        // seqan3::debug_stream << "Query" << q << "found at positions: " << occurences << "\n";
    // }
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    seqan3::debug_stream << "Elapsed time: " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " s\n";

    return 0;
}
