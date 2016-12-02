#ifndef QUERY_H
#define QUERY_H

#include <set>
#include <string>
#include <vector>
#include <list>
#include <future>
#include <iostream>

#include "BloomTree.h"
#include "Kmers.h"
#include "BF.h"

extern float QUERY_THRESHOLD;

// a QueryInfo is
//   - a nt string, the "query"
//   - the set of kmers in the query string (*not* canonicalized)
//   - the number of nodes visited until the query was resolved
//   - a vector of nodes that "hit" the query
class QueryInfo {
public:
    QueryInfo(const std::string & q) : query(q), query_kmers(kmers_in_string(q)), nodes_visited(0) {}
    QueryInfo(const std::string & q, const std::string & w) {
        nodes_visited = 0;
        query = q;
        query_kmers = kmers_in_string(q);
        std::vector<std::string> fields;
        SplitString(w, ' ', fields);
        unsigned n = 0;
        for (const auto & w : fields) {
            if (w!="") {
                // If string w has invalid letters after numbers, composite_string will return pointer (else null)
                // std::string::size_type not working here?
                std::size_t composite_string;
                //std::cerr << typeid(w).name() << std::endl;
                try{
                    float value = std::stof(w, &composite_string);
                                if (composite_string != w.size()) {
                                        std::cerr << "Invalid weight \'" << w << "\' at position " << n << std::endl;
                                        exit(3);
                                }
                                weight.emplace_back(value); //Currently zero error handling here!
    
                }
                catch(...){
                    std::cerr << "Invalid weight \'" << w << "\' at position " << n << std::endl;
                    exit(3);
                }
            }
            n++;
        }
    }
    
    // added by chen sun, 10/17/2016
    QueryInfo(): nodes_visited(0){
        query = "Experiment";
    }

    ~QueryInfo() {}

    // added by chen sun, 10/17/2016
    void InsertKmer(const std::string & q){
        auto k = jellyfish::mer_dna::k();
        if(q.size() > k){
            query_kmers.insert(jellyfish::mer_dna(q.substr(0, k)));
        }else{
            query_kmers.insert(jellyfish::mer_dna(q));
        }
    }
   
    std::string query;
    std::set<jellyfish::mer_dna> query_kmers;
    std::vector<const BloomTree*> matching;
    std::vector<float> weight;
    int nodes_visited;
};

// a QueryState is
//   - all the stuff from QueryInfo
//   - a vector of integers, kmers converted to positions in the bf
//   - the number of kmers that have previously passed/failed (and have been
//     removed from the vector of positions)
//   - the number of kmers needed to resolve the query as pass or fail
class QueryState {
public:
    QueryState(BF* bf, QueryInfo* q_info) {
        q = q_info;
        for (const auto & m : q_info->query_kmers) {
            kmer_positions.emplace_back(bf->kmer_position(m));
        }
        pass_thresh = fail_thresh = 0;
        num_passed = num_failed = 0;
    }
    QueryState(QueryState* q_state) {
        q = q_state->q;
        pass_thresh = q_state->pass_thresh;
        fail_thresh = q_state->fail_thresh;
        num_passed = q_state->num_passed;
        num_failed = q_state->num_failed;
        // kmer_positions will be intentionally empty
    }
    ~QueryState() {}

    QueryInfo* q;
    std::vector<uint64_t> kmer_positions;
    unsigned pass_thresh;
    unsigned fail_thresh;
    unsigned num_passed;
    unsigned num_failed;
};

using QuerySet = std::list<QueryInfo*>;
using QueryStateSet = std::list<QueryState*>;

void query_from_file(BloomTree* root, const std::string & fn, std::ostream & o);
void batch_query_from_file(BloomTree* root, const std::string & fn, std::ostream & o);
void batch_query_by_reduction_from_file(BloomTree* root, const std::string & fn, std::ostream & o);
void original_batch_query_from_file(BloomTree* root, const std::string & fn, std::ostream & o);
void batch_weightedquery_from_file(BloomTree* root, const std::string & fn, const std::string & wf, std::ostream & o); 
int  query_string(BloomTree* root, const std::string & q, std::vector<BloomTree*> & out);
int  query(BloomTree* root, const std::set<jellyfish::mer_dna> & q, std::vector<BloomTree*> & out);
void check_bt(BloomTree* root);
void intersect_bt(BloomTree* root);
void draw_bt(BloomTree* root, std::string outfile);
void compress_bt_rrr(BloomTree* root);
void compress_bt_roar(BloomTree* root);
void compress_bt_roar_parallel(BloomTree* root, std::vector<std::future<void> > *F = nullptr);
void uncompress_bt(BloomTree* root);

void leaf_query_from_file(BloomTree* root, const std::string & fn, std::ostream & o);

//add by chen sun, 10/17/2016
void BatchQueryByKmerSet(
    BloomTree* root,
    const std::string & fn,
    std::ostream & o);
#endif
