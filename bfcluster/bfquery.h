// contact: chensun@cse.psu.edu

#ifndef bfquery_H		// (prevent multiple inclusion)
#define bfquery_H

#include "util.h"
#include "utilities.h"
#include <vector>
#include <unordered_map>
#include <unordered_set>

#include "BF.h"
#include "BloomTree.h"
#include "Build.h"
//#include "Count.h"
#include "Query.h"



//----
// RoarQueryInfo
//  similar with QueryInfo in Query.h
//  Instead of storing kmers, store a RoarBF of query sequence
class RoarQueryInfo {
public:
    RoarQueryInfo(const std::string & q, HashPair hp, int nh, uint64_t bf_size, bool heuristic_match = true) : query(q) {
    	// generate a bit vector;
        roarbf = new RoarBF("tmp.query.bf.bv.roar", hp, nh, bf_size);
    	roarbf->roar = roaring_bitmap_create();
        kmer_number = 0;
        ones_number = 0;
        if(q == "") return;
        InsertSingleSequence(q, hp, nh, bf_size, heuristic_match);
        matching_closed = false;
        //inner_heuristic_match = heuristic_match;
        repeat_bf = NULL;
        heuristic_match_indicator = heuristic_match;
    }

    RoarQueryInfo(RoarBF* new_roarbf, bool heuristic_match){
        // not complete copy
        query = "[Warning] intermedia queryinfo";
        roarbf = new_roarbf;
        matching_closed = true;
        repeat_bf = NULL;
        heuristic_match_indicator = heuristic_match;
        // if use this construct function, it means to use split function
        // in this condition, matching is stored in RoarQueryState
        // matching function is hence closed
        // if closed, can not output matching results directly
    }

    ~RoarQueryInfo() {
    	if(roarbf != NULL){
    		delete roarbf; // call destructor to free roarbf->roar
    		roarbf = NULL;
    	}
        if(!heuristic_match_indicator)
        {        
            if(repeat_bf != NULL){
                delete repeat_bf;
                repeat_bf = NULL;
            }
        }
    }
   	
    // clear roarbf, store into memory for output
    // this is still dangerous if large queryies
    // suggest output one by one
    void ClearBf(){
        if(roarbf != NULL){
            delete roarbf;
            roarbf = NULL;
        }
    }

    void InsertSequence(const std::string & q, HashPair hp, int nh, uint64_t bf_size, bool heuristic_match = true)
    //----
    // Warning
    //      [Warning] updating ones_number is on caller's duty, because this might be called many times
    //----
    {
        std::set<jellyfish::mer_dna> query_kmers = kmers_in_string(q);
        kmer_number += query_kmers.size();
        // for each kmer, calculate its hash value and set corresponding bit as 1
        for(jellyfish::mer_dna mer : query_kmers){
            // BF::contians(const jellyfish::mer_dna & m) const: how to calculate hash value for a given kmer 
            mer.canonicalize();
            uint64_t h0 = hp.m1.times(mer);
            const size_t pos = h0 % bf_size;
            if(heuristic_match){
                roaring_bitmap_add(roarbf->roar, pos);
            }else{
                if(roarbf->operator[](pos) == 1){
                    if(multiple_hits.find(pos) != multiple_hits.end()){
                        multiple_hits[pos]++;
                    }else{
                        multiple_hits[pos] = 2;
                    }
                }else{
                    roaring_bitmap_add(roarbf->roar, pos);
                }
            }
        }
        //----
        // Warning
        //      [Warning] updating ones_number is on caller's duty.
        //----
    }

    void InsertSingleSequence(const std::string & q, HashPair hp, int nh, uint64_t bf_size, bool heuristic_match = true)
    {
        // [Warning] kmer_number will be updated, ones_number will not
        InsertSequence(q, hp, nh, bf_size, heuristic_match);
        if(heuristic_match)
        // update ones_number by hand
        {
            ones_number = (int)(roarbf->count_ones());
        }
        //else no need to use ones_number
    }

    void InsertSequenceSet(std::vector<std::string> & kmer_set, HashPair hp, int nh, uint64_t bf_size, bool heuristic_match = true)
    {
        for(auto kmer: kmer_set){
            InsertSequence(kmer, hp, nh, bf_size, heuristic_match);
        }
        if(heuristic_match){
            ones_number = (int)(roarbf->count_ones());
        }
    }

    void CheckKmer(bool heuristic_match)
    // this is for unit test
    {
        std::cerr << "kmer_number: " << kmer_number << std::endl;
        int counts = (int)(roarbf->count_ones());
        if(ones_number != counts) ones_number = counts;
        if(! heuristic_match){
            for(auto it = multiple_hits.begin(); it != multiple_hits.end(); ++it){
                counts += it->second - 1;
            }
            if(kmer_number != counts){
                std::cerr << "[RoarQueryInfo] Error: kmer number do not match, can not pass Check Kmer" << std::endl;
            }
        }
        std::cerr << "counts: " << counts << std::endl;
    }

    void CreateRepeatBF(HashPair hp, int nh, uint64_t bf_size, bool heuristic_match){
        if(heuristic_match){
            repeat_bf = NULL;
        }else{
            repeat_bf = new RoarBF("repeat.kmer.bf.bv.roar", hp, nh, bf_size);
            repeat_bf->roar = roaring_bitmap_create();
            for(auto it = multiple_hits.begin(); it != multiple_hits.end(); ++it){
                roaring_bitmap_add(repeat_bf->roar, it->first);
            }
        }
    }

    std::string query;
    int kmer_number;

    // in heuristic matching, if InsertSequence is called, you should update ones_number on your own
    int ones_number;
    
    int nodes_visited;
    //std::set<jellyfish::mer_dna> query_kmers;
    bool matching_closed;
    std::vector<const BloomTree*> matching;
    std::unordered_map<size_t, int> multiple_hits;
    //std::vector<float> weight;
    RoarBF* roarbf;
    RoarBF* repeat_bf;
    bool heuristic_match_indicator;

};

// a QueryState is
//   - all the stuff from QueryInfo
//   - a vector of integers, kmers converted to positions in the bf
//   - the number of kmers that have previously passed/failed (and have been
//     removed from the vector of positions)
//   - the number of kmers needed to resolve the query as pass or fail
class RoarQueryState {
public:
    RoarQueryState(RoarQueryInfo* q_info, unsigned _pass_thresh, unsigned _fail_thresh): pass_thresh(_pass_thresh), fail_thresh(_fail_thresh)
    // [Warning] It's creator's duty to set pass_thresh and fail_thresh
    {
        q = q_info;
        num_passed = num_failed = 0;
    }

    RoarQueryState(RoarQueryState* q_state): pass_thresh(q_state->pass_thresh), fail_thresh(q_state->fail_thresh) {
        q = NULL;
        id = q_state->id; 
        num_passed = q_state->num_passed;
        num_failed = q_state->num_failed;
    }


    int IntersectWithSome(RoarQueryState* child_query_state, RoarBF* some_roarbf, bool heuristic_match = true){
        // update roarbf with new intersection bit vector
        RoarBF* intersect_bf = dynamic_cast<RoarBF*>(some_roarbf->intersection_with(q->roarbf->filename, q->roarbf));
        RoarQueryInfo* child_query_info = new RoarQueryInfo(intersect_bf, heuristic_match);
        if(child_query_state->q == NULL){
            child_query_state->q = child_query_info;
        }else{
            std::cerr << "[RoarQueryState] intersect into an existing bloom filter" << std::endl;
            return -1;
        }

        int result = 0;
        if(heuristic_match){
            child_query_info->ones_number = (int)(intersect_bf->count_ones());
            result = q->ones_number - child_query_info->ones_number;
        }else{
            // update kmer number
            //std::cout << "calculate kmer match in some node: " << some_roarbf->filename << std::endl;
            child_query_info->repeat_bf = dynamic_cast<RoarBF*>(some_roarbf->intersection_with(q->repeat_bf->filename, q->repeat_bf));
            //std::cout << q->repeat_bf->count_ones() << "," << child_query_info->repeat_bf->count_ones() << std::endl;
            child_query_info->kmer_number = (int)(intersect_bf->count_ones());
            // delete while iterate
            
            for(auto it = q->multiple_hits.begin(); it != q->multiple_hits.end(); it++){
                const size_t pos = it->first;
                //std::cout << "check position: " << pos << std::endl; 
                if(child_query_info->repeat_bf->contains(pos)){
                    child_query_info->kmer_number += it->second-1;
                    child_query_info->multiple_hits[pos] = it->second;
                }
            }
            //std::cout << "original_kmer_number:" << original_kmer_number << " | " << "kmer_number:" << kmer_number << std::endl;
            result = q->kmer_number - child_query_info->kmer_number;
            //std::cout << "finish calculate kmer match in some node: " << some_roarbf->filename << std::endl;
        }
        return result;
    }

    int MatchAll(RoarBF* all_roarbf, bool heuristic_match = true){
        if(q == NULL) return -1;
        RoarBF* intersect_bf = dynamic_cast<RoarBF*>(all_roarbf->intersection_with("tmp.interset.bf.bv.roar", q->roarbf));
        
        int counts = intersect_bf->count_ones();
        if(!heuristic_match){
            // copied directly
            //std::cout << "calculate kmer match in all node: " << all_roarbf->filename << std::endl;

            RoarBF* intersect_repeat_bf = dynamic_cast<RoarBF*>(all_roarbf->intersection_with(q->repeat_bf->filename, q->repeat_bf));
            //std::cout << q->repeat_bf->count_ones() << "," << intersect_repeat_bf->count_ones() << std::endl;
            for(auto it = q->multiple_hits.begin(); it != q->multiple_hits.end(); it++){
                const size_t pos = it->first;
                if(all_roarbf->contains(pos)){
                    counts += it->second - 1;
                }
            }
            delete intersect_repeat_bf;
            //std::cout << "finish calculate kmer match in all node: " << all_roarbf->filename << std::endl;
        }
        delete intersect_bf;
        return counts;
    }

    RoarQueryState* QueryAllsomeNode(RoarBF* all_roarbf, RoarBF* some_roarbf, bool heuristic_match = true)
    // generate a new RoartQueryState with the same id
    // how to update num_passed and num_failed in all some tree
    // num_passed += MatchAll
    // num_failed += IntersectWithSome - MatchAll
    {
        RoarQueryState* child_query_state = new RoarQueryState(this);
        // child_query_state->q is still NULL, can not call MatchAll
        int counts = MatchAll(all_roarbf, heuristic_match);
        assert(counts >= 0);
        child_query_state->num_passed += counts;
        if(child_query_state->num_passed >= child_query_state->pass_thresh){
            // no need to intersect with some
            return child_query_state;
        }
        int tmp_fail = IntersectWithSome(child_query_state, some_roarbf, heuristic_match) - counts;
        assert(tmp_fail >= 0);
        child_query_state->num_failed += tmp_fail;
        return child_query_state;
    }

    ~RoarQueryState() {
        if(q != NULL){
            delete q;
            q = NULL;
        }
    }

    int id; // initialize as the index in query_state_set

	RoarQueryInfo* q; // this pointer points to the top of stack
    
    //once set, do not change
    const unsigned pass_thresh;
    const unsigned fail_thresh;

    unsigned num_passed;
    unsigned num_failed;
};

// because RoarQueryInfo stores the full roarbf, so its better not use RoarQuerySet
using RoarQuerySet = std::list<RoarQueryInfo*>;
using RoarQueryStateSet = std::vector<RoarQueryState*>;

//----
// class BFQuery 
//	rewrite Count.cc, generate bf from query
//	rewrite Query.cc, query not by kmer but by bf
class BFQuery{
public:
	//---- constructor & destructor
	BFQuery(double threshold = 0.9);

	~BFQuery();

	//---- public member variables

	//---- public member functions

    // ----single file query
	int CheckBfSize(std::string bf_filename, HashPair hp, int nh);
	
	int CheckRoarBfSize(std::string, HashPair, int);

    void PrintQueryInfo(RoarQueryInfo* q, std::ostream & out);

    int HeuristicQuery(BloomTree* root, RoarQueryInfo* q);

    bool HeuristicQueryPass(BloomTree* root, RoarQueryInfo* q);
	
	int Query(BloomTree* root, RoarQueryInfo* q);

    int Query(BloomTree* root, RoarBF* roarbf, std::vector<BloomTree*> & out);

    bool QueryPass(BloomTree* root, RoarQueryInfo* q);

    bool QueryPass(BloomTree* root, RoarBF * roarbf);

    // ----batch querying

    void BatchQuery(BloomTree* root, RoarQuerySet& qs, bool heuristic_match = true);

    void BatchQueryFromFile(const std::string& bloomtree_filename, 
                        const std::string& query_filename, 
                        const std::string& hash_filename,
                        const std::string& output_filename,
                        bool heuristic_match = true);

    void PrintQuerySet(const RoarQuerySet & qs, std::ostream & out) ;

	std::vector<std::string> ReadQueries(std::string query_filename);

    void SequentialQuery(const std::string& bloomtree_filename, 
                        const std::string& query_filename, 
                        const std::string& hash_filename,
                        const std::string& output_filename,
                        bool heuristic_match = true);

    // ---- experiment querying

    void ExperimentQueryFromFile(const std::string& bloomtree_filename, 
                        const std::string& query_filename, 
                        const std::string& hash_filename,
                        const std::string& output_filename,
                        bool heuristic_match = true);

    RoarQueryInfo* ReadExperiment(std::string query_filename, 
                                HashPair* hash_pair, 
                                int number_of_hashes, 
                                uint64_t bf_size, 
                                bool heuristic_match = true);

    // ---- split allsome querying

    //RoarQueryState* SplitQueryByStatus(BloomTree* root, RoarQueryState* q_state, bool heuristic_match = true);

    // already decide that q_state matches all leaves, iterate all leaves into matching_result
    void SplitBatchQueryAllLeavesMatch(BloomTree* root, int q_state_id);

    void SplitBatchQuery(BloomTree* root, RoarQueryStateSet active_queries, bool heuristic_match = true);

    void SplitBatchQueryStart(BloomTree* root, RoarQuerySet & qs, bool heuristic_match = true);

	void PrintQueryStateSet(const RoarQueryStateSet &, std::ostream &);	

protected:
	// const values and class values
	enum OPERATION { COUNT, PRIME, UPDATE };

	// jellyfish default counting values
    const uint64_t hash_size = 10000000;
    const uint32_t num_reprobes = 126;
    const uint32_t counter_len = 7;
    const bool canonical = true;

    //---- protected member variables
    double QUERY_THRESHOLD; 

    std::vector<std::string> query_sequence_list;
    RoarQuerySet query_info_set;
    RoarQueryStateSet query_state_set;
    BloomTree* bt_root;
    HashPair* hash_pair;
    int number_of_hashes;
    uint64_t bf_size;
    bool enable_debug = true;

    unsigned split_batch_query_count;

    std::unordered_map<int, std::vector<const BloomTree*> > matching_table;
	//---- protected member functions

};

#endif
