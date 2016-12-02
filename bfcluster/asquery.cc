#include "asquery.h"
#include "util.h"
#include <cassert>
#include <tclap/CmdLine.h>

float QUERY_THRESHOLD = 0.9;

// ** THIS IS NOW PARTIALLY DEPRICATED. ONLY WORKS WITH HARDCODED SIMILARITY TYPE
void assert_is_union(BloomTree* u) {
    BloomTree::protected_cache(true);
    BF* ubf = u->child(0)->bf()->union_with("union", u->child(1)->bf());
    BloomTree::protected_cache(false);

    // if the similaritity isn't 100%, we stop
    std::ostringstream oss;
    uint64_t sim = ubf->similarity(u->bf(),0);
    if (sim != ubf->num_filter_bits()) {
        std::cerr << "Filter at " << u->name() << " is not the union of its two children!" << std::endl;
        std::cerr << "Sim= " << sim << "Size= " << ubf->num_filter_bits() << std::endl;
        std::cerr << "Children:" << u->child(0)->name() << " " << u->child(1)->name() << std::endl;
        std::cerr << std::endl;
        exit(1);
    } else {
        std::cerr << "Filter at " << u->name() << " looks good." << std::endl;
        std::cerr << "Children:" << u->child(0)->name() << " " << u->child(1)->name() << std::endl;
        std::cerr << std::endl;
    }
    delete ubf;
}

void check_bt(BloomTree* root) {
    if (root == nullptr) return;

    std::cerr << "Checking node " << root->name() << std::endl;
    unsigned long nb_ones = root->bf()->count_ones();
    std::cerr << "Number of ones: " << nb_ones << std::endl;
    std::cerr << "Saturation: " << (double)nb_ones/root->bf()->num_filter_bits() << std::endl;

    if (root->child(0) && root->child(1)) {
        assert_is_union(root);
    }

    check_bt(root->child(0));
    check_bt(root->child(1));
}

void compute_intersection_across_all_leaves_descendants(BloomTree* u) {
    // strip the .roar extension
    size_t lastindex = u->name().find_last_of(".");
    std::string rawname = u->name().substr(0, lastindex);

    BloomTree::protected_cache(true);
    BF *bf1;
    bool bf1_is_leaf = !(u->child(0)->child(0) || u->child(0)->child(1));
    if (bf1_is_leaf)
        bf1 = u->child(0)->bf();
    else
        bf1 = u->child(0)->bf_inter_leaves();

    BF *bf2;
    bool bf2_is_leaf = !(u->child(1)->child(0) || u->child(1)->child(1));
    if (bf2_is_leaf)
        bf2 = u->child(1)->bf();
    else
        bf2 = u->child(1)->bf_inter_leaves();

    if ( (!bf1) || (!bf2) )
    {
        std::cerr << "!!!Problem with accessing children of" << u->name() << std::endl;
        return;
    }


    BF* ubf = bf1->intersection_with(rawname + ".inter_leaves.roar", bf2);
    ubf->save();
    BloomTree::protected_cache(false);

    unsigned long nb_ones = ubf->count_ones();
    std::cerr << "Node " << u->name() << std::endl;
    std::cerr << "Intersection of all descendants: " << nb_ones << std::endl;

    delete ubf;
}

void intersect_bt(BloomTree* root) {
    if (root == nullptr) return;

    intersect_bt(root->child(0));
    intersect_bt(root->child(1));

    if (root->child(0) && root->child(1)) {
        compute_intersection_across_all_leaves_descendants(root);
    }
}


void draw_bt_recur(BloomTree* root, std::ostream& out) {
    if (root == nullptr) return;

    std::string current = quote(test_basename(root->name(), ".bf.bv"));

    if (root->child(0)) {
        out << current << " -> " << quote(test_basename(root->child(0)->name(), ".bf.bv")) << " ; " << std::endl;
        draw_bt_recur(root->child(0), out);
    }
    if (root->child(1)) {
        out << current << " -> " << quote(test_basename(root->child(1)->name(), ".bf.bv")) << " ; " << std::endl;
        draw_bt_recur(root->child(1), out);
    }
}

void draw_bt(BloomTree* root, std::string outfile) {
    std::ofstream out(outfile.c_str());
    DIE_IF(!out, "Couldn't open output file");
    out << "digraph BloomTree {" << std::endl;
    draw_bt_recur(root, out);
    out << "}" << std::endl;
}

void compress_bt_rrr(BloomTree* root) {
    if (root == nullptr) return;

    root->bf()->compress_rrr();
    if (root->is_split()) {
        root->bf_all()->compress_rrr();
    }

    if (root->child(0)) {
        compress_bt_rrr(root->child(0));
    }
    if (root->child(1)){
        compress_bt_rrr(root->child(1));
    }
}

void uncompress_bt(BloomTree* root) {
    if (root == nullptr) return;

    root->bf()->uncompress();
    if (root->is_split()) {
        root->bf_all()->uncompress();
    }

    if (root->child(0)) {
        uncompress_bt(root->child(0));
    }
    if (root->child(1)){
        uncompress_bt(root->child(1));
    }
}

void compress_bt_roar(BloomTree* root) {
    if (root == nullptr) return;

    root->bf()->compress_roar();
    if (root->is_split()) {
        root->bf_all()->compress_roar();
    }

    if (root->child(0)) {
        compress_bt_roar(root->child(0));
    }
    if (root->child(1)){
        compress_bt_roar(root->child(1));
    }
}

inline bool file_exists (const std::string& name) {
    return ( access( name.c_str(), F_OK ) != -1 );
}

void compress_bt_roar_parallel(BloomTree* root, std::vector<std::future<void> > *F) {
    BloomTree::caching = false; // because we're going to load many bloom filters in parallel

    if (root == nullptr) return;
    bool is_root = (F == nullptr);

    if (is_root) {
        F = new std::vector<std::future<void> > ;
    }

    F->emplace_back(
            BloomTree::pool.enqueue(
                [root] () {
                    if (file_exists(root->name() + ".roar")) return;
                    std::cout << "compressing roar for " << root->name() << std::endl;
                    //std::lock_guard<std::mutex> lock(BloomTree::btLock);
                    root->bf()->compress_roar();
                    if (root->is_split()) {
                        root->bf_all()->compress_roar();
                    }
                }
            )
        );

//    if (!is_root) return;

    if (root->child(0)) {
        compress_bt_roar_parallel(root->child(0), F);
    }
    if (root->child(1)) {
        compress_bt_roar_parallel(root->child(1), F);
    }

    if (is_root) {
        std::cerr << "Waiting for subtrees to finish." << std::endl;
        // join on all the threads
        for (auto & f : *F) {
            f.wait();
        }
    }
}

/*============================================*/

bool query_passes(BloomTree* root, const std::set<jellyfish::mer_dna> & q) {
    assert(q.size() > 0);
    auto bf = root->bf();
    unsigned pass_thresh = ceil(QUERY_THRESHOLD * q.size());
    unsigned fail_thresh = q.size() - pass_thresh;
    unsigned good = 0, bad = 0;
    //DEBUG:
    //std::cerr << "query_passes(kmer set)"
    //          << " q.size=" << q.size()
    //          << " pass_thresh=" << pass_thresh
    //          << " fail_thresh=" << fail_thresh
    //          << std::endl;
    for (const auto & m : q) {
        if (bf->contains(m)) {
            if (++good >= pass_thresh) return true;
        } else {
            if (++bad > fail_thresh) return false;
        }
        //DEBUG: std::cerr << "checking: " << m.to_str()
        //                 << " " << good << ":" << bad << std::endl;
    }
    return (good >= pass_thresh);  // this statement is never reached
}

// return true if the filter at this node contains > QUERY_THRESHOLD kmers
bool query_passes(BloomTree* root, QueryInfo* q) {
    //assert(q.size() > 0);
    auto bf = root->bf();

    // UNweighted case
    if (q->weight.empty()) {
        unsigned pass_thresh = ceil(QUERY_THRESHOLD * q->query_kmers.size());
        unsigned fail_thresh = q->query_kmers.size() - pass_thresh;
        unsigned good = 0, bad = 0;
        //DEBUG:
        //std::cerr << "query_passes(QueryInfo)"
        //          << " q->query_kmers.size=" << q->query_kmers.size()
        //          << " pass_thresh=" << pass_thresh
        //          << " fail_thresh=" << fail_thresh
        //          << std::endl;
        for (const auto & m : q->query_kmers) {
            if (bf->contains(m)) {
                if (++good >= pass_thresh) {
                    //DEBUG
                    //std::cout << "good " << good << "/" << (good+bad) << std::endl;
                    return true;
                }
            } else {
                if (++bad > fail_thresh) {
                    //DEBUG
                    //std::cout << "bad " << bad << "/" << (good+bad) << std::endl;
                    return false;
                }
            }
            //DEBUG: std::cerr << "checking: " << m << " " << m.to_str()
            //                 << " " << good << ":" << bad << std::endl;
        }
        return (good >= pass_thresh);  // this statement is never reached

    // WEIGHTED case
    } else {
        // nota bene: threshold ought to be computed as float, and ought to be
        // based on the sum of the weights; to implement a fail_thresh we'd
        // need to know the sum of the weights
        float pass_thresh = QUERY_THRESHOLD * q->query_kmers.size();
        float good = 0;
        unsigned n = 0;
        for (const auto & m : q->query_kmers) {
            if (n >= q->weight.size()) {
                std::cerr << "Number of weights (" << q->weight.size() <<") less than query kmers (" << q->query_kmers.size() << ")."  << std::endl;
                exit(3);
            }
            if (bf->contains(m)) {
                good += q->weight[n];
                if (good >= pass_thresh) return true;
            }
            n++;
            //DEBUG: std::cerr << "checking: " << m << " " << m.to_str() << " " << good << std::endl;
        }
        return (good >= pass_thresh);  // this statement IS reachable
    }
}

bool original_query_passes(BloomTree* root, QueryInfo* q) {
    //--- this should match the version in the Solomon paper ---
    float weight = 1.0;
    auto bf = root->bf();
    float c = 0;
    unsigned n = 0;
    bool weighted = 0;
    if (q->weight.empty()) {
        weighted = 0;
    } else {
        weighted = 1;
    }

    for (const auto & m : q->query_kmers) {
        //DEBUG: std::cout << "checking: " << m.to_str();
        if (weighted) {
            if (q->weight.size() > n) {
                weight=q->weight[n];
            } else {
                std::cerr << "Number of weights (" << q->weight.size() <<") less than query kmers (" << q->query_kmers.size() << ")."  << std::endl;
                exit(3);
            }
        }
        if (bf->contains(m)) c += weight;
        n++;
        //DEBUG: std::cout << c << std::endl;
    }
    //std::cerr << "comparing " << c << " vs " << (QUERY_THRESHOLD * q->query_kmers.size()) << std::endl;
    return (c >= QUERY_THRESHOLD * q->query_kmers.size());
}

// 'try' a query at a node and return the return the status of the query
QueryState* query_by_reduction_status(BloomTree* root, QueryState* q_state) {
    // assert (q_state->q_info->weight.empty());
    auto bf = root->bf();
    unsigned fail_thresh = q_state->fail_thresh;
    unsigned num_passed = 0;
    QueryState* still_active = new QueryState(q_state);

    for (const auto & pos : q_state->kmer_positions) {
        if (bf->contains(pos)) {
            num_passed++;
            // $$$ where does the space for this set get de-allocated?
            still_active->kmer_positions.emplace_back(pos);
            // note that we can't shortcut the exit, like this:
            //   if (num_passed >= pass_thresh) break;
            // because all the positives remaining in the list must be copied
        } else {
            still_active->num_failed++;
            // note that the following can never be true if num_passed has
            // reached pass_thresh
            if (still_active->num_failed > fail_thresh) break;
        }
        //DEBUG: std::cerr << "checking: " << pos
        //                 << " " << num_passed << ":" << still_active->num_failed << std::endl;
    }
    return still_active;
}

// 'try' a query at a node and return the return the status of the query
QueryState* split_query_by_reduction_status(BloomTree* root, QueryState* q_state) {
    // assert (q_state->q_info->weight.empty());
    auto bf_all = root->bf_all();
    auto bf_some = root->bf();
    unsigned pass_thresh = q_state->pass_thresh;
    unsigned fail_thresh = q_state->fail_thresh;
    QueryState* still_active = new QueryState(q_state);

    for (const auto & pos : q_state->kmer_positions) {
        if (bf_all->contains(pos)) {
            still_active->num_passed++;
            if (still_active->num_passed >= pass_thresh) break;
        } else if (bf_some->contains(pos)) {
            // $$$ where does the space for this set get de-allocated?
            still_active->kmer_positions.emplace_back(pos);
        } else {
            still_active->num_failed++;
            if (still_active->num_failed > fail_thresh) break;
        }
        //DEBUG: std::cerr << "checking: " << pos
        //                 << " " << still_active->num_passed << ":" << still_active->num_failed << std::endl;
    }
    return still_active;
}

// recursively walk down the tree, proceeding to children only
// if their parent passes the query threshold;
int query(
    BloomTree* root,
    const std::set<jellyfish::mer_dna> & q,
    std::vector<BloomTree*> & out
) {
    int nodes_visited = 1;
    root->increment_usage();
    if (query_passes(root, q)) {
        //DEBUG: std::cout << "passed at " << root->name() << std::endl;
        int children = 0;
        if (root->child(0)) {
            nodes_visited += query(root->child(0), q, out);
            children++;
        }
        if (root->child(1)) {
            nodes_visited += query(root->child(1), q, out);
            children++;
        }
        if (children == 0) {
            out.push_back(root);
        }
    }
    return nodes_visited;
}

// same as query() but the string is first converted into a set of kmers.
int query_string(
    BloomTree* root,
    const std::string & q,
    std::vector<BloomTree*> & out
) {
    return query(root, kmers_in_string(q), out);
}

// read 1 query per line, execute it, and print to the output stream o the
// results in the format:
//      *QUERY number_results
//      # number_nodes_visited
//      BF names
//
void query_from_file(
    BloomTree* root,
    const std::string & fn,
    std::ostream & o
) {
    std::vector<BloomTree*> out;
    std::string line;

    std::ifstream in(fn);
    DIE_IF(!in.good(), "Couldn't open query file.");
    while (getline(in, line)) {
        line = Trim(line);
        if (line.size() < jellyfish::mer_dna::k()) continue;

        o << "*" << line;
        int nodes_visited = query_string(root, line, out);
        o << " " << out.size() << std::endl;
        o << "# " << nodes_visited << " nodes visited" << std::endl;

        for (const auto& n : out) {
            o << n->name() << std::endl;
        }

        out.clear();
    }
}

/******
 * Batch querying
 ******/

void print_query_results(const QuerySet & qs, std::ostream & out) {
    for (auto& q : qs) {
        out << "*" << q->query << " " << q->matching.size() << std::endl;
        out << "# " << q->nodes_visited << " nodes visited" << std::endl;
        for (const auto& n : q->matching) {
            out << n->name() << std::endl;
        }
    }
}


void dump_positions(const std::vector<uint64_t> positions) {
    std::cerr << "  ";
    for (const auto & pos : positions) {
        std::cerr << " " << pos;
    }
    std::cerr << std::endl;
}


void query_batch(BloomTree* root, QuerySet & qs) {
    // how many children do we have?
    bool has_children = root->child(0) || root->child(1);

    // construct the set of queries that pass this node
    QuerySet pass;
    unsigned n = 0;
    for (auto & q : qs) {
        q->nodes_visited++;
        if (query_passes(root, q)) {
            if (has_children) {
                pass.emplace_back(q);
            } else {
                q->matching.emplace_back(root);
                n++;
            }
        }
    }

    // $(node name) $(internal / leaf) $(number of matches)
    if (has_children) {
        std::cout << "[search] " << root->name() << " internal " << pass.size() << "/" << qs.size() << std::endl;
    } else {
        std::cout << "[search] " << root->name() << " leaf "     << n           << "/" << qs.size() << std::endl;
    }

    if (pass.size() > 0) {
        // if present, recurse into left child
        if (root->child(0)) {
            query_batch(root->child(0), pass);
        }

        // if present, recurse into right child
        if (root->child(1)) {
            query_batch(root->child(1), pass);
        }
    }
}

void original_query_batch(BloomTree* root, QuerySet & qs) {
    //--- this should match the version in the Solomon paper ---
    // how many children do we have?
    bool has_children = root->child(0) || root->child(1);

    // construct the set of queries that pass this node
    QuerySet pass;
    unsigned n = 0;
    for (auto & q : qs) {
        q->nodes_visited++;
        if (original_query_passes(root, q)) {
            if (has_children) {
                pass.emplace_back(q);
            } else {
                q->matching.emplace_back(root);
                n++;
            }
        }
    }

    // $(node name) $(internal / leaf) $(number of matches)
    if (has_children) { //Changing format
        std::cout << root->name() << " internal " << pass.size() << std::endl;
    } else {
        std::cout << root->name() << " leaf " << n << std::endl;
    }

    if (pass.size() > 0) {
        // if present, recurse into left child
        if (root->child(0)) {
            original_query_batch(root->child(0), pass);
        }

        // if present, recurse into right child
        if (root->child(1)) {
            original_query_batch(root->child(1), pass);
        }
    }
}

/*============================================*/

static unsigned batch_query_count;  // this only counts nodes visited, maybe less than all nodes

void batch_query_by_reduction(BloomTree* root, QueryStateSet active_queries) {
    bool has_children = root->child(0) || root->child(1);

    batch_query_count++;

    QueryStateSet undetermined;
    unsigned n = 0;
    for (auto & q_state : active_queries) {
        q_state->q->nodes_visited++;
        QueryState* query_result = query_by_reduction_status(root, q_state);
        if (query_result->num_failed > query_result->fail_thresh) {
            // this query is not considered further
            delete query_result;
        } else if (has_children) {
            undetermined.emplace_back(query_result);
        } else {
            delete query_result;
            q_state->q->matching.emplace_back(root);
            n++;
        }
    }

    // $(node name) $(internal / leaf) $(number of matches)
    if (has_children) {
        std::cout << "[search " << batch_query_count << "] " << root->name() << " internal " << undetermined.size() << "/" << active_queries.size() << std::endl;
    } else {
        std::cout << "[search " << batch_query_count << "] " << root->name() << " leaf "     << n                   << "/" << active_queries.size() << std::endl;
    }

    if (undetermined.size() > 0) {
        // if present, recurse into left child
        if (root->child(0)) {
            batch_query_by_reduction(root->child(0), undetermined);
        }

        // if present, recurse into right child
        if (root->child(1)) {
            batch_query_by_reduction(root->child(1), undetermined);
        }

        for (auto & q_state : undetermined) {
            delete q_state;
        }
    }
}

void batch_query_by_reduction_start(BloomTree* root, QuerySet & qs) {
    // construct a state for each query, reducing kmers to bf-positions
    QueryStateSet active_queries;
    for (auto & q : qs) {
        QueryState* q_state = new QueryState(root->bf(),q);
        q_state->pass_thresh = ceil(QUERY_THRESHOLD * q->query_kmers.size());
        q_state->fail_thresh = q->query_kmers.size() - q_state->pass_thresh;
        active_queries.emplace_back(q_state);
        //DEBUG:
        //std::cerr << "query to bf positions" << std::endl;
        //std::cerr << q->query << std::endl;
        //dump_positions(q_state->kmer_positions);
    }

    // perform the queries
    batch_query_count = 0;
    batch_query_by_reduction(root, active_queries);

    // delete the query states (but not the query results)
    for (auto & q_state : active_queries) {
        delete q_state;
    }
}

/*============================================*/

static unsigned split_batch_query_count;  // this only counts nodes visited, maybe less than all nodes

void split_batch_query_all_leaves_match(BloomTree* root, QueryState* q_state) {
    bool has_children = root->child(0) || root->child(1);

    if (has_children) {
        if (root->child(0)) split_batch_query_all_leaves_match(root->child(0), q_state);
        if (root->child(1)) split_batch_query_all_leaves_match(root->child(1), q_state);
    } else {
        q_state->q->matching.emplace_back(root);
    }
}

void split_batch_query_by_reduction(BloomTree* root, QueryStateSet active_queries) {
    bool has_children = root->child(0) || root->child(1);

    split_batch_query_count++;

    QueryStateSet undetermined;
    unsigned n = 0;
    for (auto & q_state : active_queries) {
        q_state->q->nodes_visited++;
        QueryState* query_result = split_query_by_reduction_status(root, q_state);
        if (query_result->num_passed >= query_result->pass_thresh) {
            if (has_children) {
                if (root->child(0)) split_batch_query_all_leaves_match(root->child(0), q_state);
                if (root->child(1)) split_batch_query_all_leaves_match(root->child(1), q_state);
            } else {
                q_state->q->matching.emplace_back(root);
                n++;
            }
            delete query_result;
        } else if (query_result->num_failed > query_result->fail_thresh) {
            // this query is not considered further
            delete query_result;
        } else if (has_children) {
            undetermined.emplace_back(query_result);
        } else {
            DIE("Internal Error, split tree apparently built incorrectly");
        }
    }

    // $(node name) $(internal / leaf) $(number of matches)
    if (has_children) {
        std::cout << "[search " << split_batch_query_count << "] " << root->name() << " internal " << undetermined.size() << "/" << active_queries.size() << std::endl;
    } else {
        std::cout << "[search " << split_batch_query_count << "] " << root->name() << " leaf "     << n                   << "/" << active_queries.size() << std::endl;
    }

    if (undetermined.size() > 0) {
        // if present, recurse into left child
        if (root->child(0)) {
            split_batch_query_by_reduction(root->child(0), undetermined);
        }

        // if present, recurse into right child
        if (root->child(1)) {
            split_batch_query_by_reduction(root->child(1), undetermined);
        }

        for (auto & q_state : undetermined) {
            delete q_state;
        }
    }
}

void split_batch_query_by_reduction_start(BloomTree* root, QuerySet & qs) {
    // construct a state for each query, reducing kmers to bf-positions
    QueryStateSet active_queries;
    for (auto & q : qs) {
        QueryState* q_state = new QueryState(root->bf(),q);
        q_state->pass_thresh = ceil(QUERY_THRESHOLD * q->query_kmers.size());
        q_state->fail_thresh = q->query_kmers.size() - q_state->pass_thresh;
        active_queries.emplace_back(q_state);
        //DEBUG:
        //std::cerr << "query to bf positions" << std::endl;
        //std::cerr << q->query << std::endl;
        //dump_positions(q_state->kmer_positions);
    }

    // perform the queries
    split_batch_query_count = 0;
    split_batch_query_by_reduction(root, active_queries);

    // delete the query states (but not the query results)
    for (auto & q_state : active_queries) {
        delete q_state;
    }
}

/*============================================*/

void query_leaves (BloomTree* root, QuerySet & qs) {
    // how many children do we have?
    bool has_children = root->child(0) || root->child(1);

    // construct the set of queries that pass this node
    // But only for leaf nodes
    unsigned n=0;
    if (!has_children) {
        for (auto & q : qs) {
            q->nodes_visited++;
            if (query_passes(root, q)) {
                q->matching.emplace_back(root);
                n++;
            }
        }
    }

    // $(node name) $(internal / leaf) $(number of matches)
    if (!has_children) {
        std::cout << "[search] " << root->name() << " leaf " << n << "/" << qs.size() << std::endl;
    }

    // if present, recurse into left child
    if (root->child(0)) {
        query_leaves(root->child(0), qs);
    }

    // if present, recurse into right child
    if (root->child(1)) {
        query_leaves(root->child(1), qs);
    }

}

void batch_query_from_file(
    BloomTree* root,
    const std::string & fn,
    std::ostream & o
) {
    // read in the query lines from the file.
    std::string line;
    QuerySet qs;
    std::ifstream in(fn);
    DIE_IF(!in.good(), "Couldn't open query file.");
    std::size_t n = 0;
    while (getline(in, line)) {
        line = Trim(line);
        // RSH: reject fasta header lines (a common error case)
        if (starts_with(line,">")) continue;
        if (line.size() < jellyfish::mer_dna::k()) continue;
        QueryInfo* q = new QueryInfo(line);
        qs.emplace_back(q);
        n++;
    }
    in.close();
    std::cerr << "Read " << n << " queries." << std::endl;

    // batch process the queries
    if (root->is_split())
        split_batch_query_by_reduction_start(root, qs);
    else
        query_batch(root, qs);
    print_query_results(qs, o);

    // free the query info objects
    for (auto & p : qs) {
        delete p;
    }
}



void original_batch_query_from_file(
    BloomTree* root,
    const std::string & fn,
    std::ostream & o
) {
    //--- this should match the version in the Solomon paper ---
    DIE_IF(root->is_split(), "query-original doesn't support split-node trees");

    // read in the query lines from the file.
    std::string line;
    QuerySet qs;
    std::ifstream in(fn);
    DIE_IF(!in.good(), "Couldn't open query file.");
    std::size_t n = 0;
    while (getline(in, line)) {
        line = Trim(line);
        if (line.size() < jellyfish::mer_dna::k()) continue;
        qs.emplace_back(new QueryInfo(line));
        n++;
    }
    in.close();
    std::cerr << "Read " << n << " queries." << std::endl;

    // batch process the queries
    original_query_batch(root, qs);
    print_query_results(qs, o);

    // free the query info objects
    for (auto & p : qs) {
        delete p;
    }
}

void batch_query_by_reduction_from_file(
    BloomTree* root,
    const std::string & fn,
    std::ostream & o
) {
    // read in the query lines from the file.
    std::string line;
    QuerySet qs;
    std::ifstream in(fn);
    DIE_IF(!in.good(), "Couldn't open query file.");
    std::size_t n = 0;
    while (getline(in, line)) {
        line = Trim(line);
        // RSH: reject fasta header lines (a common error case)
        if (starts_with(line,">")) continue;
        if (line.size() < jellyfish::mer_dna::k()) continue;
        QueryInfo* q = new QueryInfo(line);
        qs.emplace_back(q);
        n++;
    }
    in.close();
    std::cerr << "Read " << n << " queries." << std::endl;

    // batch process the queries
    if (root->is_split())
        split_batch_query_by_reduction_start(root, qs);
    else
        batch_query_by_reduction_start(root, qs);
    print_query_results(qs, o);

    // free the query info objects
    for (auto & p : qs) {
        delete p;
    }
}

void batch_weightedquery_from_file(
    BloomTree* root,
    const std::string & fn,
    const std::string & wf,
    std::ostream & o
) {
    // read in the query lines from the file.
    std::string line;
    std::string wfline;
    QuerySet qs;
    std::ifstream in(fn);
    std::ifstream wfin(wf);
    DIE_IF(root->is_split(), "Weighted queries for this tree file aren't currently support supported.");
    DIE_IF(!in.good(), "Couldn't open query file.");
    DIE_IF(!wfin.good(), "Couldn't open weight file.");
    std::size_t n = 0;
    while (getline(in, line)) {
        getline(wfin, wfline);
        line = Trim(line);
        wfline = Trim(wfline);
        // RSH: rejecting fasta header lines would be more difficult here
        if (line.size() < jellyfish::mer_dna::k()) continue;
        qs.emplace_back(new QueryInfo(line, wfline));
        n++;
    }
    in.close();
    std::cerr << "Read " << n << " queries." << std::endl;

    // batch process the queries
    query_batch(root, qs);
    print_query_results(qs, o);

    // free the query info objects
    for (auto & p : qs) {
        delete p;
    }
}


void leaf_query_from_file(
    BloomTree* root,
    const std::string & fn,
    std::ostream & o
) {
    std::string line;
    QuerySet qs;
    DIE_IF(root->is_split(), "Tree file doesn't support leaf queries.");
    std::ifstream in(fn);
    DIE_IF(!in.good(), "Couldn't open query file.");
    std::size_t n=0;
    while (getline(in, line)) {
        line = Trim(line);
        // RSH: reject fasta header lines (a common error case)
        if (starts_with(line,">")) continue;
        if (line.size() < jellyfish::mer_dna::k()) continue;
        qs.emplace_back(new QueryInfo(line));
        n++;
    }
    in.close();
    std::cerr << "Read " << n << " queries." << std::endl;

    // batch process the queries on ONLY the leaves
    query_leaves(root, qs);
    print_query_results(qs, o);

    for (auto & p : qs) {
        delete p;
    }
}

// add by chen sun, 10/17/2016

void BatchQueryByKmerSet(
    BloomTree* root,
    const std::string & fn,
    std::ostream & o)
{
    // read in the query lines from the file.
    std::string line;
    QuerySet qs;
    std::ifstream in(fn);
    QueryInfo* q = new QueryInfo();
    DIE_IF(!in.good(), "Couldn't open query file.");
    std::size_t n = 0;
    while (getline(in, line)) {
        line = Trim(line);
        // RSH: reject fasta header lines (a common error case)
        if (starts_with(line,">")) continue;
        if (line.size() < jellyfish::mer_dna::k()) continue;
        q->InsertKmer(line);
        n++;
    }
    in.close();
    qs.emplace_back(q);
    std::cerr << "Read " << n << " queries." << std::endl;

    // batch process the queries
    if (root->is_split())
        split_batch_query_by_reduction_start(root, qs);
    else
        query_batch(root, qs);
    print_query_results(qs, o);

    // free the query info objects
    for (auto & p : qs) {
        delete p;
    }  
}

typedef struct Args {
    std::string bloomtree_filename;
    std::string query_filename;
    std::string hash_filename;
    std::string output_filename;
    double matching_threshold;
    bool heuristic_match;
    bool multiple_match;
    bool enable_debug;
}Args;

bool TclapParser(Args & args, int argc, char** argv){
    std::string version = "0.9";
    try {
        std::string desc = "squery \n";
        TCLAP::CmdLine cmd(desc, ' ', version);
        TCLAP::ValueArg<std::string> arg_bloomtree_filename("b", "bloomtree", "bloom tree topology file", true, "", "file");
        TCLAP::ValueArg<std::string> arg_query_filename("q", "query", "query FASTA file", true, "", "file");
        //TCLAP::ValueArg<std::string> arg_output_dir("d", "output_dir", "output directory, default is current working directory", false, ".", "string");
        TCLAP::ValueArg<std::string> arg_output_filename("o", "output", "output file", true, "", "file");

        TCLAP::ValueArg<double> arg_matching_threshold("t", "threshold", "matching threshold (default=0.9)", false, 0.9, "double");
        
        TCLAP::SwitchArg arg_enable_debug("d", "debug", "output debug information", cmd, false);

        TCLAP::SwitchArg arg_heuristic_match("E", "exact_match", "calculate exact number of matching kmers (defaultly not)", cmd, true);

        TCLAP::SwitchArg arg_multiple_match("m", "multiple_match", "consider entries in FASTA file separately (defaultly consider all kmers as one query)", cmd, false);

        cmd.add(arg_matching_threshold);
        cmd.add(arg_output_filename);
        //cmd.add(arg_output_dir);
        cmd.add(arg_query_filename);
        //cmd.add(arg_hash_filename);
        cmd.add(arg_bloomtree_filename);
        cmd.parse(argc, argv);

        args.bloomtree_filename = arg_bloomtree_filename.getValue();

        args.query_filename = arg_query_filename.getValue();
        args.output_filename = arg_output_filename.getValue();
        args.matching_threshold = arg_matching_threshold.getValue();
        args.heuristic_match = arg_heuristic_match.getValue();
        args.enable_debug = arg_enable_debug.getValue();
        args.multiple_match = arg_multiple_match.getValue();
    }
    catch (TCLAP::ArgException &e)
    {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << "\n";
        abort();
    }
    return true;
}

//void file_exists (const std::string& filename) {
//    struct stat buffer;   
//    if (stat (filename.c_str(), &buffer) != 0) 
//        printf( "cannot access %s\n", filename.c_str() );
//}

void dir_exists(const std::string& pathname){
    struct stat info;
    if( stat( pathname.c_str(), &info ) != 0 )
        printf( "cannot access %s\n", pathname.c_str() );
    else if(! (info.st_mode & S_IFDIR)) 
        printf( "%s is no directory\n", pathname.c_str() );
}

int main(int argc, char* argv[]){
    
    BloomTree::caching = true;

    Args args;
    TclapParser(args, argc, argv);

    file_exists(args.bloomtree_filename);
    file_exists(args.query_filename);

    QUERY_THRESHOLD = args.matching_threshold;
    std::cout << "searching threshold: " << QUERY_THRESHOLD << std::endl;

    BloomTree* root = read_bloom_tree(args.bloomtree_filename);

    std::ofstream out(args.output_filename);

    BatchQueryByKmerSet( root,
                    args.query_filename,
                    out);
	
	//batch_query_from_file(root, args.query_filename, out);

    out.close();

    return 0;
}
