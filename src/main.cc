#include "Query.h"
#include "Build.h"
#include "BloomTree.h"
#include "BF.h"
#include "util.h"
#include "Count.h"

#include <string>
#include <cstdlib>
#include <getopt.h>

#include <jellyfish/file_header.hpp>

/* TODO:
 */

// various commandline filenames
std::string command;
std::string bloom_tree_file;
std::string query_file;
std::string leaves_file;
std::string fasta_file;
std::string out_file;
std::string jfbloom_file;
std::string bvfile1, bvfile2, bvfile3;
std::string nt_string;
int sim_type=0;
std::string bloom_storage;
int leaf_only;
std::string weighted="";
unsigned cutoff_count=3;
int bv_dump_spacing = 10;

std::string hashes_file;
unsigned nb_hashes = 1;
uint64_t bf_size;

unsigned num_threads = 16;
//unsigned parallel_level = 3; // no parallelism by default

const char * OPTIONS = "t:p:f:l:c:w:s:z:";

static struct option LONG_OPTIONS[] = {
    {"max-filters", required_argument, 0, 'f'},
    {"threads", required_argument, 0, 'p'},
    {"query-threshold", required_argument, 0, 't'},
    {"k", required_argument, 0, 'k'},
    {"sim-type", required_argument, 0, 's'},
    {"leaf-only", required_argument,0,'l'},
    {"cutoff", required_argument,0,'c'},
    {"weighted", required_argument,0,'w'},
    {"spacing", required_argument,0,'z'},
    {0,0,0,0}
};

void print_usage() {
    std::cerr
        << "Usage: bt [query|convert|build] ...\n"
        << "   Construction:\n"
        << "    \"hashes\" [-k 20] hashfile\n"
        << "    \"count\" [--cutoff 3] [--threads 16] hashfile bf_size fasta_in filter_out.bf.bv\n"
        << "    \"build\" [--sim-type 0] hashfile filterlistfile outfile\n"
        << "    \"compress-rrr\" bloomtreefile outfile\n"
        << "    \"compress-rrr-single\" hashfile filter_in.bf.bv\n"
        << "    \"compress-roar\" bloomtreefile outfile\n"
        << "    \"compress-roar-single\" hashfile filter_in.bf.bv\n"
        << "    \"uncompress\" bloomtreefile\n"
        << "    \"split\" bloomtreefile outfile\n"
        << "    \"rebuild\" bloomtreefile\n"
        << "    \"intersection-leaves\" bloomtreefile\n"
        << "    \"append\" hashfile filter1_in.bf.bv filter2_in.bf.bv filter_out.bf.bv\n"
        << "    \"compress-rrr-double\" hashfile filter1_in.bf.bv filter2_in.bf.bv filter_out.bf.bv.rrr\n"
        << "    \"compress-wraprrr-double\" hashfile filter1_in.bf.bv filter2_in.bf.bv\n"
        << "   Misc:\n"
        << "    \"check\" bloomtreefile\n"
        << "    \"draw\" bloomtreefile out.dot\n"
        << "    \"convert\" jfbloomfilter outfile\n"
        << "    \"sim\" [--sim-type 0] bloombase bvfile1 bvfile2\n"
        << "    \"bv-stats\" bloomtreefile\n"
        << "    \"bv-dump\" [--spacing 10] bloomtreefile\n"
        << "    \"iotest\" bloomtreefile\n"
        << "   Query:\n"
        << "    \"query\" [--max-filters 1] [--query-threshold " << QUERY_THRESHOLD << "] [--leaf-only 0] [--weighted weightfile] bloomtreefile queryfile outfile\n"
        << "    \"query-redux\" [--max-filters 1] [--query-threshold " << QUERY_THRESHOLD << "] [--leaf-only 0] [--weighted weightfile] bloomtreefile queryfile outfile\n"
        << "    \"query-original\" [--max-filters 1] [--query-threshold " << QUERY_THRESHOLD << "] [--leaf-only 0] [--weighted weightfile] bloomtreefile queryfile outfile\n"
        << std::endl;
    exit(3);
}

// construct a new set of hashes for the current k
// also record k and the number of hash applications
void construct_hashes(std::string & hashesfile, int nh) {
    HashPair hp;
    jellyfish::file_header fh;
    fh.matrix(hp.m1, 1);
    fh.matrix(hp.m2, 2);
    fh.key_len(jellyfish::mer_dna::k() * 2);
    fh.nb_hashes(nh);
    std::ofstream hashesout(hashesfile.c_str());
    fh.write(hashesout);
    hashesout.close();
}


int process_options(int argc, char* argv[]) {
    int k = 20;
    //int sim_type = 0;
    int a;
    while ((a=getopt_long(argc, argv, OPTIONS, LONG_OPTIONS, 0)) != -1) {
        switch(a) {
            case 't':
                QUERY_THRESHOLD = atof(optarg);
                break;
            case 'p':
                num_threads = unsigned(atoi(optarg));
                //parallel_level = unsigned(atoi(optarg));
                break;
            case 'f':
                BF_INMEM_LIMIT = unsigned(atoi(optarg));
                break;
            case 'l':
                leaf_only = atoi(optarg);
            case 'k':
                k = atoi(optarg);
                break;
            case 's':
                sim_type = atoi(optarg);
                break;
            case 'c':
                cutoff_count = unsigned(atoi(optarg));
                break;
            case 'w':
                weighted = optarg;
                break;
            case 'z':
                bv_dump_spacing = atoi(optarg);
                if (bv_dump_spacing < 1) bv_dump_spacing = 1;
                break;
            default:
                std::cerr << "Unknown option." << std::endl;
                print_usage();
        }
    }

    jellyfish::mer_dna::k(k);
    std::cerr << "Kmer size = " << jellyfish::mer_dna::k()
              << " (but we will use the k from a hashfile if we read one)" << std::endl;

    if (optind >= argc) print_usage();
    command = argv[optind];
    if ((command == "query") || (command == "query-redux") || (command == "query-original")) {
        if (optind >= argc-3) print_usage();
        bloom_tree_file = argv[optind+1];
        query_file = argv[optind+2];
        out_file = argv[optind+3];

    } else if (command == "check" || command == "intersection-leaves" || command == "iotest") {
        if (optind >= argc-1) print_usage();
        bloom_tree_file = argv[optind+1];

    } else if (command == "draw") {
        if (optind >= argc-2) print_usage();
        bloom_tree_file = argv[optind+1];
        out_file = argv[optind+2];

    } else if (command == "sim") {
        if (optind >= argc-3) print_usage();
        jfbloom_file = argv[optind+1];
        bvfile1 = argv[optind+2];
        bvfile2 = argv[optind+3];
     // sim_type = argv[optind+4];

    } else if (command == "convert") {
        if (optind >= argc-2) print_usage();
        jfbloom_file = argv[optind+1];
        out_file = argv[optind+2];

    } else if (command == "build") {
        if (optind >= argc-3) print_usage();
        hashes_file = argv[optind+1];
        leaves_file = argv[optind+2];
        out_file = argv[optind+3];
        //bloom_storage = argv[optind+4];
        // sim_type = argv[optind+4];

    } else if (command == "hashes") {
        // we allow the unadvertised setting of nb_hashes, for backward
        // compatibility with existing scripts; but were prefer nb_hashes=1,
        // since any other value will break the all-some split paradigm
        if (optind >= argc-1) print_usage();
        if (optind == argc-2) {
            hashes_file = argv[optind+1];
            nb_hashes = 1;
        } else if (optind == argc-3) {
            hashes_file = argv[optind+1];
            nb_hashes = atoi(argv[optind+2]);
        }
    } else if (command == "count") {
        if (optind >= argc-4) print_usage();
        hashes_file = argv[optind+1];
        bf_size = atol(argv[optind+2]);
        fasta_file = argv[optind+3];
        out_file = argv[optind+4];

    } else if (command == "compress" || command == "compress-rrr" || command == "compress-roar") {
        if (optind >= argc-2) print_usage();
        bloom_tree_file = argv[optind+1];
        out_file = argv[optind+2];

    } else if (command == "uncompress") {
        if (optind >= argc-1) print_usage();
        bloom_tree_file = argv[optind+1];

    } else if (command == "compress-rrr-single" || command == "compress-roar-single") {
        if (optind >= argc-2) print_usage();
        hashes_file = argv[optind+1];
        bvfile1 = argv[optind+2];

    } else if (command == "split") {
        if (optind >= argc-2) print_usage();
        bloom_tree_file = argv[optind+1];
        out_file = argv[optind+2];

    } else if (command == "rebuild") {
        if (optind >= argc-1) print_usage();
        bloom_tree_file = argv[optind+1];

    } else if (command == "bv-stats") {
        if (optind >= argc-1) print_usage();
        bloom_tree_file = argv[optind+1];

    } else if (command == "bv-dump") {
        if (optind >= argc-1) print_usage();
        bloom_tree_file = argv[optind+1];

    } else if (command == "append") {
        if (optind >= argc-4) print_usage();
        hashes_file = argv[optind+1];
        bvfile1 = argv[optind+2];
        bvfile2 = argv[optind+3];
        bvfile3 = argv[optind+4];

    } else if (command == "compress-rrr-double") {
        if (optind >= argc-4) print_usage();
        hashes_file = argv[optind+1];
        bvfile1 = argv[optind+2];
        bvfile2 = argv[optind+3];
        bvfile3 = argv[optind+4];

    } else if (command == "compress-wraprrr-double") {
        if (optind >= argc-3) print_usage();
        hashes_file = argv[optind+1];
        bvfile1 = argv[optind+2];
        bvfile2 = argv[optind+3];

    }

    return optind;
}



int main(int argc, char* argv[]) {
    std::cerr << "Starting Bloom Tree" << std::endl;
    BloomTree::caching = true; // enable caching of bloom filters by default (default SBT behavior)

    process_options(argc, argv);

    if (command == "query") {
        std::cerr << "Loading bloom tree topology: " << bloom_tree_file
            << std::endl;
        BloomTree* root = read_bloom_tree(bloom_tree_file);

        std::cerr << "In memory limit = " << BF_INMEM_LIMIT << std::endl;

        std::cerr << "Querying..." << std::endl;
        std::ofstream out(out_file);
        if (leaf_only == 1) {
            leaf_query_from_file(root, query_file, out);
        } else if (weighted!="") {
            std::cerr << "Weighted query \n";
            batch_weightedquery_from_file(root, query_file, weighted, out);
        } else {
            batch_query_from_file(root, query_file, out);
        }

    } else if (command == "query-redux") {
        std::cerr << "Loading bloom tree topology: " << bloom_tree_file
            << std::endl;
        BloomTree* root = read_bloom_tree(bloom_tree_file);

        std::cerr << "In memory limit = " << BF_INMEM_LIMIT << std::endl;

        std::cerr << "Querying..." << std::endl;
        std::ofstream out(out_file);
        if (leaf_only == 1) {
            DIE("query-redux doesn't support leaf queries");
        } else if (weighted!="") {
            DIE("query-redux doesn't support weighted queries");
        } else {
            batch_query_by_reduction_from_file(root, query_file, out);
        }

    } else if (command == "query-original") {
        std::cerr << "Loading bloom tree topology: " << bloom_tree_file
            << std::endl;
        BloomTree* root = read_bloom_tree(bloom_tree_file);

        std::cerr << "In memory limit = " << BF_INMEM_LIMIT << std::endl;

        std::cerr << "Querying..." << std::endl;
        std::ofstream out(out_file);
        if (leaf_only == 1) {
            DIE("query-original doesn't support leaf queries");
        } else if (weighted!="") {
            DIE("query-original doesn't support weighted queries");
        } else {
            original_batch_query_from_file(root, query_file, out);
        }

    } else if (command == "draw") {
        std::cerr << "Drawing tree in " << bloom_tree_file << " to " << out_file << std::endl;
        BloomTree* root = read_bloom_tree(bloom_tree_file, false);
        draw_bt(root, out_file);

    } else if (command == "check") {
        BloomTree* root = read_bloom_tree(bloom_tree_file);
        std::cerr << "Checking tree" << std::endl;
        check_bt(root);

    } else if (command == "iotest") {
        BloomTree* root = read_bloom_tree(bloom_tree_file);
        std::cerr << "Checking filesystem IO stats" << std::endl;
        iotest(root);

    } else if (command == "intersection-leaves") {
        BloomTree* root = read_bloom_tree(bloom_tree_file);
        std::cerr << "Computing intersection of all descendant leaves" << std::endl;
        intersect_bt(root);

    } else if (command == "sim") {
        // read hash functions
        int num_hash;
        HashPair* hashes = get_hash_function(jfbloom_file, num_hash);

        // read bloom filters
        std::cerr << "Loading BFs:" << bvfile1 << " " << bvfile2 << std::endl;
        BF* bf1 = load_bf_from_file(bvfile1, *hashes, num_hash);
        BF* bf2 = load_bf_from_file(bvfile2, *hashes, num_hash);
        bf1->load();
        bf2->load();

        UncompressedBF *u_bf1 = uncompress_bf(bf1);
        UncompressedBF *u_bf2 = uncompress_bf(bf2);

        std::cerr << "Computing Sim..." << std::endl;
        uint64_t test = u_bf1->similarity(u_bf2, sim_type);
        std::cerr << "Done " << std::endl;
        std::cout << test << std::endl;
        std::tuple<uint64_t, uint64_t> sim = u_bf1->b_similarity(u_bf2);
        std::cout << "num_filter_bits: " <<  u_bf1->num_filter_bits() << " and: " << std::get<0>(sim) << " or: " << std::get<1>(sim) << std::endl;

        //uint64_t sim = bf1->similarity(bf2);
        //std::cout << bf1->size() << " " << sim << std::endl;
        //std::cout << "Similarity: " << sim << std::endl;
        //std::cout << "Difference: " << bf1->size() - sim << std::endl;
        //std::cout << "Ones: " << bf1->count_ones() << " " << bf2->count_ones() << std::endl;
        //std::cout << "Size: " << bf1->size() << std::endl;

        delete bf1;
        delete bf2;

    } else if (command == "convert") {
        std::cerr << "Converting..." << std::endl;
        convert_jfbloom_to_rrr(jfbloom_file, out_file);

    } else if (command == "hashes") {
        // construct a new hashpair
        construct_hashes(hashes_file, nb_hashes);
    } else if (command == "count") {
        int nh;
        HashPair* hp = get_hash_function(hashes_file, nh);
        std::cerr << "Cutoff Count: " << cutoff_count << std::endl;
        count(fasta_file, out_file, *hp, nh, bf_size, num_threads, cutoff_count);

    } else if (command == "build") {
        std::cerr << "Building..." << std::endl;
        std::vector<std::string> leaves = read_filter_list(leaves_file);
        //build_bt_from_jfbloom(leaves, out_file, parallel_level);
        dynamic_build(hashes_file, leaves, out_file, sim_type); //std::stoi(sim_type));

    } else if (command == "compress" || command == "compress-rrr") {
        std::cerr << "Compressing.." << std::endl;
        BloomTree* root = read_bloom_tree(bloom_tree_file, false);
        std::ifstream in(bloom_tree_file.c_str());
        std::string header;
        getline(in, header);
        std::vector<std::string> fields;
        SplitString(header, ',', fields);
        compress_bt_rrr(root);
        if (!root->is_split()) {
            write_compressed_bloom_tree(out_file, "rrr", root, fields[1]);
        } else {
            write_split_bloom_tree(out_file, "rrr", root, fields[2]);
        }
    } else if (command == "uncompress") {
        std::cerr << "Uncompressing.." << std::endl;
        BloomTree* root = read_bloom_tree(bloom_tree_file, false);
        std::ifstream in(bloom_tree_file.c_str());
        std::string header;
        getline(in, header);
        std::vector<std::string> fields;
        SplitString(header, ',', fields);
        uncompress_bt(root);
    } else if (command == "compress-roar") {
        std::cerr << "Compressing using Roaring bitarrays.." << std::endl;
        BloomTree* root = read_bloom_tree(bloom_tree_file, false);
        std::ifstream in(bloom_tree_file.c_str());
        std::string header;
        getline(in, header);
        std::vector<std::string> fields;
        SplitString(header, ',', fields);
        if (num_threads > 1) {
            compress_bt_roar_parallel(root);
        } else {
            compress_bt_roar(root);
        }
        if (!root->is_split()) {
            write_compressed_bloom_tree(out_file, "roar", root, fields[1]);
        } else {
            write_split_bloom_tree(out_file, "roar", root, fields[2]);
        }
    } else if (command == "compress-rrr-single") {
        std::cerr << "Compressing a single file (" << bvfile1 << ") using RRR.." << std::endl;
        int nh;
        HashPair* hp = get_hash_function(hashes_file, nh);
        auto bf = UncompressedBF(bvfile1, *hp, nh);
        bf.load();
        bf.compress_rrr();
    } else if (command == "compress-roar-single") {
        std::cerr << "Compressing a single file (" << bvfile1 << ") using Roaring bitarrays.." << std::endl;
        int nh;
        HashPair* hp = get_hash_function(hashes_file, nh);
        auto bf = UncompressedBF(bvfile1, *hp, nh);
        bf.load();
        bf.compress_roar();
    } else if (command == "split") {
        std::cerr << "Splitting.." << std::endl;
        BloomTree* root = read_bloom_tree(bloom_tree_file, false);
        std::ifstream in(bloom_tree_file.c_str());
        std::string header;
        getline(in, header);
        std::vector<std::string> fields;
        SplitString(header, ',', fields);
        SplitBloomTree* split_root = split_bt(root);
        write_split_bloom_tree(out_file, root->compression_name(), split_root, fields[1]);
    } else if (command == "rebuild") {
        std::cerr << "Rebuilding.." << std::endl;
        BloomTree* root = read_bloom_tree(bloom_tree_file, false);
        rebuild_bt(root);
    } else if (command == "bv-stats") {
        BloomTree* root = read_bloom_tree(bloom_tree_file, false);
        dump_bloom_tree_bit_vector_stats(root);
    } else if (command == "bv-dump") {
        BloomTree* root = read_bloom_tree(bloom_tree_file, false);
        dump_bloom_tree_bit_vectors(root,bv_dump_spacing);
    } else if (command == "append") {
        std::cerr << "Appending two single files together (" << bvfile1 << "," << bvfile2 << ")" << std::endl;
        BloomTree::caching = true; // make sure bloom filter caching is on
        if (BF_INMEM_LIMIT < 3)    // .. and that we will cache enough nodes
            BF_INMEM_LIMIT = 3;
        int nh;
        HashPair* hp = get_hash_function(hashes_file, nh);
        auto bf1 = UncompressedBF(bvfile1, *hp, nh);
        auto bf2 = UncompressedBF(bvfile2, *hp, nh);
        bf1.load();
        bf2.load();
        auto bf3 = bf1.append(bvfile3, (const BF*) &bf2);
        bf3->save();
    } else if (command == "compress-rrr-double") {
        std::cerr << "Compressing two single files together (" << bvfile1 << "," << bvfile2 << ") using RRR.." << std::endl;
        BloomTree::caching = true; // make sure bloom filter caching is on
        if (BF_INMEM_LIMIT < 3)    // .. and that we will cache enough nodes
            BF_INMEM_LIMIT = 3;
        int nh;
        HashPair* hp = get_hash_function(hashes_file, nh);
        auto bf1 = UncompressedBF(bvfile1, *hp, nh);
        auto bf2 = UncompressedBF(bvfile2, *hp, nh);
        bf1.load();
        bf2.load();
        auto bf3 = bf1.append(bvfile3, (const BF*) &bf2);
        bf3->compress_rrr();
    } else if (command == "compress-wraprrr-double") {
        std::cerr << "Compressing two single files together (" << bvfile1 << "," << bvfile2 << ") using WrapRRR.." << std::endl;
        BloomTree::caching = true; // make sure bloom filter caching is on
        if (BF_INMEM_LIMIT < 3)    // .. and that we will cache enough nodes
            BF_INMEM_LIMIT = 3;
        int nh;
        HashPair* hp = get_hash_function(hashes_file, nh);
        auto bf1 = UncompressedBF(bvfile1, *hp, nh);
        auto bf2 = UncompressedBF(bvfile2, *hp, nh);
        bf1.load();
        bf2.load();
        bf1.compress_wraprrr((const BF*) &bf2);
    }

    std::cerr << "(main program) Done." << std::endl;
}

