#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string>

#include <tclap/CmdLine.h>
#include "bfcluster.h"

typedef struct Args {
    std::string hash_file;
    std::string leave_bf_list;
    std::string bloom_filter_path;
    std::string output_dir;
    std::string output_file;
    bool greedy_clustering;
    bool enable_debug;
}Args;

bool TclapParser(Args & args, int argc, char** argv){
	std::string version = "0.9";
	try {
		std::string desc = "bfcluster \n";
		TCLAP::CmdLine cmd(desc, ' ', version);
		TCLAP::ValueArg<std::string> arg_hash_file("f", "hash_file", "hash file", true, "", "file");
		TCLAP::ValueArg<std::string> arg_bloom_filter_path("p", "path", "directory of bloom filter files", true, "", "string");
		TCLAP::ValueArg<std::string> arg_leave_bf_list("l", "bf_list", "bloom filter filename list", true, "", "file");
		//TCLAP::ValueArg<std::string> arg_output_dir("d", "output_dir", "output directory, default is current working directory", false, ".", "string");
		TCLAP::ValueArg<std::string> arg_output_dir("o", "output_directory", "output directory", true, "", "string");

        TCLAP::ValueArg<std::string> arg_output_file("b", "bf_tree", "bloom tree topology file", false, "bf_tree_topology.txt", "string");
        
        TCLAP::SwitchArg arg_enable_debug("d", "debug", "output debug information", cmd, false);

        std::string greedy_clustering_string = "use greedy clustering instread of hierarchical clustering. \n";
        //TCLAP::SwitchArg arg_greedy_clustering("g", "greedy", greedy_clustering_string, cmd, false);

        cmd.add(arg_output_file);
        cmd.add(arg_output_dir);
        //cmd.add(arg_output_dir);
        cmd.add(arg_bloom_filter_path);
        cmd.add(arg_leave_bf_list);
        cmd.add(arg_hash_file);
		cmd.parse(argc, argv);

		args.hash_file = arg_hash_file.getValue();
		args.leave_bf_list = arg_leave_bf_list.getValue();
		args.bloom_filter_path = arg_bloom_filter_path.getValue();
        args.output_dir = arg_output_dir.getValue();
        args.output_file = arg_output_file.getValue();
        //args.greedy_clustering = arg_greedy_clustering.getValue();
        args.enable_debug = arg_enable_debug.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << "\n";
		abort();
	}
	return true;
}

void file_exists (const std::string& filename) {
    struct stat buffer;   
    if (stat (filename.c_str(), &buffer) != 0) 
        printf( "cannot access %s\n", filename.c_str() );
}

void dir_exists(const std::string& pathname){
    struct stat info;
    if( stat( pathname.c_str(), &info ) != 0 )
        printf( "cannot access %s\n", pathname.c_str() );
    else if(! (info.st_mode & S_IFDIR)) 
        printf( "%s is no directory\n", pathname.c_str() );
}

int main(int argc, char* argv[])
{
    Args args;
    TclapParser(args, argc, argv);

    if(args.output_dir[args.output_dir.size()-1] != '/')
        args.output_dir.push_back('/');

    if(args.bloom_filter_path[args.bloom_filter_path.size()-1] != '/')
        args.bloom_filter_path.push_back('/');

    file_exists(args.hash_file);
    file_exists(args.leave_bf_list);
    dir_exists(args.bloom_filter_path);
    dir_exists(args.output_dir);

    BFCluster* bfcluster = new BFCluster(args.hash_file, 
                                        args.leave_bf_list, 
                                        args.bloom_filter_path,
                                        'u',
                                        args.enable_debug);
    

    std::string topology_file = args.output_dir + args.output_file;

    //if(args.greedy_clustering){
        bfcluster->GreedyClustering(args.output_dir, topology_file);
    //}else{
    //    bfcluster->HierarchicalClustering(args.output_dir, topology_file);
    //}

    //return 0;

    bfcluster->BuildUncompressedFilters();

    delete bfcluster;
    return 0;
}