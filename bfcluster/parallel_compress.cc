#include "parallel_compress.h"
#include "tbb/task_scheduler_init.h"
#include "tbb/parallel_for.h"
#include <tclap/CmdLine.h>

bool enable_debug = false;


void ParallelApplyCompression(std::string bv_list[], size_t n, std::string _hash_filename, int _roar_compress){
	tbb::parallel_for(tbb::blocked_range<size_t>(0,n), ApplyCompression(bv_list, _hash_filename, _roar_compress), tbb::auto_partitioner() );
}

inline std::string GetLinuxFilePath(const std::string& str){
	int found = str.find_last_of("/");
	return str.substr(0,found+1);
}

inline std::string GetLinuxFileName(const std::string& str){
	int found = str.find_last_of("/");
	return str.substr(found+1);
}

inline std::string RemoveStars(const std::string& str){
	int found = str.find_last_of("*");
 	return str.substr(found+1);
}

/*split function*/
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        if (!item.empty()) {
            elems.push_back(item);
        }
    }
    return elems;
}

/*This split function only support char as delim, string as delim please boost split function*/
std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

std::vector<std::string> ReadBvFileList(const std::string & bv_list_filename, std::string & hash_filename){
    std::ifstream input(bv_list_filename);
    std::vector<std::string> bloom_filter_name_list;
    std::string bv_list_filepath = GetLinuxFilePath(bv_list_filename);

    if(!input.good()){
        std::cerr << "Error: error opening '"<<bv_list_filename<<"'. Bailing out." << std::endl;
        return bloom_filter_name_list;
    }
    if(enable_debug){
    	std::cout << "----reading bv list from file " << bv_list_filename << std::endl;
    }

    std::string line;
    bool first_line = true;
    //bool split_bv = false;
    while(std::getline(input, line).good()){
    	if(!line.empty()){
    		if(first_line){ // only one line contians ",", and it's in first line
                size_t n = std::count(line.begin(), line.end(), ',');
                if(n > 2){
                    std::cerr << "[pcompress] allsome split bloom filter detected." << std::endl;
                    //split_bv = true;
                }
    			int found = line.find_last_of(",");
    			if(found >= 0){
    				hash_filename = line.substr(found+1);
    				line = line.substr(0, found);
    			}else if(enable_debug){
    				std::cerr << "Error: can not find hash file info" << std::endl;
    			}
    			first_line = false;
    		}
    		std::string bloom_filter_name_line = RemoveStars(line);
    		std::vector<std::string> columns = split(bloom_filter_name_line, ',');
            for (std::string bloom_filter_name: columns){
                int found = bloom_filter_name.find_last_of("/");
        		if(found < 0){
        			bloom_filter_name = bv_list_filepath + bloom_filter_name;
        		}
        		bloom_filter_name_list.push_back(bloom_filter_name);
            }
    		if(enable_debug){
    			std::cout << line << std::endl;
    		}
    	}
    }

    if(enable_debug){
    	std::cout << "----finish reading bf list" << std::endl;
    }

    if(bloom_filter_name_list.size() == 0){
    	std::cerr << "Error: empty valid bloom filter " << std::endl;
    	return bloom_filter_name_list;
    }

    return bloom_filter_name_list;
}

typedef struct Args {
    std::string bv_list_filename;
    std::string compress_type;
    int thread_num;
}Args;

bool TclapParser(Args & args, int argc, char** argv){
    std::string version = "0.9";
    try {
        std::string desc = "squery \n";
        TCLAP::CmdLine cmd(desc, ' ', version);
        TCLAP::ValueArg<std::string> arg_bv_list_filename("l", "bv_list", "bv list filename", true, "", "file");
        TCLAP::ValueArg<std::string> arg_compress_type("c", "compress_type", "compress type(rrr/roar, rrr by default)", false, "rrr", "string");
        TCLAP::ValueArg<double> arg_thread_num("t", "thread", "thread number, default=8", false, 8, "int");
        cmd.add(arg_thread_num);
        cmd.add(arg_compress_type);
        cmd.add(arg_bv_list_filename);
        cmd.parse(argc, argv);

        args.bv_list_filename = arg_bv_list_filename.getValue();
        args.compress_type = arg_compress_type.getValue();
        args.thread_num = arg_thread_num.getValue();
    }
    catch (TCLAP::ArgException &e)
    {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << "\n";
        abort();
    }
    return true;
}

int main(int argc, char* argv[]){

    Args args;
    TclapParser(args, argc, argv);
    tbb::task_scheduler_init init(args.thread_num);
    std::cerr << "[pcompress] using " << args.thread_num << " threads." << std::endl;
	std::string bv_list_filename = args.bv_list_filename;
	// if file not provide path, then it should be in the same directory of bv_list_filename 

	int roar_compress = 0;
    if(args.compress_type == "roar"){
        roar_compress = 1;
        std::cerr << "[pcompress] performing roar compressing..." << std::endl;
    }

	std::string hash_filename;
	std::vector<std::string> bf_file_list = ReadBvFileList(bv_list_filename, hash_filename);
    int file_number = bf_file_list.size();
	std::string * bf_file_array = &bf_file_list[0];
	ParallelApplyCompression(bf_file_array, file_number, hash_filename, roar_compress);

	return 0;
}
