#include <iostream>
#include <string>
#include <thread>
#include <chrono>
#include "BF.h"
#include "BloomTree.h"
#include "tbb/blocked_range.h"
#include "tbb/tick_count.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

class ApplyCompression{
	std::string * const my_bvlist;
	HashPair* hashes;
	int hash_number;
	bool roar_compress;
public:
	void operator()(const tbb::blocked_range<size_t> & r) const{
		std::string * bv_list = my_bvlist;
		for(size_t i = r.begin(); i != r.end(); ++i){
			CompressRRR(bv_list[i]);
		}
	}

	ApplyCompression(std::string bv_list[], const std::string & hash_filename, int _roar_compress): my_bvlist(bv_list){
		hashes = get_hash_function(hash_filename, hash_number);
		if(_roar_compress != 0){
			roar_compress = true;
		}else{
			roar_compress = false;
		}
	}

	void CompressRRR(const std::string & bv_filename) const{
		//std::cout << "compress rrr: " << bv_filename << std::endl;
		if(roar_compress){
			std::string roar_filename = bv_filename + ".roar";
			struct stat buffer; 
			if(stat (roar_filename.c_str(), &buffer) == 0){
				std::cout << roar_filename << " exist, skip compressing" << std::endl;
				return;
			}
		}
		UncompressedBF * bf = new UncompressedBF(bv_filename, *hashes, hash_number);
		bf->load();
		if(roar_compress){
			bf->compress_roar();
		}else{
			bf->compress_rrr();
		}
		delete bf;
		//std::cout << "finish compressing " << bv_filename << std::endl;
	}
};
