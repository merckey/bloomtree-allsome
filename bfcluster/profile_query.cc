#include "profile_query.h"
ProfileQuery::~ProfileQuery(){
	if(refer_roarbf != NULL){
		delete refer_roarbf;
		refer_roarbf = NULL;
	}
}

void ProfileQuery::CompareQuery(std::string hash_filename, 
	std::string roar_filename, 
	std::string query_filename, 
	int times = 1000)
{
	hash_pair = get_hash_function(hash_filename, number_of_hashes);
	refer_roarbf = new RoarBF(roar_filename, *hash_pair, number_of_hashes);
	refer_roarbf->load();
	bf_size = refer_roarbf->num_filter_bits();
	int refer_ones = refer_roarbf->count_ones();
	std::cout << "Number of ones in refer RoarBF: " << refer_ones << std::endl;

    query_sequence_list = ReadQueries(query_filename);
    
    //----
    // time profile of bf intersection
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

	IntersectBitVector(times);
	
	end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    //std::time_t end_time = std::chrono::system_clock::to_time_t(end);
    std::cout << "bf intersection query elapsed time: " << elapsed_seconds.count() << "s\n";
    //----

    //----
    // time profile of kmer query
    start = std::chrono::system_clock::now();

	QueryKmers(times);
	
	end = std::chrono::system_clock::now();
    elapsed_seconds = end-start;
    //end_time = std::chrono::system_clock::to_time_t(end);
    std::cout << "kmer query elapsed time: " << elapsed_seconds.count() << "s\n";
    //----

}

void ProfileQuery::IntersectBitVector(int times){
	for(std::string query_sequence: query_sequence_list){
        qs.emplace_back(new RoarQueryInfo(query_sequence, *hash_pair, number_of_hashes, bf_size, true));
    }
    for(int i = 0; i < times; i++){
    	for(auto q: qs){
    		    RoarBF* intersect_bf = dynamic_cast<RoarBF*>(refer_roarbf->intersection_with("tmp.interset.bf.bv.roar", q->roarbf));
			    int counts = intersect_bf->count_ones();
			    delete intersect_bf;
    	}
    }
}

void ProfileQuery::QueryKmers(int times){
	QuerySet queryset;
	for(std::string query_sequence: query_sequence_list){
		queryset.emplace_back(new QueryInfo(query_sequence));
	}
	for(int i = 0; i < times; i++){
		for(auto q: queryset){
			for (const auto & m : q->query_kmers) {
        		bool result = refer_roarbf->contains(m);
        	}
		}
	}
}

int main(int argc, char* argv[])
{
	if(argc < 4){
		std::cout << "roar_file query_file times" << std::endl;
	}
	std::string hash_filename = "/gpfs/cyberstar/pzm11/backup/sbt/compressedSBT/SRR1186053.bf";
	std::string roar_filename = argv[1];
	std::string query_filename = argv[2];
	std::string timestring = argv[3];
	int exp_times = stoi(timestring);

	ProfileQuery* pq = new ProfileQuery();
	pq->CompareQuery(hash_filename, 
					roar_filename, 
					query_filename, 
					exp_times);
	delete pq;
	return 0;
}