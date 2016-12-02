#ifndef profile_query_H		// (prevent multiple inclusion)
#define profile_query_H

#include <chrono>
#include <ctime>
#include "bfquery.h"

class ProfileQuery: public BFQuery{
public:
	~ProfileQuery();
	void CompareQuery(std::string hash_filename, std::string roar_filename, std::string query_filename, int times);
	void IntersectBitVector(int times);
	void QueryKmers(int times);
protected:
	RoarBF* refer_roarbf;

};

#endif