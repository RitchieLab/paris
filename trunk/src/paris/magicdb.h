#ifndef PARIS_MAGIC_DB
#define PARIS_MAGIC_DB
#define USE_MAGIC 1


#ifdef USE_MAGIC
// #include <soci.h>
// #include <soci-sqlite3.h>
#include <string>
#include "resultsdatabase.h"

namespace Paris {

// using namespace soci;

struct ParisResults {
	void Open(const char *dbFilename) {
		if (resultsDB) {
			std::string cnxParam = "dbname="+std::string(dbFilename)+" timeout=10";
		}
	}

	void InitTable(const char *tablename, const char *query) {
		if (resultsDB) {
			std::string sql = std::string("DROP TABLE IF EXISTS ") + std::string(tablename);
		}
	}
	static ResultsDatabase db;
	static bool resultsDB;
};



}

#endif //USE_MAGIC

#endif //PARIS_MAGIC_DB

