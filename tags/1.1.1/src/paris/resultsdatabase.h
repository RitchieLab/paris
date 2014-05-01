#ifndef __RESULTSDATABASE_H__
#define __RESULTSDATABASE_H__

#include <string>
#include <vector>
#include <sqlite3.h>
#include "database.h"

using namespace std;

class ResultsDatabase: public Database
{
public:

	void query(string sql);
	
	void InitTable(string table_name, string query_string);
	
private:

};

#endif
