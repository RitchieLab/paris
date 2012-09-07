#ifndef __DATABASE_H__
#define __DATABASE_H__

#include <string>
#include <vector>
#include <sqlite3.h>
#include "dbexcept.h"
#include "utility/types.h"

using namespace std;

class Database
{
public:
	Database();
	Database(string filename);
	~Database();
	
	bool open(string filename);
	vector<vector<string> > query(string query);
	void close();
	
protected:
	sqlite3 *database;
};

#endif
