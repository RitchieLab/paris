#include "resultsdatabase.h"

void ResultsDatabase::query(string query_string){

	sqlite3_stmt *statement;

	if(sqlite3_prepare_v2(database, query_string.c_str(), -1, &statement, 0) == SQLITE_OK)
	{
		int result = 0;
		while(true)
		{
			result = sqlite3_step(statement);
			
			if(result == SQLITE_ROW)
			{
			}
			else
			{
				break;   
			}
		}
		sqlite3_finalize(statement);
	}
	
	string error = sqlite3_errmsg(database);
	if(error != "not an error")// cout << query << " " << error << endl;
		throw DBExcept(error);

}


void ResultsDatabase::InitTable(string table_name, string query_string){

	string drop_string = "DROP TABLE IF EXISTS " + table_name;
	query(drop_string);
	
	query(query_string);

}