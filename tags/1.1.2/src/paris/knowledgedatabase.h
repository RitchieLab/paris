#ifndef __KNOWLEDGEDATABASE_H__
#define __KNOWLEDGEDATABASE_H__

#include <string>
#include <vector>
#include <sqlite3.h>
#include "database.h"
#include <map>

using namespace std;


struct DBGene{
	int geneID, start, end;
	string ensemblID, chrom;
};

struct DBKnowledge{
	int groupID, geneID;
	string name, desc, c;
};

class KnowledgeDatabase: public Database
{
public:

	void get_group(string group_type_id, uint& groupType, string& groupName);

	void load_alias(map<uint, string> & aliasLookup);

	void get_start_stop(vector<vector<int> > &values,string chrom, string popID);
	
	void get_genes(vector<DBGene>& genes ,string chrom, int popID);
	
	void get_knowledge(vector<DBKnowledge>& info, int group_type_id);
	
	void list_group_ids(ostream& os);
	
	void get_user_lookups(map<uint, string> & chromLookup, map<string,uint>& geneIDLookup);
	
	uint get_max_group_id();
	
private:

};




#endif
