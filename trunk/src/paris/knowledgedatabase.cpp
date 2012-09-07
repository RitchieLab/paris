#include "knowledgedatabase.h"
#include <sstream>
#include <iomanip>
#include <iostream>

void KnowledgeDatabase::load_alias(map<uint, string> & aliasLookup){

	sqlite3_stmt *statement;
	int id;
	string alias;
	string query = "SELECT DISTINCT alias, gene_id FROM region_alias NATURAL JOIN regions WHERE region_alias_type_id=1300 GROUP BY gene_id";

	if(sqlite3_prepare_v2(database, query.c_str(), -1, &statement, 0) == SQLITE_OK)
	{
		int result = 0;
		while(true)
		{
			result = sqlite3_step(statement);
			
			if(result == SQLITE_ROW)
			{
				id = sqlite3_column_int(statement, 1);
				alias = (char*)sqlite3_column_text(statement, 0);
				aliasLookup[id]=alias;
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

void KnowledgeDatabase::get_group(string group_type_id, uint& groupType, string& groupName){

	sqlite3_stmt *statement;
	string query = "SELECT group_type_id, group_type FROM group_type WHERE group_type_id=" + group_type_id;

	if(sqlite3_prepare_v2(database, query.c_str(), -1, &statement, 0) == SQLITE_OK)
	{
		int result = 0;
		while(true)
		{
			result = sqlite3_step(statement);
			
			if(result == SQLITE_ROW)
			{
				groupType = sqlite3_column_int(statement, 0);
				groupName = (char*)sqlite3_column_text(statement, 1);
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

void KnowledgeDatabase::list_group_ids(ostream& os){

	sqlite3_stmt *statement;
	string query = "SELECT DISTINCT group_type_id, group_type, download_date FROM group_type";

	if(sqlite3_prepare_v2(database, query.c_str(), -1, &statement, 0) == SQLITE_OK)
	{
		int result = 0;
		while(true)
		{
			result = sqlite3_step(statement);
			
			if(result == SQLITE_ROW)
			{
				int groupID = sqlite3_column_int(statement, 0);
				string name = (char*)sqlite3_column_text(statement, 1);
				string downloaddate = (char*)sqlite3_column_text(statement, 2);
				os<<setw(10)<<right<<groupID<<setw(30)<<right<<name<<setw(30)<<right<<downloaddate<<endl;
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

void KnowledgeDatabase::get_start_stop(vector<vector<int> > &values, string chrom, string popID){
	stringstream ss;
	ss << "SELECT DISTINCT start, stop FROM ld_blocks WHERE chromosome='" << chrom << "' AND population_id='" << popID << "'";
	string query = ss.str();
	
	vector<int> positions(2,0);
	sqlite3_stmt *statement;
	if(sqlite3_prepare_v2(database, query.c_str(), -1, &statement, 0) == SQLITE_OK)
	{
		int result = 0;
		while(true)
		{
			result = sqlite3_step(statement);
			
			if(result == SQLITE_ROW)
			{
				positions[0]=sqlite3_column_int(statement, 0);
				positions[1]=sqlite3_column_int(statement, 1);
				values.push_back(positions);
			}
			else
			{
				break;   
			}
		}
		sqlite3_finalize(statement);
	}
	
	string error = sqlite3_errmsg(database);
	if(error != "not an error"){// cout << query << " " << error << endl;
	cout << query << endl;
	cout << error << endl;
		throw DBExcept(error);	
	}
}


void KnowledgeDatabase::get_genes(vector<DBGene>& genes, string chrom, int popID){

	stringstream ss;
	ss << "SELECT gene_id, ensembl_id, chrom, description, start, end FROM regions NATURAL JOIN region_bounds WHERE chrom='" << chrom << "' AND population_id='" << popID << "' GROUP BY regions.gene_id";
	string query = ss.str();
	
	DBGene gene;
	sqlite3_stmt *statement;
	if(sqlite3_prepare_v2(database, query.c_str(), -1, &statement, 0) == SQLITE_OK)
	{
		int result = 0;
		while(true)
		{
			result = sqlite3_step(statement);
			
			if(result == SQLITE_ROW)
			{
				gene.geneID=sqlite3_column_int(statement, 0);
				gene.ensemblID=(char*)sqlite3_column_text(statement, 1);
				gene.chrom=(char*)sqlite3_column_text(statement, 2);
				gene.start=sqlite3_column_int(statement, 4);
				gene.end=sqlite3_column_int(statement, 5);
				genes.push_back(gene);
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

void KnowledgeDatabase::get_knowledge(vector<DBKnowledge>& info, int group_type_id){

	stringstream ss;
// 	ss << "SELECT gene_id, ensembl_id, chrom, description, start, end FROM regions NATURAL JOIN region_bounds WHERE chrom=" << chrom << " AND population_id=" << popID << " GROUP BY regions.gene_id";
	ss << "SELECT group_id, group_name, group_desc, gene_id, chrom FROM groups NATURAL JOIN group_associations NATURAL JOIN regions WHERE group_type_id="<<group_type_id <<" ORDER BY chrom, group_id";
	string query = ss.str();
	
	DBKnowledge knowledge;
	sqlite3_stmt *statement;
	if(sqlite3_prepare_v2(database, query.c_str(), -1, &statement, 0) == SQLITE_OK)
	{
		int result = 0;
		while(true)
		{
			result = sqlite3_step(statement);
			
			if(result == SQLITE_ROW)
			{
				knowledge.groupID = sqlite3_column_int(statement, 0);
				knowledge.name = (char*)sqlite3_column_text(statement, 1);
				knowledge.desc = (char*)sqlite3_column_text(statement, 2);
				knowledge.geneID = sqlite3_column_int(statement, 3);
				knowledge.c = (char*)sqlite3_column_text(statement, 4);
				info.push_back(knowledge);
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

void KnowledgeDatabase::get_user_lookups(map<uint, string> & chromLookup, map<string,uint>& geneIDLookup){

	sqlite3_stmt *statement;
	int geneID;
	string alias,chrom;
	string query = "SELECT DISTINCT alias, gene_id, chrom FROM region_alias NATURAL JOIN regions WHERE region_alias_type_id=1300 GROUP BY gene_id";
	if(sqlite3_prepare_v2(database, query.c_str(), -1, &statement, 0) == SQLITE_OK)
	{
		int result = 0;
		while(true)
		{
			result = sqlite3_step(statement);
			
			if(result == SQLITE_ROW)
			{
				geneID = sqlite3_column_int(statement, 1);
				alias = (char*)sqlite3_column_text(statement, 0);
				chrom = (char*)sqlite3_column_text(statement, 2);
				geneIDLookup[alias] = geneID;
				chromLookup[geneID] = chrom;
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


uint KnowledgeDatabase::get_max_group_id(){
	sqlite3_stmt *statement;
	uint max_id;
	string query =  "select max(group_id) from groups";
	if(sqlite3_prepare_v2(database, query.c_str(), -1, &statement, 0) == SQLITE_OK)
	{
		int result = 0;
		while(true)
		{
			result = sqlite3_step(statement);
			
			if(result == SQLITE_ROW)
			{
				max_id = sqlite3_column_int(statement, 0);
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
		
	return max_id;
}

