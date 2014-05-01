/* 
 * File:   UserKnowledgeBase.cpp
 * Author: dudeksm
 * 
 * Created on November 8, 2011, 12:12 PM
 */

#include <fstream>
#include "userknowledgebase.h"

using namespace std;
// using namespace soci;

namespace Paris{

//UserKnowledgeBase::UserKnowledgeBase() {
//}
//
//UserKnowledgeBase::UserKnowledgeBase(const UserKnowledgeBase& orig) {
//}
//
//UserKnowledgeBase::~UserKnowledgeBase() {
//}
    
    
UserKnowledgeBase::UserKnowledgeBase(uint id, const char *name):KnowledgeBase(id, name){
    
}


uint UserKnowledgeBase::LoadKnowledge(KnowledgeDatabase& knowDB, std::map<std::string, 
        Chromosome*>& chromosomes, std::ostream& os){

    // can get the gene IDs from a map
//     rowset<row> rs = (sociDB.prepare <<"SELECT DISTINCT alias, gene_id, chrom FROM region_alias NATURAL JOIN regions WHERE region_alias_type_id=1300 GROUP BY gene_id");
    map<string, uint> geneIDLookup;
    map<uint, string> chromLookup;
    knowDB.get_user_lookups(chromLookup,geneIDLookup);
    
//         for (rowset<row>::const_iterator itr = rs.begin(); itr != rs.end(); ++itr) {
//                 row const& row = *itr;
// 		string alias = row.get<string>(0);
// 		uint geneID = row.get<int>(1);
//                 string chrom = row.get<string>(2);
// 		geneIDLookup[alias] = geneID;
//         chromLookup[geneID] = chrom;
// 	}
    
//     uint max_group_id;
    // get max group ID
//     sociDB << "select max(group_id) from groups", into(max_group_id);
	uint max_group_id = knowDB.get_max_group_id();
    
    std::map<std::string, std::vector<std::string> >::iterator pathiter;
    Chromosome *chrom = chromosomes.begin()->second;
    uint count=0;
    
    for(pathiter=user_pathways.begin(); pathiter != user_pathways.end(); pathiter++){
        // create the new pathway
        uint groupID = max_group_id++;
        Pathway *pathway = new Pathway(groupID, pathiter->first.c_str(), pathiter->first.c_str());
        pathways[groupID] = pathway;
	assert(pathwayLookup.find(pathiter->first) == pathwayLookup.end());
        pathwayLookup[pathiter->first] = pathway;
        
        // add genes to the pathway
        for(vector<string>::iterator geneIter=pathiter->second.begin(); geneIter != pathiter->second.end();
                geneIter++){

                map<string, uint>::iterator geneLookIter = geneIDLookup.find(*geneIter);
                if(geneLookIter == geneIDLookup.end()){
                    cerr<<"No match for gene " << *geneIter << "\n";
                    continue;
                }
                uint geneID = geneLookIter->second;
                string c = chromLookup[geneID];
        
        
                if (chrom->ID() != c)
                        chrom = NULL;
                
                if (chromosomes.find(c) != chromosomes.end())
                        chrom = chromosomes[c];
		else {
			cerr<<"We are having trouble finding chromosome '"<<c<<"'\n";
			std::map<std::string, Chromosome*>::iterator cItr = chromosomes.begin();
			while (cItr!=chromosomes.end())
				cerr<<"'"<<cItr++->first<<"' ";
			cerr<<"\n";
		}
		if (chrom) {
                        Gene *gene = chrom->GetGene(geneID);
			assert(gene);
			gene->AddGroup(id, groupID);

        		pathways[groupID]->AddGene(gene);
			os<<id<<"\t"<<pathiter->first<<"\t"<<pathiter->first<<"\t"<<gene->id<<"\t"<<gene->EnsemblID()<<"\t"<<gene->_chromosome<<"\t"<<gene->_begin<<"\t"<<gene->_end<<"\n";
			count++;
		}
        }
    }
	return count;       
    
}

// Parses user file and stores in map for loading knowledge
void UserKnowledgeBase::ParseUserFile(std::string filename){
        
        char line[4096];
        ifstream file(filename.c_str());
        
        if (!file.good()) {
		cerr<<"Unable to open file, "<<filename<<". Aborting\n";
		exit(1);
	}
        
        vector<string> gene_names;
        string pathway_name;
		while(file.getline(line, 4096)){
		vector<string> tokens;
                
                Utility::TokenizeString(line, tokens, "\t, ");
                
                if(tokens[0].compare("PATHWAY")==0){
                    pathway_name = tokens[1];
                    user_pathways[pathway_name] = gene_names;
		}
                else{
                    user_pathways[pathway_name].push_back(tokens[0]);
                }
	}
    
}

}
