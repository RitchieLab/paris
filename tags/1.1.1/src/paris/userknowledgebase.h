/* 
 * File:   UserKnowledgeBase.h
 * Author: dudeksm
 *
 * Created on November 8, 2011, 12:12 PM
 */

#ifndef USERKNOWLEDGEBASE_H
#define	USERKNOWLEDGEBASE_H

#include "knowledgebase.h"
namespace Paris {

class UserKnowledgeBase: public KnowledgeBase {
public:
    
    UserKnowledgeBase(uint id, const char *name);
    
    virtual uint LoadKnowledge(KnowledgeDatabase& knowDB, std::map<std::string, Chromosome*>& chromosomes, std::ostream& os);
    
    void ParseUserFile(std::string filename);
    
    
//    UserKnowledgeBase();
//    UserKnowledgeBase(const UserKnowledgeBase& orig);
//    virtual ~UserKnowledgeBase();
private:
    
    std::map<std::string, std::vector<std::string> > user_pathways;
    
};

}

#endif	/* USERKNOWLEDGEBASE_H */

