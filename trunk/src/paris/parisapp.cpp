/* 
 * File:   parisapp.cpp
 * Author: torstees
 * 
 * Created on January 5, 2010, 1:29 PM
 */

#include "parisapp.h"
#include "utility/strings.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include "strings.h"
#include <algorithm>
#include "parislogs.h"
#include "bin.h"
#include "magicdb.h"

using namespace std;
using namespace Utility;

namespace Paris {

//ParisResults ParisResults::db;
ResultsDatabase ParisResults::db;
bool ParisResults::resultsDB = false;
int ParisApp::UserDefinedGroupType=7;

ParisApp::ParisApp() {
}

ParisApp::ParisApp(const ParisApp& orig) {
}

ParisApp::~ParisApp() {
}


void ParisApp::InitKnowledge(const char *dbFilename) {

	//Test the DB connection
	if (!Utility::FileExists(dbFilename)) {
		cerr<<"The database, "<<dbFilename<<", could not be found. Unable to continue.\n";
		exit(1);
	}
	try {
		knowDB.open(dbFilename);
		string sql = "SELECT version_id, ensembl_version, hapmap_version FROM version";
		vector<vector<string> > results = knowDB.query(sql);
		cout<<"\n------------------------- Dependency Versions ----------\n";
		cout<<setw(35)<<right<<"dbSNP: "<<results[0][0]<<"\n";
		cout<<setw(35)<<right<<"Ensembl: "<<results[0][1]<<"\n";
		cout<<setw(35)<<right<<"Hap Map LD: "<<results[0][2]<<"\n";

	} catch (DBExcept &e) {
		cerr<<"Problems were encountered trying to open the database, "<<dbFilename<<". Error: "<<e.what()<<"\n";
	}

}



uint ParisApp::InitBins(uint idealBinSize) {
	std::set<Feature*, SortByFeatureSize> binpool;
	std::set<Feature*> singleFeaturePool;

	//if (FileLogs::WriteBinDetails)
	FileLogs::logger.binDetails<<"Chrom,Gene ID,Start,Stop,RS ID,pvalue\n";

	std::map<std::string, Chromosome*>::iterator itr = chromosomes.begin();
	std::map<std::string, Chromosome*>::iterator end = chromosomes.end();

	while (itr != end) {
		itr->second->InitBins(binpool, singleFeaturePool, FileLogs::logger.binDetails);
		itr++;
	}
	uint activeBinCount = binpool.size();

	std::set<Feature*>::iterator sfItr = singleFeaturePool.begin();
	std::set<Feature*>::iterator sfEnd = singleFeaturePool.end();
	vector<Feature*> curBin;
	curBin.reserve(activeBinCount/4);
	uint spf = 0;

	while (sfItr != sfEnd) {
		Feature *f = *sfItr++;
		f->BinIndex(0);
		curBin.push_back(f);
		if (f->_begin == f->_end)
			spf++;

	}
	// it is (now) possible for a single-point feature to have multiple SNPs in it
	sfItr = binpool.begin();
	sfEnd = binpool.end();
	while (sfItr != sfEnd) {
		Feature *f = *sfItr++;
		if (f->_begin == f->_end)
			spf++;
	}

	//We want to mix the bins up, so that our bins aren't grouped the
	//same for each random seed
	vector<Feature*> binShuffler;
	binShuffler.insert(binShuffler.begin(), binpool.begin(), binpool.end());
	std::random_shuffle(binShuffler.begin(), binShuffler.end(), Bin::rnd);
	binpool.clear();
	binpool.insert(binShuffler.begin(), binShuffler.end());

	cerr<<"Single Point Features: "<<spf<<"\n";

	std::set<Feature*, SortByFeatureSize>::iterator fItr = binpool.begin();
	std::set<Feature*, SortByFeatureSize>::iterator fEnd = binpool.end();

	bins.clear();
	bins[0] = curBin;

	//By doing this, we make sure it's the smallest deviation from the desired size
	uint binCount = (uint)(((float)(activeBinCount))/(float)idealBinSize+0.5);
	if (binCount == 0)
		binCount = 1;
	uint binSize = (uint)(((float)(activeBinCount)+0.5) / binCount);
	int sparePoints = (activeBinCount)%binCount;



	cerr<<setw(10)<<"Bin"<<setw(10)<<"Feature Count"<<setw(12)<<"Avg Bin Size"<<"\n";
	cerr<<setw(10)<<"0"<<setw(10)<<curBin.size()<<setw(12)<<1<<"\n";

	for (uint i=1; i<=binCount; i++) {
		uint minFeatureSize = (uint)-1;
		uint maxFeatureSize = 0;
		curBin.clear();
		uint localBinSize = binSize;
		if (sparePoints-- > 0)
			localBinSize++;
		curBin.reserve(localBinSize);

		double sumFeatureSize = 0.0;

		//stringstream blockBuffer;
		//int superBlock = 0;

		for (uint n=0; n<localBinSize && fItr != fEnd; n++) {
			Feature *f = *fItr++;
			if (f->FeatureSize() < minFeatureSize)
				minFeatureSize = f->FeatureSize();
			maxFeatureSize = f->FeatureSize();
			sumFeatureSize += (double)maxFeatureSize;


			/*
			 if (f->FeatureSize() > 200) {
				Chromosome *chr = chromosomes[f->_chromosome];
				newFile<<"<TABLE  CELLSPACING=1 CELLPADDING=3 BORDER=1 RULES=ALL FRAME=HSIDES><TR><TH>Chromosome</TH><TH>Feature Beginning<TH><TH>Feature End</TH><TH>Feature Size</TH><TH>Link To SNPs</TH></TR>\n";

				newFile<<"\t<TR><TD>"<<f->_chromosome<<"</TD><TD>"<<f->_begin<<"</TD><TD>"<<f->_end<<"</TD><TD>"<<f->FeatureSize()<<"</TD><TD><A HREF=supersized-blocks.html#block_"<<superBlock<<"> Block "<<superBlock<<"</A></TD></TR>\n";
				blockBuffer<<"\n<A name=\"block_"<<superBlock<<"\"><H2>Block "<<superBlock<<"</H2></A><TABLE CELLSPACING=1 CELLPADDING=3 BORDER=1 RULES=ALL FRAME=HSIDES><TR><TH>Chromosome</TH><TH>RS ID</TH><TH>Position</TH><TH>P Value</TH><TH>Link</TH></TR>\n";
				superBlock++;
				map<uint, float> pvalues = f->GetPValues();
				map<uint, float>::iterator itr = pvalues.begin();
				map<uint, float>::iterator end = pvalues.end();
				while (itr != end) {
					blockBuffer<<"<TR><TD>"<<f->_chromosome<<"</TD><TD>"<<chr->GetSNP(itr->first)<<"</TD><TD>"<<itr->first<<"</TD><TD>"<<itr->second<<"</TD><TD><A HREF='http://uswest.ensembl.org/Homo_sapiens/Variation/Summary?source=dbSNP;v=rs"<<chr->GetSNP(itr->first)<<"'>rs"<<chr->GetSNP(itr->first)<<"</A></TD></TR>\n";
					itr++;
				}
				blockBuffer<<"</TABLE>\n";
			}*/
			
			f->BinIndex(i);
			curBin.push_back(f);
		}

		//newFile<<blockBuffer.str()<<"</HTML>\n";
		bins[i] = curBin;
		cerr<<setw(10)<<i<<setw(10)<<curBin.size()<<setw(12)<<(sumFeatureSize/(double)curBin.size())<<" ["<<minFeatureSize<<" - "<<maxFeatureSize<<"]"<<"\n";
	}


	double sumFeatureSize = 0.0;
	while (fItr != fEnd) {
		Feature *f = *fItr++;
		sumFeatureSize += (double)f->FeatureSize();
		f->BinIndex(binCount);
		bins[binCount].push_back(f);
	}
	if (bins[binCount].size() > binSize)
		cerr<<setw(10)<<binCount<<setw(10)<<bins[binCount].size()<<setw(12)<<(sumFeatureSize/(double)bins[binCount].size())<<"\n";



	return bins.size();
}

void ParisApp::InitData(vector<ParisApp::DataPoint>& data, uint binSize) {
	vector<ParisApp::DataPoint>::iterator itr = data.begin();
	vector<ParisApp::DataPoint>::iterator end = data.end();
	cerr<<"Populating Data\n";

	uint observedPoints = 0,
		 ignoredPoints = 0;
	while (itr != end) {
		chromosomes[itr->chrom]->AddValue(itr->rsID, itr->pvalue);
		if (itr->pvalue > 0.0 || !Feature::IgnorePValueOfZero) {
			observedPoints++;
			//chromosomes[itr->chrom]->AddValue(itr->rsID, itr->pvalue);
		}
		else
			ignoredPoints++;
		itr++;
	}
	cerr<<"Observed :"<<observedPoints<<"\tIgnored Points: "<<ignoredPoints<<"\n";
	cerr<<"Binning\n";
	//Perform Binning
	InitBins(binSize);
}

std::set<Pathway*> ParisApp::GetPathways(std::multimap<uint, uint>& groups) {
	std::multimap<uint, uint>::iterator itr = groups.begin();
	std::multimap<uint, uint>::iterator end = groups.end();
	std::set<Pathway*> pathways;
	while (itr != end) {
		if (knowledge.find(itr->first) != knowledge.end())
			pathways.insert(knowledge[itr->first]->GetPathway(itr->second));
		itr++;
	}
	return pathways;
}

void ParisApp::InvestigatePathway(uint permutationCount, float dataSig, float pathSig, const char *pathwayName, bool showAllPathways) {

	if (string(pathwayName).length() > 0) {
		std::map<uint, KnowledgeBase*>::iterator kbItr = knowledge.begin();
		std::map<uint, KnowledgeBase*>::iterator kbEnd = knowledge.end();
		cerr<<"Investigating Pathway "<<pathwayName<<":\n";
		Pathway *pathway = NULL;
		KnowledgeBase *kb = NULL;

		map<uint, Analyzer::Result> resultBuffer;

		while (kbItr != kbEnd && pathway == NULL) {
			kb = kbItr->second;
			pathway = kb->GetPathway(pathwayName);
			//kbItr->second->InvestigatePathways(bins, permutationCount, significance, scores, pathwayName);
			//kbItr->second->ReportResults(scores, report, true);
			kbItr++;
		}

		if (pathway) {
			ofstream snpReport(FileLogs::logger.GetFilename((pathway->Name() + "-snp-report").c_str(), "csv").c_str());
			pathway->GenerateSnpReport(snpReport,chromosomes, true);
			snpReport.close();
			ofstream file(FileLogs::logger.GetFilename(pathway->Name().c_str(), "html").c_str());
			file<<
	"<HTML><HEAD><TITLE>PARIS Pathway Breakdown "<<reportName<<" "<<pathwayName<<"</TITLE></HEAD>\n"<<
	"<BODY LINK=#4F2412 VLINK=#4F2412>\n";
	


			std::set<Gene*> genes = pathway->GetGenes();
			std::set<Gene*>::iterator gItr = genes.begin();
			std::set<Gene*>::iterator gEnd = genes.end();
			Analyzer analyzer;
			Analyzer::Result pathwayResult;
			if (resultBuffer.find(pathway->GroupID()) == resultBuffer.end()) {
				pathwayResult = analyzer.Analyze(pathway->GroupID(), genes, bins, permutationCount, false);
				resultBuffer[pathway->GroupID()] = pathwayResult;
			}
			else
				pathwayResult = resultBuffer[pathway->GroupID()];
			file<<"<H2>Pathway Investigation : ";
			if (reportName.length() > 0)
				file<<reportName<<" - ";
			if (kb->GetID() == 2)
				file<<"<A HREF='http://www.genome.jp/dbget-bin/www_bget?"<<pathwayName<<"' TEXT=A3C586>"<<pathwayName<<"</A> </H2><P>"<<pathway->Description()<<"<P>\n";
			else
				cerr<<"ID "<<kb->GetID()<<" </H2><P>"<<kb->GetName()<<"\n";
			file<<"<P><TABLE border=1 cellspacing=0>"
				 <<"<TR bgcolor=#CCCCCC border=1 styles='color:#A3C586;'><TH>Feature Count</TH>"
				 <<"<TH>Simple Features</TH><TH>Simple Features*</TH>"
				 <<"<TH>Complex Features</TH><TH>Complex Features*</TH>"
				 <<"<TH>P Value n/"<<permutationCount<<"</TH></TR>\n";
			file<<"\t<TD>"<<pathwayResult.complexFeatures+pathwayResult.simpleFeatures+pathwayResult.sigSimple+pathwayResult.sigComplex
				 <<"</TD><TD>"<<pathwayResult.simpleFeatures<<"</TD><TD>"<<pathwayResult.sigSimple
				 <<"</TD><TD>"<<pathwayResult.complexFeatures<<"</TD><TD>"<<pathwayResult.sigComplex<<"</TD><TD>"<<pathwayResult.GetPValue()<<"</TD></TR></TABLE>\n";
			file<<"\t<H2>Gene Breakdown for Pathway "<<pathwayName<<"</H2>\n<CENTER>\n";
			file<<"\t<P><TABLE border=1 cellspacing=0><TR bgcolor=#CCCCCC border=1><TH>Gene Name</TH><TH>Ensembl ID</TH><TH>Feature Count</TH><TH>Simple Feature</TH>"<<
				 "<TH>Simple Feature*</TH><TH>Complex Features</TH><TH>Complex Features*</TH><TH>Gene PValue n/"<<permutationCount<<"</TH><TH>Associated Pathways</TH></TR>\n";

			while (gItr != gEnd) {
				Gene *gene = *gItr;
				std::set<Gene*> singleGene;
				singleGene.insert(gene);
				Analyzer::Result result = analyzer.Analyze(gene->id, singleGene, bins, permutationCount, false);
				if (result.PValue() <= pathSig)
					file<<"\t<TR bgcolor=#E9E0DB border=1>";
				else
					file<<"\t<TR border=1>";
				file<<"<TD>"<<gene->Alias()<<"</TD><TD><A HREF='http://www.ensembl.org/Homo_sapiens/Gene/Summary?g="<<gene->EnsemblID()<<"'>"<<gene->EnsemblID()<<"</A></TD><TD>"<<gene->FeatureCount()
					<<"</TD><TD>"<<result.simpleFeatures<<"</TD><TD>"<<result.sigSimple<<"</TD>"
					<<"<TD>"<<result.complexFeatures<<"</TD><TD>"<<result.sigComplex<<"</TD>"
					<<"<TD>"<<result.GetPValue()<<"</TD><TD>";
					std::multimap<uint, uint> groups = gene->GetGroups();
					std::set<Pathway*> pathways = GetPathways(groups);
					std::set<Pathway*>::iterator itr = pathways.begin();
					std::set<Pathway*>::iterator end = pathways.end();
					int count = 0;


					while (itr != end) {
						if (count++ > 0)
							file<<" ";
						std::set<Gene*> pwGenes = (*itr)->GetGenes();

						Analyzer::Result pwResult;
						if (resultBuffer.find((*itr)->GroupID()) == resultBuffer.end()) {
							pwResult = analyzer.Analyze((*itr)->GroupID(), pwGenes, bins, permutationCount, false);
							resultBuffer[(*itr)->GroupID()] = pwResult;
						}
						else
							pwResult = resultBuffer[(*itr)->GroupID()];
						if (pathway->GroupID() != (*itr)->GroupID()) {
							if (pwResult.PValue() < pathSig){
								file<<"<B>"<<(*itr)->Name()<<"</B>";
							}
							else if (showAllPathways)
								file<<(*itr)->Name();
						}
						itr++;
					}
					file<<"</TD></TR>\n";

				gItr++;
			}
			file<<"</TABLE><P><I>*pvalues of "<<dataSig<<" or less.<P>Shaded rows indicate genes whose pvalues are "<<pathSig<<" or less.</P>";
			file << "<P><I>NA indicates where bins were not large enough for permutation testing.</BODY></HTML>\n";
		}
	}
	//Evaluate each group for statistical value (including the ptests)
}

void ParisApp::RunAnalysis(uint permutationCount, float datasetSignificance, float pathSignificance) {
	std::map<uint, KnowledgeBase*>::iterator kbItr = knowledge.begin();
	std::map<uint, KnowledgeBase*>::iterator kbEnd = knowledge.end();
	if (kbItr != kbEnd) {
		if (reportName.length() > 0)
			FileLogs::logger.finalReport<<reportName<<" ";
		kbItr->second->ResultsHeader(FileLogs::logger.finalReport);
		while (kbItr != kbEnd) {
			std::multiset<Analyzer::Result> scores;
			if (FileLogs::WriteDetailedFeatureReport)
				kbItr->second->DetailedReport(chromosomes, FileLogs::logger.detailedFeatureReport);
			kbItr->second->RunAnalysis(bins, permutationCount, pathSignificance, scores);
			kbItr->second->ReportResults(scores, FileLogs::logger.finalReport);
			kbItr++;
		}
		knowledge.begin()->second->ResultsFooter(FileLogs::logger.finalReport, datasetSignificance, pathSignificance);
	}

	//Evaluate each group for statistical value (including the ptests)
}
void ParisApp::InitKB(const char* popID, uint geneExpansion, Utility::StringArray& groups, string user_filename) {
		if (ParisResults::resultsDB) {
			string resultsDB = reportPrefix + "-results.sqlite";
			ParisResults::db.open(resultsDB);
			ParisResults::db.InitTable("knowledgebase",	"CREATE TABLE knowledgebase (kb_type INTEGER, kb_name VARCHAR(40))");
			ParisResults::db.InitTable("pathways",			"CREATE TABLE pathways (pathway_id INTEGER, kb_type INTEGER, pathway_name VARCHAR, pathway_description VARCHAR)");
			ParisResults::db.InitTable("genes",				"CREATE TABLE genes (gene_id INTEGER, ensembl_id VARCHAR, chrom VARCHAR, start INTEGER, end INTEGER)");
			ParisResults::db.InitTable("features",			"CREATE TABLE features (feature_id INTEGER, chrom STRING, start INTEGER, end INTEGER)");
			ParisResults::db.InitTable("feature_snps",	"CREATE TABLE feature_snps (feature_id INTEGER, rsid INTEGER, pos INTEGER, pvalue DECIMAL)");
			ParisResults::db.InitTable("gene_to_feature", "CREATE TABLE gene_to_feature (gene_id INTEGER, feature_id INTEGER)");
			ParisResults::db.InitTable("pathway_to_gene", "CREATE TABLE pathway_to_gene (pathway_id INTEGER, gene_id INTEGER)");
			ParisResults::db.InitTable("snps",				"CREATE TABLE snps (chrom STRING, rsid INTEGER, pos INTEGER)");
		}

		cerr<<"Initializing knowledge for population: "<<popID<<"\n";
		std::map<string, Chromosome*>::iterator itr = chromosomes.begin();
		std::map<string, Chromosome*>::iterator end = chromosomes.end();
		uint featureID = 0;

		//Load the alias lookup table
		map<uint, string> aliasLookup;

		knowDB.load_alias(aliasLookup);
		while (itr != end) {
			//First step, is to load the features
			itr->second->LoadFeatures(knowDB, popID, featureID);
			//Next, we want to Load the genes, and then merge the features into the genes
			itr->second->LoadGenes(knowDB, 0, geneExpansion, aliasLookup);
			itr++;
		}

		Utility::StringArray::iterator gitr = groups.begin();
		Utility::StringArray::iterator gend = groups.end();

		while (gitr != gend) {
			uint groupType=0;
			string groupName="";
				//Finally, we want to set up the group information and associated genes with the various groupings
			KnowledgeBase *kb;
			if(atoi(gitr->c_str()) != UserDefinedGroupType){
				knowDB.get_group(gitr->c_str(), groupType, groupName);
				if (ParisResults::resultsDB){
					stringstream ss;
					ss << "INSERT INTO knowledgebase VALUES(" << groupType << "," << groupName << ")";
					ParisResults::db.query(ss.str());
				}
				kb=new KnowledgeBase(groupType, groupName.c_str());

				cerr<<"Loading Knowledgebase "<<groupType<<" "<<groupName<<"\n";
				}
				else{
					groupType=UserDefinedGroupType;
                	UserKnowledgeBase *ub = new UserKnowledgeBase(groupType, "User-defined");
                    ub->ParseUserFile(user_filename);
                	kb=ub; 
				}
				kb->LoadKnowledge(knowDB, chromosomes, FileLogs::logger.kbReport);
				knowledge[groupType] = kb;
				++gitr;
		}
	
}

void ParisApp::ListGroupIDs() {
	knowDB.list_group_ids(cerr);
}

uint ParisApp::InitSNPs(std::multimap<string, uint>& allSNPs, const char *fn) {
	ifstream file(fn, ios::binary);
	if (!file.good()) {
		cerr<<"Unable to open file, "<<fn<<". Aborting\n";
		exit(1);
	}
	uint count = 0;
	std::set<uint> snps_seen;
	
	while (file.good()) {
		char label[3];
		file.read(label, 2);

		if (!file.eof()) {
			label[2]='\0';

			string chrID = label;
			chrID.erase(remove_if(chrID.begin(), chrID.end(), ::isspace), chrID.end());


			int snpCount=0, maxPosition=0;
			file.read((char*)&snpCount, 4);
			file.read((char*)&maxPosition, 4);

			std::multimap<string, uint>::iterator chromItr = allSNPs.lower_bound(chrID);
			std::multimap<string, uint>::iterator chromEnd = allSNPs.upper_bound(chrID);

			std::set<uint> snps;
			while (chromItr != chromEnd)  {
				snps.insert(chromItr->second);
				chromItr++;
			}

			Chromosome *newChrom = new Chromosome(label);

			if (file.good()) {
				cerr<<".";cerr.flush();
				for (int i=0; i<snpCount; i++) {
					int rs=0, pos=0;
					file.read((char*)&rs, 4);
					file.read((char*)&pos, 4);

					bool bad = false;
					if (true) { // TODO: optional global ambiguity check
						if (snps_seen.insert(rs).second == false) {
						//	cerr<<"bad snp: "<<rs<<"\n";
							bad = true;
							std::map<std::string, Chromosome*>::iterator cItr = chromosomes.begin();
							std::map<std::string, Chromosome*>::iterator cEnd = chromosomes.end();
							while (cItr != cEnd)
								cItr++->second->AmbiguousSNP(rs);
						}
					}

					// why load *all* SNPs on the chromosome if the user provided none? --atf 2014-09-23
					if (/*snps.size() == 0 ||*/ (snps.find(rs) != snps.end())) {
						//EST-RS>0 if (rs > 0) {
							newChrom->AddSNP(rs, pos);
							count++;
							if (bad)
								newChrom->AmbiguousSNP(rs);
						//}
					}
				}
			}
			chromosomes[chrID] = newChrom;
		}
	}

	return count;
}

void ParisApp::SetReportPrefix(const char *prefix) {
	reportPrefix = prefix;
}

void ParisApp::ReportName(const char *name) {
	reportName = name;
}
std::string ParisApp::ReportName() {
	return reportName;
}

}
