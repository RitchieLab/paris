/* 
 * File:   main.cpp
 * Author: torstees
 * 
 * Created on January 6, 2010, 12:54 PM
 */

#include "main.h"
#include <iostream>
#include <iomanip>
#include <string>
#include "timestamp.h"
#include "utility/strings.h"
#include "parislogs.h"
using namespace std;
using namespace Utility;

namespace Paris {


Main::Main() : action(ParisAction::NoAction) {
}

Main::Main(const Main& orig) : action(orig.action) {
}

Main::~Main() {
}


void Main::PrintBanner()  {
	cerr<<"          *******************************************************\n";
	stringstream ss;
	ss<<APPMAJOR<<"."<<APPMINOR<<"."<<APPBUGFIX<<"b ("<<BUILD_NUMBER<<") "<<BUILD_TYPE;
	cerr<<"          * paris        "<<setw(35)<<right<<ss.str()<<"    *\n"
		 <<"          *                                                     *\n"
	    <<"          *              "<<setw(35)<<right<<BUILD_DATE<<"    *\n";
#ifdef USE_MPI
	cerr<<"* This application is compiled to run on parallel computing systems using MPI\n";
##else
	cerr<<"* (serial)\n";
#endif
	cerr<<"          * Brian Yaspan, Jonathan Haines, Marylyn Ritchie and  *\n"
		 <<"          * Eric Torstenson                                     *\n"
	    <<"          * Comments/errors to software@ritchielab.psu.edu      *\n"
		 <<"          *******************************************************\n";
}

void Main::PrintHelp() {
	PrintBanner();
#ifdef USE_MPI
	cerr<<"usage: paris <configuration file> [ [command] ...] [ [parameter] ...]\n";
#else
	cerr<<"usage: paris <configuration file> [options]\n";
#endif
	cerr<<"Optional Commands Include:\n";
	cerr<<"  -S [--sample-config]                 -- Print sample configuration to std-out\n";
	cerr<<"  -I [--investigate-pathways] filename -- Generate gene reports for each gene found\n"
		 <<"                                          in the list of groups specified by filename\n";
	cerr<<"  -K [--list-knowledge]                -- Lists the knowledge bases with their keys\n"
		 <<"                                          (keys are used in the 'INCLUDE_GROUPS'\n"
		 <<"                                          configuration parameter)\n";
	cerr<<"\nOptional Parameters Include:\n";
	cerr<<"  -P [--list-populations]              -- Lists all available Population based LD\n"
		 <<"                                          boundary options\n";
}


StringArray Main::LoadGroupFile(const char *groups) {
	Utility::LineParser lp;
	Utility::FileToArray fta;
	lp.Parse(groups, &fta);
	return fta.strings;
}

int Main::ParseCmd(int curr, int argc, char **argv) {
	int nextCmd = curr+1;
/*	if (strcmp(argv[curr], "-G")==0 || strcmp(argv[curr], "--list-groups")==0) {
		action = ParisAction::ListGroups;
		if (nextCmd < argc) {
			cfg.AppendValue("GROUP_SEARCH_CRITERIA", argv[nextCmd++]);
		}
	}
*/	if (strcmp(argv[curr], "-K")==0 || strcmp(argv[curr], "--list-knowledge")==0) {
		action = ParisAction::ListGroups;
	}
	else if (strcmp(argv[curr], "-I")==0 || strcmp(argv[curr], "--investigate-pathways")==0) {
		action = ParisAction::InvestigatePathways;
		if (nextCmd < argc) {
			cfg.AppendValue("PATHWAY_INVESTIGATION_FILENAME", argv[nextCmd++]);
		}
		else {
			action = ParisAction::ParseError;
			cerr<<"-I (--investigate-groups) must be followed by the filename which holds the groups to be investigated\n";
		}
	}
	else if (strcmp(argv[curr], "-P")==0 || strcmp(argv[curr], "--list-population-ids")==0) {
		action = ParisAction::ListPopulationIDs;
	}

	else if (strcmp(argv[curr], "-S")==0 || strcmp(argv[curr], "--sample-config")==0)
		action = ParisAction::PrintSampleConfig;
	else {
		action = ParisAction::ParseError;
		cerr<<"Unknown argument: "<<argv[curr]<<"\n";
		return -1;
	}

	return nextCmd;
}

bool Main::ParseCmdLine(int argc, char **argv) {

	//Test the DB connection
#ifdef USE_MPI
	MPI::Init(argc, argv);
#endif
	if (argc < 2) {
		PrintHelp();
		return false;
	}
	PrintBanner();
	int i=1;
	if (argv[1][0] != '-')
		LoadConfiguration(argv[i++]);
	//Work out any other cmd line arguments
	for (; i<argc && i>0;) {
		i=ParseCmd(i, argc, argv);
	}
	if (action == ParisAction::ParseError) {
		return false;
	}
	if (action == ParisAction::PrintSampleConfig) {
		cfg.Init();
		cfg.Write(cout);
		return false;
	}
	parisApp.SetReportPrefix(cfg.GetString("REPORT_PREFIX").c_str());
	parisApp.ReportName(cfg.GetString("REPORT_NAME").c_str());

	parisApp.InitKnowledge(cfg.GetLine("SETTINGS_DB").c_str());

	//doLoadRegionAliases = cfg.GetBoolean("LOAD_ALL_ALIASES");

	cfg.ReportConfiguration(cout);

	return true;
}

AppConfiguration *Main::LoadConfiguration(const char *cfgFilename) {
	cfg.Init();
	cfg.SetValue("REPORT_PREFIX", Utility::ExtractBaseFilename(cfgFilename));
	if (cfgFilename)
		cfg.Parse(cfgFilename);
	cfg.ExecuteConfiguration();
	if (cfgFilename)
		configFilename=cfgFilename;
	else
		configFilename="";
	return &cfg;
}


string Main::GetReportPrefix() {
	string prefix = cfg.GetLine("REPORT_PREFIX");
	if (prefix == "")
		prefix = configFilename;
	return prefix;
}


void Main::RunCommands() {
	parisApp.SetReportPrefix(cfg.GetLine("REPORT_PREFIX").c_str());
	switch (action) {
		case ParisAction::PrintSampleConfig:
		{
			cfg.Write(cout);
			return;
		}
		case ParisAction::ListGroups:
			{
				parisApp.ListGroupIDs();
//				parisApp.ListGroupIDs(keywords);
			} return;
/*
		case ParisAction::ListPopulationIDs:
			parisApp.ListPopulationIDs();
			return;
		case ParisAction::ListMetaGroups:
			{
				InitGroupData();
				parisApp.ListMetaGroups(cout);
			}
			return;
		 */
		default: {}
	}

	multimap<string, uint> snps;
	vector<ParisApp::DataPoint> data;
	LoadDataPoints(data, snps);
	string ldPopulation					= cfg.GetLine("POPULATION");
	uint geneExpansion					= cfg.GetInteger("GENE_BOUNDARY_EXTENSION");
	StringArray knowledgeBases;
	cfg.GetLines("INCLUDE_KNOWLEDGE", knowledgeBases);
	string user_filename = cfg.GetLine("USER_PATHWAY_FILE");
	parisApp.InitKB(ldPopulation.c_str(), geneExpansion, knowledgeBases, user_filename);

	uint binSize							= cfg.GetInteger("BIN_SIZE");
	uint pCount								= cfg.GetInteger("P_COUNT");
	float pathSig							= cfg.GetDouble("PATHWAY_SIG_THRESH");
	float dataSig							= cfg.GetDouble("RESULTS_SIG_THRESH");
	
	parisApp.InitData(data, binSize);
	if (action == ParisAction::InvestigatePathways) {
		StringArray groups				= LoadGroupFile(cfg.GetString("PATHWAY_INVESTIGATION_FILENAME").c_str());
		StringArray::iterator itr		= groups.begin();
		StringArray::iterator end		= groups.end();

		while (itr != end) {
			parisApp.InvestigatePathway(pCount, dataSig, pathSig, itr++->c_str(), cfg.GetBoolean("SHOW_ALL_ASSOCIATED_PATHWAYS"));
		}
	}
	else {
		FileLogs::logger.OpenFinalReport();
		parisApp.RunAnalysis(pCount, dataSig, pathSig);
	}

}


void Main::LoadDataPoints(vector<ParisApp::DataPoint>& data, std::multimap<string, uint >& snps) {
	string filename	= cfg.GetLine("DATA_SOURCE");		///< chr, rs, pos
	uint chrColumn		= cfg.GetInteger("COL_CHROMOSOME")-1;	///< 0
	uint rsColumn		= cfg.GetInteger("COL_RSID")-1;			///< 1
	uint pvalColumn	= cfg.GetInteger("COL_PVALUE")-1;		///< 3

	char line[4096];
	ifstream file(filename.c_str());

	if (!file.good()) {
		cerr<<"Unable to open file, "<<filename<<". Aborting\n";
		exit(1);
	}

	while (file.good()) {
		file.getline(line, 4096);

		vector<string> tokens;
		if (Utility::TokenizeString(line, tokens, "\t, ") > 2) {
			string chr	= tokens[chrColumn];
			string rs	= tokens[rsColumn];
			//string p		= tokens[2];
			string pvalue = tokens[pvalColumn];

			// make sure the rs is valid
			// atoi() will say "7:50986:T_TA" == 7 but we need to be more careful than that
			if (rs.find("rs") == 0 || rs.find("RS") == 0)
				rs=rs.erase(0,2);
			std::string::iterator rsItr = rs.begin();
			while (rsItr != rs.end() && std::isdigit(*rsItr))
				rsItr++;
			if (!rs.empty() && rsItr == rs.end()) {
				uint rsID = atoi(rs.c_str());
				if(chr=="23")
					chr="X";
				else if(chr=="24")
					chr="Y";
				else if(chr=="25" || chr=="26") // 26 is also MT --atf 2014-09-23
					chr="MT";
				if ((chr == "X" || chr == "Y" || chr == "MT" || atoi(chr.c_str()) > 0)) { // loci are loaded from MT, so why not user data? --atf 2014-09-23
//	EST-RS>0		if (rs > 0.00) {
						data.push_back(ParisApp::DataPoint(rsID, chr.c_str(), atof(pvalue.c_str()) ));
						snps.insert(pair<string, uint>(chr, rsID));
//					}
				}
			}
		}
	}
	cout<<"\n"<<setw(35)<<filename<<" : "<<snps.size()<<" SNPs ";cout.flush();

	uint matches = parisApp.InitSNPs(snps, cfg.GetLine("VARIATION_FILENAME").c_str());
// 	cout<<" ("<<snps.size()<<" matches in our database )\n"; 
	cout<<" ("<<matches << " matches in our database )\n";
}


}
#ifndef TEST_APP
int main(int argc, char **argv) {
	string cfgFilename;

	Paris::Main app;					///<The application object
	if (!app.ParseCmdLine(argc, argv))
		exit(1);
	//Performs any commands
	app.RunCommands();

  	return EXIT_SUCCESS;
}
#endif //TEST_APP
