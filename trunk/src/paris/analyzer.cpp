/* 
 * File:   analyzer.cpp
 * Author: torstees
 * 
 * Created on February 8, 2010, 10:30 AM
 */

#include "analyzer.h"
#include "feature.h"
#include <iostream>
#include <iomanip>

using namespace std;

namespace Paris {

Utility::Random &Analyzer::rnd = Utility::Random::globalGenerator;
std::pair<float, float> Analyzer::refinementThreshold = std::pair<float, float>(0.04, 0.06);
Analyzer::Analyzer() {
}

Analyzer::Analyzer(const Analyzer& orig) {
}

Analyzer::~Analyzer() {
}

Analyzer::Result Analyzer::Analyze(uint reportID, std::set<Gene*>& genes, std::map<uint, std::vector<Feature*> >& bins, uint permutationCount, bool verbose) {
	std::set<Gene*>::iterator itr = genes.begin();
	std::set<Gene*>::iterator end = genes.end();

	std::multiset<uint> binIDs;						///< This will tell us which bins to draw from
	std::set<uint> featureIDs;							///< Prevent from resampling one from the original source

	uint sigCount = 0;

	set<Feature*> usedFeatures;
	
	Analyzer::Result result(reportID);

	//Set up a pvalue of 1.0 to short circuit pathways with 0 significant features
	result.totPerms = permutationCount;
	result.sigPerms = permutationCount;
	result.geneCount = genes.size();

	while (itr != end) {
		Gene *gene = *itr;
		//binIDs will hold N copies of each of the bin IDs present, allowing us to
		//query how many of each bin we need to test
		gene->GetFeatureMap(binIDs, featureIDs);

		//Set up the number of features that have 1 or more significant pvalues in them
		uint sigMembers = gene->CountSignificantMembers(usedFeatures, result.simpleFeatures, result.complexFeatures, result.sigSimple, result.sigComplex);
		if (sigMembers > 0)
			result.sigGenes++;

		sigCount += sigMembers;
		itr++;
	}

	if (verbose)
		cerr<<"\tGenes: "<<setw(5)<<genes.size()<<"\tFeatures: "<<setw(5)<<binIDs.size()<<"\tSignificant: "<<setw(5)<<sigCount<<"..."; cerr.flush();

	//Short circuit to avoid testing situations that can never be less than 1.0
	if (sigCount > 0) {
		// check that bins are large enough for permutation testing for this set of features
		if(!CheckBinSize(binIDs, bins)){
			if(verbose)
				cerr << " no permutations run -- bins too small" << endl;
			// reset Result to show that no permutations were done
			result.totPerms=0;
			return result;
		}
	
		uint significantPTests=0;
		for (uint i=0; i<permutationCount; i++) {
			uint permutation = 0;
			//Draw without replacement-we don't want to allow these to be used in
			//permutations
			std::set<uint> untouchables = featureIDs;

			std::multiset<uint>::iterator binID = binIDs.begin();
			std::multiset<uint>::iterator binEnd = binIDs.end();
			while (binID != binEnd) {
				uint bIdx = *binID;
				vector<Feature*>& features = bins[bIdx];
				uint fIdx = this->rnd((int)features.size());
				Feature *f = features[fIdx];
				//Make sure this is "new" to the permutation
				if (untouchables.find(f->id) == untouchables.end()){
					if (f->CountSignificantMembers() > 0)
						permutation++;
					//Make sure that we can't draw this again
					untouchables.insert(f->id);
					// Only move to next bin draw when a feature was selected that had not been selected before
					binID++;
				}
			}

			//Permutation is only significant if it has more sig. features than the original
			if (permutation > sigCount)
				significantPTests++;
		}
		result.sigPerms = significantPTests;
		if (verbose)
			cerr<<" Sig PTests: "<<setw(8)<<significantPTests<<" -> ";

	}
	else
		if (verbose)
		  cerr<<" Sig PTests: "<<setw(8)<<0.0<<" -> ";
	if (verbose)
		cerr<<result.PValue()<<"\n";
		
	return result;
}

bool Analyzer::CheckBinSize(std::multiset<uint> &binIDs,
	std::map<uint, std::vector<Feature*> >& bins){
	
	std::multiset<uint>::iterator binEnd = binIDs.end();
	
	for(std::multiset<uint>::iterator binIter = binIDs.begin(); binIter != binEnd;
		binIter = binIDs.upper_bound(*binIter)){
		int count=int(binIDs.count(*binIter));
		// as we are sampling without replacement and excluding the original features
		// we need the bin to be twice the size of the count
		if(int(bins[*binIter].size()) - count*2 <= 0){
			return false;
		}
	}
	return true;
}

}

