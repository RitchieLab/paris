/* 
 * File:   chromosome.h
 * Author: torstees
 *
 * Chromosome is responsible for setting up the features based on LD-Blocks
 * This structure is the sole party responsible for clean-up of feature*
 *
 * Created on January 5, 2010, 2:35 PM
 */

#ifndef _CHROMOSOME_H
#define	_CHROMOSOME_H
#include <vector>
#include <iostream>
#include <string>
#include "feature.h"
#include "gene.h"
#include <assert.h>
#include <map>
#include "knowledgedatabase.h"
// #include <soci.h>
// #include <soci-sqlite3.h>

namespace Paris {

// The zone size is a tradeoff between storage size (larger zones -> fewer
// zone:feature pair records) and lookup speed (smaller zones -> fewer to
// scan for matches against a given input position/region).
// For LD feature regions, 100kb errs towards lower memory usage (~1.03
// records per LD block). --atf 2014-09-24
#define FEATURE_ZONE_SIZE 100000


class Chromosome {

public:
	Chromosome();
	Chromosome(const char *chrID);
	Chromosome(const Chromosome& orig);

	virtual ~Chromosome();

	/**
	 * @brief Add a SNP to the chromosome.
	 * This is necessary for us to know what positions are present in the dataset
	 * Snps are stored position->rsID
	 */
	void AddSNP(uint rsID, uint pos);

	/**
	 * @brief Notify that a SNP (RS#) has been seen again on a later chromosome.
	 */
	bool AmbiguousSNP(uint rsID);

	/**
	 * @brief Add a feature to the chromosome, and attribute it to the gene, geneID
	 */
	void AddFeature(uint geneID, uint start, uint stop);

	/**
	 * @brief Add a gene to the chromosome
	 */
	void AddGene(const char *ensID, uint id, uint start, uint end, std::map<uint, std::string>& aliasLookup);
	/**
	 * @brief Loads all features associated with this chromosome
	 */
// 	uint LoadFeatures(soci::session& sociDB, const char *popID, uint &featureID);
	uint LoadFeatures(KnowledgeDatabase& knowDB, const char *popID, uint &featureID);

	/**
	 * @brief Load genes from the database and store them in genes
	 */
// 	uint LoadGenes(soci::session& sociDB, uint popID, uint geneExpansion, std::map<uint, std::string>& aliasLookup);
	uint LoadGenes(KnowledgeDatabase& knowDB, uint popID, uint geneExpansion, std::map<uint, std::string>& aliasLookup);

	/**
	 * @brief Merge Features into the genes.
	 * This must be done after both features and genes are loaded
	 */
	void MergeFeaturesIntoGenes(std::ostream& os);

	/**
	 * @brief Add values from the results
	 * This should be done after the genes and features have been loaded
	 */
	void AddValue(uint snpIndex, float pvalue);

	/**
	 * @brief return a pointer to the gene associated with geneID or NULL
	 */
	Gene *GetGene(uint geneID);

	/**
	 * @brief Initialize the bins (sorted by feature size)
	 */
	void InitBins(std::set<Feature*, SortByFeatureSize>& bins, std::set<Feature*>& singleFeatureBins, std::ostream& os);

	/**
	 * @brief returns the number of snps associated with the current chromosome
	 */
	uint SnpCount();

	/**
	 * @brief Returns number of features found in the chromosome
	 */
	uint FeatureCount();

	/**
	 * @brief returns all features that overlap with the region between start and stop
	 */
	std::set<Feature*> GetFeatures(uint start, uint stop);

	/**
	 * @brief returns rsID for SNP at position, pos
	 */
	uint GetSNP(uint pos);
	/**
	 * @brief returns SNPs between start and stop
	 */
	std::multimap<uint, uint> GetSnps(uint start, uint stop);

	/**
	 * @Brief Returns all SNPs (position -> rsID)
	 */
	std::multimap<uint, uint> &GetSNPs();

	/**
	 * @brief return the chromosome's ID
	 */
	std::string ID();

	/**
	 * @brief return the length of the chromosome
	 */
	uint Length();

private:
	std::string chrID;									///< Numerical ID
	std::vector<Feature*> features;					///< features associated with the chromosome
	std::multimap<uint, Feature*> featureZone;	///< zone(pos/100k) -> Feature Lookup
	std::map<uint, Gene*> genes;						///< genes

	std::multimap<uint, uint> positionSNPs;			///< position -> rsid
	std::multimap<uint, uint> snpPositions;			///< rsid -> position
	std::set<uint> badSNPs;							///< rsids which should be ignored
	std::set<uint> badPositions;					///< rsids which should be ignored
	uint length;											///< Length of the chromosome
};




}

#endif	/* _CHROMOSOME_H */

