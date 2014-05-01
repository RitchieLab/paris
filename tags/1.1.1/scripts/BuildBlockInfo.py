#!/usr/analysis/bin/EPD/bin/python2.6


import sys, numpy
import struct
import os
import sqlite3

tolerance			= 1


def InitDB(filename):
	cnx = sqlite3.connect(filename)
	dest = cnx.cursor()
	dest.execute("CREATE TABLE IF NOT EXISTS ld_blocks (chromosome VARCHAR(2), start INTEGER, STOP INTEGER, population_id VARCHAR(3));")
	cnx.commit()
	return cnx


class Block:
	def __init__(self, blockID, rsBeg, rsEnd):
		self.blockID 	= blockID
		self.begin 		= rsBeg
		self.end			= rsEnd

class Chromosome:
	def __init__(self, chromID):
		self.chromID 		= chromID
		self.length			= 0
		self.SNPs 			= dict()		# RS -> position
		self.dupeRSID 	= []
		self.blocks			= dict()

	def AddSNP(self, rsID, position):
		if rsID in self.SNPs:
			self.dupeRSID.append(rsID)
		else:
			if position > self.length:
				self.length = position
			self.SNPs[rsID] = position

	def PurgeDupes(self):
		for snp in self.dupeRSID:
			self.SNPs.erase(snp)
		self.rsIDs = self.SNPs.keys()

	def Reset(self):
		"""This allows us to do them all at once!"""
		self.blocks 		= dict()


	def PurifySNPs(self, differences, locations):
		"""Finds distance outliers and removes them from the final location list"""
		global tolerance
		a					= numpy.array(differences)
		mean				= numpy.mean(a)
		
		distances			= []
		#build the array of distances from the mean of each subset
		for i in xrange(0, len(differences)):
			t				= [a for a in differences]
			v				= t[i]
			t.remove(v)
			
			tt				= numpy.array(t)
			localMean		= numpy.mean(tt)
			dist			= mean - localMean
			distances.append(dist)
		
		od					= numpy.array(distances)
		oMean				= numpy.mean(od)
		oStd				= numpy.std(od)
		
		oMin				= oMean - (tolerance * oStd)
		oMax				= oMean + (tolerance * oStd)
		
		validMembers		= []
		validDiff			= []
		for i in xrange(0, len(differences)):
			dist			= distances[i]
			if dist >= oMin and dist <= oMax:
				#print i, len(locations), len(differences)
				validMembers.append(locations[i])
				validDiff.append(differences[i])
		return validMembers, validDiff
		
		
	# our blocks are defined only by SNPs inside our set of SNPs
	def DefineBlockRegion(self, blockID, snpSet):
		if len(snpSet) > 1:
			snps	= list(snpSet)
			snps.sort()
			#print snpSet, snps
			first 	= snps[0]
			second 	= snps[-1]
			
			
			if first in self.blocks:
				#there is no point in adding this one again, if we have already gotten it, unless we
				#find ours is larger than what we previously saw (right hand part of a sliding window)
				if second <= self.blocks[first].end:
					print "Skipping duplicate block: %s %s (%s)" % (first, second, self.blocks[first].end)
					return
			#print "Defining Block: %s %s-%s" % (self.SNPs[first], self.SNPs[first], self.SNPs[second])
			self.blocks[first] = Block(blockID, first, second)

	def ParseBlockFileRegionBased(self, blockFile, markerFile, chrom):
		global tolerance
		markers 					= []	#index  -> rsID map
		hapmapLoc					= []	#index  -> bp location -- We will use this to determine which SNPs to drop, if any
		hmSnpsInBlocks				= 0
		snpsRetainedInBlocks		= 0
		#print>>sys.stderr, "Loading file, %s" % (blockFile)
		blockCount = 0
		if os.path.exists(markerFile) and os.path.exists(blockFile):
			for line in open(markerFile):
				if line[0] != '#':
					words = line.split()
					markers.append(int(words[1].replace("rs", "")))
					hapmapLoc.append(int(words[2]))
			for line in open(blockFile):
				words = line.split()
				distances 			= []	# Differences in sequence
				snps				= []	# bp locations of the snps on the local sequence
				if len(words) > 1 and words[0] == "BLOCK":
					blockCount += 1
					for snp in words[3:]:
						rsid = markers[int(snp)-1]
						hLoc = hapmapLoc[int(snp)-1]
						if rsid in self.SNPs:
							d = abs(self.SNPs[rsid]-hLoc)
							distances.append(d)
							snps.append(self.SNPs[rsid])
					size = 0
					hmSnpsInBlocks+= len(snps)
					while len(snps) != size:
						size = len(snps)
						snps, distances		= self.PurifySNPs(distances, snps)
					snpsRetainedInBlocks += len(snps)
					if len(snps)> 0:
						if max(distances) - min(distances) > 500000:
							for i in xrange(len(distances)):
								print>>sys.stderr,  "%s\t - %s\t%s" % (i, snps[i], distances[i])
							print>>sys.stderr, "Oversized Block: ", max(distances) - min(distances)
						else:
							self.DefineBlockRegion(words[1], set(snps))
		#try:
		#	os.remove(blockFile)
		#except:
		#	pass
		#try:
		#	os.remove(markerFile)
		#except:
		#	pass
		#print>>sys.stderr, "%s blocks found. (%s)" % (blockCount, len(self.blocks))
		percKept				= 0.0
		if snpsRetainedInBlocks > 0 and hmSnpsInBlocks > 0:
			percKept			= float(snpsRetainedInBlocks)/float(hmSnpsInBlocks)
		print>>sys.stderr, "%s,%s,%s,%s,%s,%s" % (blockFile, blockCount, len(self.blocks), snpsRetainedInBlocks, hmSnpsInBlocks, percKept)
		#print>>sys.stderr, "%s out of %s (%s) snps were found inside blocks\nTolerance was %s standard deviations from mean" % (snpsRetainedInBlocks, hmSnpsInBlocks, percKept, tolerance)
	
	def ParseBlockFile(self, blockFile, markerFile):
		markers = []
		print "Loading file, %s" % (blockFile)
		blockCount = 0
		if os.path.exists(markerFile) and os.path.exists(blockFile):
			for line in open(markerFile):
				if line[0] != '#':
					words = line.split()
					markers.append(int(words[1].replace("rs", "")))
			for line in open(blockFile):
				words = line.split()
				if len(words) > 1 and words[0] == "BLOCK":
					snps = []
					blockCount += 1
					for snp in words[3:]:
						snps.append(markers[int(snp)-1])

					self.DefineBlock(words[1], set(snps))
		try:
			os.remove(blockFile)
		except:
			pass
		try:
			os.remove(markerFile)
		except:
			pass
		print "%s blocks found. (%s)" % (blockCount, len(self.blocks))

	def WriteToFile(self, file, popID):
		for b in self.blocks:
			block = self.blocks[b]
		#	print>>file, "%s,%s,%s,%s" % (self.chromID, popID, block.begin, block.end)

	def Commit(self, cur, popID):
		cur.execute("DELETE FROM ld_blocks WHERE chromosome='%s' and population_id='%s'" % (self.chromID, popID))
		#print "Committing %s blocks into chromosome %s, pop %s" % (len(self.blocks), self.chromID, popID)
		for b in self.blocks:
			block = self.blocks[b]
			#print "INSERT INTO ld_blocks VALUES ('%s', %s, %s, %s);" % (self.chromID, block.begin, block.end, popID)
			cur.execute("INSERT INTO ld_blocks VALUES ('%s', %s, %s, '%s');" % (self.chromID, block.begin, block.end, popID))

	def GatherBlockDetails(self, eth, stride, buffer):
		start = 0
		end 	= stride

		steps = stride - buffer
		length = self.length / 1000 #haploview limits are in kb

		while start < length:
			syscall = "java -Xms2000m -Xmx6500m -jar ~/Downloads/Haploview.jar -n -hapmapDownload -chromosome %s -startpos %s -endpos %s -panel %s -blockoutput GAB -release 22 -check -out BBI -memory 6128" % (self.chromID, start, end, eth)
			#print syscall
			rv = os.system(syscall)
			if rv == 0:
				self.ParseBlockFile("BBI.GABRIELblocks", "BBI.CHECK")

			start+=steps
			end = start + stride

		print "Chromosome %s completed. %s blocks found." % (self.chromID, len(self.blocks))

				
	# our blocks are defined only by SNPs inside our set of SNPs
	def DefineBlock(self, blockID, rsIDs):
		localRsIDs = set(self.SNPs.keys())
		validSNPs = rsIDs & localRsIDs

		if len(validSNPs) > 1:
			first 	= None
			second 	= None

			for snp in rsIDs:
				if snp in validSNPs:
					if first:
						second = snp
					else:
						first = snp
			if self.SNPs[first] in self.blocks:
				#there is no point in adding this one again, if we have already gotten it, unless we
				#find ours is larger than what we previously saw (right hand part of a sliding window)
				if self.SNPs[second] <= self.blocks[self.SNPs[first]].end:
					print "Skipping duplicate block: %s %s" % (self.SNPs[second], self.blocks[self.SNPs[first]].end)
					return
			#print "Defining Block: %s %s-%s" % (self.SNPs[first], self.SNPs[first], self.SNPs[second])
			self.blocks[self.SNPs[first]] = Block(self.SNPs[first], self.SNPs[first], self.SNPs[second])
def MakeAscii(str):
        s = ""
        for i in str:
                try:
                        i.encode("iso-8859-8")
                except UnicodeDecodeError:
                        i = " "
                        continue

                s += i
        return s
def ParseVariationsNew(db):
	c 					= db.cursor()
	c.execute("SELECT version FROM versions WHERE element = 'variations'")
	filename			= c.fetchone()[0]
	#for now, let's assume it's a hard path
	file				= open(MakeAscii(filename))

	print "Opening variations file, ", filename
	chromosomes = dict()
	data = file.read(4)
	completed = False
	while not completed:
		#read chromosome and the number of SNPs
		data = ""
		chr = ""
		role = ""
		chr = file.read(2)
		if chr != "":
			data = file.read(8)
			(snpCount, offset) = struct.unpack('II', data)
			chr = chr.strip()
			chromosome = Chromosome(chr)


			for i in xrange(0,snpCount):
				data = file.read(9)
				#print "-->> %s" % (len(data))
				(rsID, position, role) = struct.unpack("IIc", data)
				#print "%s\t%s\t%s" % (rsID, position, "%s" % role)
				chromosome.AddSNP(rsID, position)
			chromosomes[chr] = chromosome
			print "Chr: %s : %s, %s, %s" % (chr, snpCount, offset,chromosome.length)
		else:
			completed = True

	return chromosomes
	

def ParseVariationsOrig(filename):
	print "Opening variations file, ", filename
	file = open(filename)
	chromosomes = dict()

	completed = False
	while not completed:
		#read chromosome and the number of SNPs
		data = ""
		chr = ""
		chr = file.read(2)
		if chr != "":
			data = file.read(8)
			(snpCount, offset) = struct.unpack('II', data)
			chr = chr.strip()
			chromosome = Chromosome(chr)


			for i in xrange(0,snpCount):
				data = file.read(8)
				(rsID, position) = struct.unpack("II", data)
				chromosome.AddSNP(rsID, position)
			chromosomes[chr] = chromosome
			print "Chr: %s : %s, %s, %s" % (chr, snpCount, offset,chromosome.length)
		else:
			completed = True

	return chromosomes



if __name__ == '__main__':
	filename = None
	if len(sys.argv) > 3:
		dbFilename		 			= sys.argv[1]
		db = InitDB(dbFilename)
		
		chromosomes 				= ParseVariationsNew(db)
		#chromosomes				= ParseVariationsOrig("variations.bn")
		
		eth 						= sys.argv[3:]
		
		pathToData					= sys.argv[2]
				
		file = ("BlockData.txt", "w")
		jobs = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22']
		
		cur = db.cursor()
		for pop in eth:				#['CEU', 'YRI', 'JPT', 'CHB']:
			#print "Gathering LD Blocks for population: ", pop
			print>>sys.stderr, "chromosome,block_count,num_blocks,total_snps_used,total_snps_observed,perc_used"
			for chromosome in jobs:
				chromosomes[chromosome].Reset()
				#These are in kb...
				filename = "%s/Chromosome%s%s" % (pathToData, chromosome, pop)
				chromosomes[chromosome].ParseBlockFileRegionBased("%s.GABRIELblocks" % (filename), "%s.CHECK" % (filename), chromosome)
				#chromosomes[chromosome].GatherBlockDetails(pop, 50000,750)
				chromosomes[chromosome].WriteToFile(cur, pop)
				chromosomes[chromosome].Commit(cur, pop)
				db.commit()
			chromosomes['X'].Reset()
			chromosomes['X'].ParseBlockFileRegionBased("%s/ChromosomeXCEU_0_2689.GABRIELblocks" % pathToData, "%s/ChromosomeXCEU_0_2689.CHECK"%pathToData, "X")
			chromosomes['X'].WriteToFile(cur, pop)
			chromosomes['X'].Commit(cur, pop)

			chromosomes['X'].Reset()
			chromosomes['X'].ParseBlockFileRegionBased("%s/ChromosomeXCEU_2690_154000.GABRIELblocks"%pathToData, "%s/ChromosomeXCEU_2690_154000.CHECK"%pathToData, "X")
			chromosomes['X'].WriteToFile(cur, pop)
			chromosomes['X'].Commit(cur, pop)
			db.commit()	

	else:
		print >> sys.stderr, "Usage: db-filename path-to-data eth [eth]..."
		
