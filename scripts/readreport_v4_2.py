import sys,os,glob
import numpy as np
from scipy.stats.distributions import binom
import pysam
import operator
import time
import argparse
# module load pysam/intel/0.10.0
#python readreport_v3_2.py --infile bamfiles/filename.sorted.bam --ref reference/ref.name.fa

start_time = time.time()
parser = argparse.ArgumentParser() #argparse always need this first statment and then the following add arguments below
parser.add_argument('--ref','-r',required=True,help='give full path to reference file. Needs full path if not local!') #ref so this is a primer that inserts the reference file in the proper location
parser.add_argument('--infile','-i',help='input single bamfile. Needs full path if not local!')
parser.add_argument('--segment','-s',help='input single segment')
parser.add_argument('--qual','-q',type=int,default=25,help='phred quality cutoff (default is at 20)')
parser.add_argument('--cutoff','-c',type=float,default=0.01,help='minor variant frequency cutoff (default is at 0.01)')
parser.add_argument('--strain','-T',type=str,default='strain',help='need strain')
# parser.add_argument('--skipqual','-k',action='store_true',default=False,help='skip quality (usually for minion or pacbio)')
# parser.add_argument('--allowdel','-d',action='store_true',default=False,help='allow deletions (default no!)')
# parser.add_argument('--verbose','-v',action='store_true',default=False,help='Prints a read report file for every position')
parser.add_argument('--covercutoff','-C',type=int,default=200,help='coverage cutoff for consensus (default is at 200)')
args = parser.parse_args()

"""Below is an instance of 'self' variables. first make a class with object in this cass it is
class seqObject"""
class seqObject:
	def __init__(self, sampname,segment,ntpos,majornt,majorfreq,minornt,minorfreq,bcheck,
		A_nt,C_nt,G_nt,T_nt,gap_nt,totalcount):
		self.sampname = sampname
		self.segment = segment
		self.ntpos = ntpos
		self.majornt = majornt
		self.majorfreq = majorfreq
		self.minornt = minornt
		self.minorfreq = minorfreq
		self.bcheck = bcheck
		self.A = A_nt
		self.C = C_nt
		self.G = G_nt
		self.T = T_nt
		self.gap = gap_nt
		self.totalcount = totalcount
	def updateAA(self,aapos,majoraa,majorcodon):
		self.appos = aapos
		self.majoraa = majoraa
		self.majorcodon = majorcodon
	def printer(self,outfile):
		printlist = [self.sampname,self.segment,self.ntpos,self.majornt,
		self.majorfreq,self.minornt,self.minorfreq,self.bcheck,
		self.A,self.C,self.G,self.T,self.gap,self.totalcount,
		self.appos,self.majoraa,self.majorcodon]
		printlist = [str(x) for x in printlist]
		print>>outfile,','.join(printlist)
	def printer_no3(self,outfile):
		printlist = [self.sampname,self.segment,self.ntpos,self.majornt,
		self.majorfreq,self.minornt,self.minorfreq,self.bcheck,
		self.A,self.C,self.G,self.T,self.gap,self.totalcount]
		printlist = [str(x) for x in printlist]
		print>>outfile,','.join(printlist)

def seqUpdater(cigartuple,read,readidx,readq): #as long as you don't use updater again - good. cigartuple, read, readidx,readq
	"""look up cigar the aligned sequence may have additional bases that aren't in the reference or may be missing
	bases- CIGAR is a string to indicate if the reference and sequence match or mismatch due to deltions, insertions,
	soft-clipping or hard clipping. a CIGAR may look like (3M1I3M1D5M) which would stand for 3 match, 1 insertion, 3 match, 1 deletion, 5 match
	the tuple for cigartuple would look like [(0,3), (1,3), (0,3), (2,1), (0,5)]"""
	updatedcigarseq = []
	updatedseq = []
	updatedidx = []
	updatedreadq = []
	idxctr = 0
	for CIGIDX, (identifier,length) in enumerate(cigartuple):
		if identifier == 0: #match
			updatedcigarseq.extend('M'*length)
			for i in range(length):
				updatedseq.append(read[idxctr])
				updatedreadq.append(readq[idxctr])
				idxctr+=1
		elif identifier == 1: #insertion
			updatedcigarseq.extend('I'*length)
			for i in range(length):
				updatedseq.append(read[idxctr]) #with an insertion will insert idxctr
				updatedreadq.append(readq[idxctr])
				idxctr+=1
		elif identifier == 2: #deletion
			updatedcigarseq.extend('D'*length) #deletion will insert a little gap look up cigartuple
			for i in range(length):
				updatedseq.append('-')
				updatedreadq.append('-')

		elif identifier == 4: #softclip
			updatedcigarseq.extend('S'*length)
			for i in range(length):
				updatedseq.append(read[idxctr])
				updatedreadq.append(readq[idxctr])
				idxctr+=1

		elif identifier == 3: #skipped region added 4/9/2018
			pass

		elif identifier == 5: #hardclip. I put this code back. 8/30/2017
			pass
		else:
			sys.exit("cigartuple: invalid number encountered! %s in %s " % (identifier,cigartuple))
	#Because the Cigartuple doesn't come in a sequence like MMMIMMMMMDMMM - Tim wrote it out so that it matches with the sequence length
	# print cigartuple
	idxctr = 0
	last_delnt = 0
	last_insnt = 0
	for i,j,q in zip(updatedcigarseq,updatedseq,updatedreadq):
		if i == 'D':
			last_delnt+=1
			# print 'deletion',i,j,last_delnt,q
			updatedidx.append(last_delnt)
		elif i == 'I':
			# print 'insertion',i,j,last_insnt,q
			updatedidx.append(last_insnt)
			idxctr+=1
		else:
			# print i,j,readidx[idxctr],q
			updatedidx.append(readidx[idxctr])
			last_delnt = readidx[idxctr]
			last_insnt = readidx[idxctr]
			idxctr+=1
	assert len(updatedcigarseq) == len(updatedseq) == len(updatedidx) == len(updatedreadq) #check that they are all the same length assert function checks if True

	return (updatedcigarseq,updatedseq,updatedidx,updatedreadq)

"""Dictionaries: dictionary collection of many values (similar to list) but indexes for dictionaries
can use many different data types- not just integers. Indexes for dictionaries are called 'keys'. dictionaries are
typed with {} braces. Dict are not ordered like lists. B/c they are not ordered they can't be spliced
like lists can"""
def analyzer(isReverse,updatedOut): #takes output from above and puts in. isReverse is checking to make sure forward or reverse
	# updatedcigarseq,updatedseq,updatedidx,updatedreadq
	cig = updatedOut[0]
	seq = updatedOut[1]
	ntpos = updatedOut[2]
	qual = updatedOut[3]

	tempinsdict = {}
	if isReverse:
		# REVERSE_DICT
		for c,nt,pos,q in zip(cig,seq,ntpos,qual): #allowed to use zip because of same length take
			if pos == None and c == 'S':
				pass



			else:
				if c == 'I':
					if tempinsdict.has_key(pos):
						tempinsdict[pos].append(nt)
					else:
						tempinsdict[pos] = [nt]
				elif q >= args.qual: #quality passes threshold
					if REVERSE_DICT[pos].has_key(nt): #position populated
						REVERSE_DICT[pos][nt] = REVERSE_DICT[pos][nt] + 1
					else:
						REVERSE_DICT[pos][nt] = 1
					if CONSENSUS_DICT[pos].has_key(nt):
						CONSENSUS_DICT[pos][nt] = CONSENSUS_DICT[pos][nt] + 1
					else:
						CONSENSUS_DICT[pos][nt] = 1
		if tempinsdict:
			for ntpos in tempinsdict:
				fullnt = ''.join(tempinsdict[ntpos])
				if INSERTION_DICT.has_key(ntpos):
					pass
				else:
					INSERTION_DICT[ntpos] = {}

				if INSERTION_DICT[ntpos].has_key(fullnt):
					INSERTION_DICT[ntpos][fullnt] = INSERTION_DICT[ntpos][fullnt] + 1
				else:
					INSERTION_DICT[ntpos][fullnt] = 1
	else:
		#FORWARD_DICT
		for c,nt,pos,q in zip(cig,seq,ntpos,qual):
			if pos == None and c == 'S': #c == soft clip
				pass

			else:
				if c == 'I': #insertion
					if tempinsdict.has_key(pos):
						tempinsdict[pos].append(nt)
					else:
						tempinsdict[pos] = [nt]
				elif q >= args.qual:
					if FORWARD_DICT[pos].has_key(nt):
						FORWARD_DICT[pos][nt] = FORWARD_DICT[pos][nt] + 1
					else:
						FORWARD_DICT[pos][nt] = 1

					if CONSENSUS_DICT[pos].has_key(nt):
						CONSENSUS_DICT[pos][nt] = CONSENSUS_DICT[pos][nt] + 1
					else:
						CONSENSUS_DICT[pos][nt] = 1
		if tempinsdict:
			for ntpos in tempinsdict:
				fullnt = ''.join(tempinsdict[ntpos])
				if INSERTION_DICT.has_key(ntpos):
					pass
				else:
					INSERTION_DICT[ntpos] = {}

				if INSERTION_DICT[ntpos].has_key(fullnt):
					INSERTION_DICT[ntpos][fullnt] = INSERTION_DICT[ntpos][fullnt] + 1
				else:
					INSERTION_DICT[ntpos][fullnt] = 1



def binomCheck(ntpos): #checking the nt position for both forward and reverse mate pairs
	fordict = FORWARD_DICT[ntpos]
	revdict = REVERSE_DICT[ntpos]


	if not fordict or not revdict: 	 #if either forward dict or reverse dict is empty..
		accept = False

	else:
		topF =sorted(fordict, key=fordict.get, reverse=True)[:2] #dictionaries can be in any particular order, where lists indexes start at 0
		topR =sorted(revdict, key=revdict.get, reverse=True)[:2]

		if len(topF) == 1 or len(topR) == 1:
			accept = False
			print 'unequal minor variant count in forward/reverse %d' % ntpos #%s= string %d = number variable ntpos has already been defined
			print 'forward',fordict
			print 'reverse',revdict
			NOTESLIST.append('take a closer look at, only one minorvar %d' % ntpos)
			NOTESLIST.append(fordict)
			NOTESLIST.append(revdict)
		else:
			f_majornt = topF[0] #this will be the major because it will be the first highest
			f_minornt = topF[1] #this will be the minor variant because it will be the second highest, remember that python starts numbering at 0

			r_majornt = topR[0]
			r_minornt = topR[1]

			if f_majornt != r_majornt or f_minornt != r_minornt:
				print 'binom not equal'
				NOTESLIST.append('take a closer look at %d' % ntpos ) #this will be added to the notelist for what went wrong exactly
				NOTESLIST.append([f_majornt,r_majornt,f_minornt,r_minornt])
				accept = False
			else:
				forwardMajorCount = fordict[f_majornt]
				forwardMinorCount = fordict[f_minornt]

				reverseMajorCount = revdict[r_majornt]
				reverseMinorCount = revdict[r_minornt]
				ALPHA = 0.05 #to check and make sure that it is signficiant, or in a significant number of reads

				pforward = 1 - binom.cdf( forwardMinorCount, (forwardMajorCount + forwardMinorCount), args.cutoff) #calculating the p value
				preverse = 1 - binom.cdf( reverseMinorCount, (reverseMajorCount + reverseMinorCount), args.cutoff)
				if pforward <= ALPHA/2 and preverse <= ALPHA/2:
					accept = True
				else:
					accept = False
	return accept

def returnCodon(ntobject,seqlist): #0-index
	aminoacid = {'ttt': 'F', 'tct': 'S', 'tat': 'Y', 'tgt': 'C',
	'ttc': 'F', 'tcc': 'S', 'tac': 'Y', 'tgc': 'C',
	'tta': 'L', 'tca': 'S', 'taa': '*', 'tga': '*',
	'ttg': 'L', 'tcg': 'S', 'tag': '*', 'tgg': 'W',
	'ctt': 'L', 'cct': 'P', 'cat': 'H', 'cgt': 'R',
	'ctc': 'L', 'ccc': 'P', 'cac': 'H', 'cgc': 'R',
	'cta': 'L', 'cca': 'P', 'caa': 'Q', 'cga': 'R',
	'ctg': 'L', 'ccg': 'P', 'cag': 'Q', 'cgg': 'R',
	'att': 'I', 'act': 'T', 'aat': 'N', 'agt': 'S',
	'atc': 'I', 'acc': 'T', 'aac': 'N', 'agc': 'S',
	'ata': 'I', 'aca': 'T', 'aaa': 'K', 'aga': 'R',
	'atg': 'M', 'acg': 'T', 'aag': 'K', 'agg': 'R',
	'gtt': 'V', 'gct': 'A', 'gat': 'D', 'ggt': 'G',
	'gtc': 'V', 'gcc': 'A', 'gac': 'D', 'ggc': 'G',
	'gta': 'V', 'gca': 'A', 'gaa': 'E', 'gga': 'G',
	'gtg': 'V', 'gcg': 'A', 'gag': 'E', 'ggg': 'G'}
	ntpos = ntobject.ntpos
	ntpos = int(ntpos) - 1 #offset by one to increment by one

	if ntpos % 3 == 0:# first position
		fp = ntpos#,ntpos+1,ntpos+2]
		p2 = ntpos+1
		p3 = ntpos+2

		p1nt = seqlist[fp].majornt #this will print out the codon, and capatalize the nt position of the codon that you are looking at
		p2nt = seqlist[p2].majornt.lower()
		p3nt = seqlist[p3].majornt.lower()

		# altcodon = [minor.lower(),p2nt,p3nt]
	elif ntpos % 3 == 1: #second
		p1 = ntpos-1#,ntpos,ntpos+1]
		fp = ntpos
		p3 = ntpos+1

		p1nt = seqlist[p1].majornt.lower()
		p2nt = seqlist[fp].majornt
		p3nt = seqlist[p3].majornt.lower()

		# altcodon = [p1nt,minor.lower(),p3nt]

	elif ntpos % 3 == 2: #third
		# fp = ntpos-2#,ntpos-1,ntpos]
		p1 = ntpos-2
		p2 = ntpos-1
		fp = ntpos

		p1nt = seqlist[p1].majornt.lower()
		p2nt = seqlist[p2].majornt.lower()
		p3nt = seqlist[fp].majornt
		# altcodon = [p1nt,p2nt,minor.lower()]
	# print codon
	aapos = ntpos/3 + 1 #amino acid position number
	# ref_codon = refseq[fp:fp+3]
	majorcodon = ''.join([p1nt,p2nt,p3nt])

	############## DELETE BELOW ME ############################# 8/30/2017
	# if '-' in majorcodon: #check to see if there is a gap in the major codon
	# 	majoraa = ''
	# else:
	# 	majoraa = aminoacid[majorcodon.lower()]
	############## DELETE ABOVE ME #############################


	#added a catch-all 8/30/2017
	try:
		majoraa = aminoacid[majorcodon.lower()]

	except KeyError:
		#This will happen because of a deletion (-) or no coverage (N) in the codon
		majoraa = ''




	ntobject.updateAA(aapos,majoraa,majorcodon)


def ensure_dir(f):
	d = os.path.dirname(f)
	if not os.path.exists(d):
		os.makedirs(d)


def read_fasta(fp):
	name, seq = None, []
	for line in fp:
		line = line.rstrip()
		if line.startswith(">"):
			if name: yield (name, ''.join(seq))
			name, seq = line, []
		else:
			seq.append(line)
	if name: yield (name, ''.join(seq))


def open_fasta(filename):
	# filename = '../FILES/reference/'+someref
	segdict = {}
	with open(filename) as fp:
		for name, seq in read_fasta(fp):
			segdict[name[1:]] = seq
	return segdict



def printer(outfile, sampname, SEGMENT, ntpos, bcheck, hasMinor): #this is how the output data looks 'pretty'
	ref_nt = REF_DICT[SEGMENT][ntpos]
	tempd = CONSENSUS_DICT[ntpos]
	A_nt = 0
	C_nt = 0
	G_nt = 0
	T_nt = 0
	gap_nt = 0
	if tempd.has_key('A'):
		A_nt = tempd['A']
	if tempd.has_key('C'):
		C_nt = tempd['C']
	if tempd.has_key('G'):
		G_nt = tempd['G']
	if tempd.has_key('T'):
		T_nt = tempd['T']
	if tempd.has_key('-'):
		gap_nt = tempd['-']
	totalcount = sum(tempd.values())
	majornt = 	max(tempd,key=tempd.get)
	majorfreq = tempd[majornt]/float(totalcount)
	ntpos = ntpos + 1 #change from 0-index to 1-index
	if hasMinor:
		minornt = sorted(tempd, key=tempd.get, reverse=True)[:2][1]
		minorfreq = tempd[minornt]/float(totalcount)
		x = seqObject(sampname,SEGMENT,ntpos,majornt,majorfreq,minornt,minorfreq,bcheck,
			A_nt,C_nt,G_nt,T_nt,gap_nt,totalcount)
		# printlist = [str(x) for x in printlist]
		# print>>outfile,','.join(printlist)
	else:
		x = seqObject(sampname,SEGMENT,ntpos,majornt,majorfreq,'','','',
			A_nt,C_nt,G_nt,T_nt,gap_nt,totalcount)
		# printlist = [sampname,SEGMENT,ntpos,majornt,majorfreq,'','','',
		# 	A_nt,C_nt,G_nt,T_nt,gap_nt,totalcount]
		# printlist = [str(x) for x in printlist]
		# print>>outfile,','.join(printlist)
	return(x)




if __name__ == '__main__': #this will allow this module to be imported from another module and used
	NOTESLIST = [] #abnormal prints out what went wrong continues to add to the notelist as it works through script below
	STRAIN = args.strain.upper() #arsparse gives primers (can be expanded if needed)
	# STRAIN = 'h3n2'

	##################### USER DEFINED ########################
		# MAKE SURE THE / IS AT THE END OF FULLVARLIST_DIR!!!
	FULLVARLIST_DIR = 'FILES/fullvarlist/' #change output 8/30/2017
	# ensure_dir('../FILES/consensus/')
	ensure_dir('FILES/') #I don't think you can make a subdirectory without making the parent directory first 8/30/2017
	ensure_dir(FULLVARLIST_DIR) #this looks for the already made directory, and if not made then it makes it on its own
	bampath = '../FILES/bamfiles/' #specify the path and it will go in and do a batch file
	##########################################

	REF_FILE = args.ref #this is in the primers above when we go --ref
	REF_DICT = open_fasta(REF_FILE) #open_fasta is defined above

	SAMPLELIST = []
	if args.infile is None: #Batch files if you don't specify it will go into bamfile directory
		for infile in glob.glob( os.path.join(bampath, '*.bam') ): #glob is a pathway
			SAMPLELIST.append(infile) #appending all file names to samplelist
	else:
		infile = args.infile #Examining a single input file!
		# SAMPLENAME = infile.split('/')[-1].split('.')[0]
		SAMPLELIST.append(infile) #append is used to add one element samplist is appending the infile (bamfile)
		# samplebamdict[SAMPLENAME] = pysam.AlignmentFile(infile, "rb",check_header=False,check_sq=False )

	if args.segment is not None: #All segments relies on segment from reference
		REF_DICT = dict((key,value) for key, value in REF_DICT.iteritems() if key == args.segment)
		"""Dictionaries are designed so that they have a key and each key has a value assigned to them. dict.iteritems()
		will return an iterator () over the dict key,value pairs - I am guessing that the iter part of this means that it is going through all
		key,value pairs"""

	for sampctr, SAMPLEFILE in enumerate(SAMPLELIST):
		SAMPLENAME = SAMPLEFILE.split('/')[-1].split('.')[0] #makes assumption of name
		BAMFILE = pysam.AlignmentFile(SAMPLEFILE, "rb",check_header=False,check_sq=False)
		sys.stdout.write('**Analyzing sample %s [%d out of %d samples]** \n' %  (SAMPLENAME,sampctr+1,len(SAMPLELIST)))
		for segidx,SEGMENT in enumerate(REF_DICT): #enumerate segidx if you print will give you the index

			VARLIST_OUT = open("%s%s.%s.%s.%s.snplist.csv" % (FULLVARLIST_DIR,SAMPLENAME,STRAIN,SEGMENT,str(args.cutoff)),'w') #writing the file name for the snplist
			HEADER = 'name,segment,ntpos,major,majorfreq,minor,minorfreq,binocheck,A,C,G,T,-,totalcount,aapos,majoraa,majorcodon,minoraa,minorcodon' #header for the snplist
			print>>VARLIST_OUT,HEADER	#print out header to that file
			SEQLIST = []
			sys.stdout.write('Examining segment: %s  [%d/%d]\n' % (SEGMENT,segidx+1,len(REF_DICT) ))
			sys.stdout.write('\tBegin populating forward and reverse dictionaries ...')
			SEGLEN = len(REF_DICT[SEGMENT]) #length of ref seg
			FORWARD_DICT = {} #these are cap making fresh for dict building dictionaries see more above
			REVERSE_DICT = {}
			INSERTION_DICT = {}
			CONSENSUS_DICT = {}
			for idx in range(SEGLEN):
				FORWARD_DICT[idx] = {}
				REVERSE_DICT[idx] = {}
				CONSENSUS_DICT[idx] = {}
			counter = 0
			try:
				for read in BAMFILE.fetch(SEGMENT):
					if read.is_unmapped:
						pass
					else:
						counter+=1 #otherwise go into mate

						cigartup= read.cigartuples #each read take cigar string, indx position, quality, forward/reverse
						# cigarstr = read.cigarstring
						unfiltreadidx = read.get_reference_positions(full_length=True)
						unfiltread = list(read.query_sequence)
						unfiltreadqual = read.query_qualities
						isReverse = read.is_reverse #tells if forward or reverse read

						###### analyzer will populate dictionaries #######
						analyzer(isReverse,seqUpdater(cigartup,unfiltread,unfiltreadidx,unfiltreadqual))
						# createCombineDict()
				sys.stdout.write('[done]\n')
				sys.stdout.write('\tAnalyzing frequencies and binomial checks ...')
				for ntpos in sorted(CONSENSUS_DICT.iterkeys()): #looking at each nt position going down
					totalcount = sum(CONSENSUS_DICT[ntpos].values())
					posdict = {k: v / float(totalcount) for k, v in CONSENSUS_DICT[ntpos].iteritems()}
					if posdict: #If posdict empty, this means there is no coverage
						#max() will throw an error if dictionary is empty.
						maxfreq = max(posdict.iteritems(), key=operator.itemgetter(1))[1]

						if maxfreq >= (1 - args.cutoff): # indicates no "significant" minor variant
							hasMinor = False
							x = printer(VARLIST_OUT,SAMPLENAME,SEGMENT,ntpos,False,hasMinor)
							SEQLIST.append(x)
						elif totalcount <= (args.covercutoff): # indicates no "significant" minor variant
							hasMinor = False
							x = printer(VARLIST_OUT,SAMPLENAME,SEGMENT,ntpos,False,hasMinor)
							SEQLIST.append(x)

						else: #indicates a minor variant is present
							hasMinor = True
							#proceed with binomial check
							binomPass = binomCheck(ntpos)
							x = printer(VARLIST_OUT,SAMPLENAME,SEGMENT,ntpos,binomPass,hasMinor)
							SEQLIST.append(x)
					else:
						#Added new character (N) if there is NO coverage. 8/30/2017
						x = seqObject(SAMPLENAME,SEGMENT,ntpos,'N','0.00','','','',
							'0','0','0','0','0','0')
						SEQLIST.append(x)

				sys.stdout.write('[done]\n')
				sys.stdout.write('\tBuilding codons and amino acids ...')

				if len(SEQLIST) % 3 != 0:
					sys.stdout.write('[skip, not divisible by 3]\n')
					for x in SEQLIST:
						x.printer_no3(VARLIST_OUT)
				else:
					for x in SEQLIST:
						returnCodon(x,SEQLIST)
						#print x.ntpos, x.minornt
						#can insert a print such as print x.ntpos
						x.printer(VARLIST_OUT)
					sys.stdout.write('[done]\n')
				# sys.stdout.write('-= Analysis finished for segment =-')

			except ValueError:
				print 'This mapping failed for %s' % SEGMENT


print("--- %s seconds ---" % (time.time() - start_time))
