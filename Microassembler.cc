#include "Microassembler.hh"

/******************************************************************
** Microassembler.cc
**
** Tool for localized assembly of genomic regions using
** the de Bruijn paradigm to detect genetic variants
** 
**  Authors: Giuseppe Narzisi & Michael C. Schatz
**    Date: December 11, 2013
**
*******************************************************************/

void Microassembler::printConfiguration(ostream & out)
{
	out << "VERBOSE: "          << VERBOSE << endl;

	out << "bamfile: "          << BAMFILE << endl;
	out << "reffile: "          << REFFILE << endl;
	out << "readset: "          << READSET << endl;
	out << "prefix: "           << PREFIX  << endl;

	out << "min K: "            << minK << endl;
	out << "max K: "            << maxK << endl;
	out << "MAX_TIP_LEN: "      << MAX_TIP_LEN << endl;

	out << "QV_RANGE: "         << QV_RANGE << endl;
	out << "MIN_QV: "           << MIN_QV << endl;
	out << "MIN_QUAL: "         << (char) MIN_QUAL << endl;
	out << "MIN_MAP_QUAL: "     << MIN_MAP_QUAL << endl;

	out << "INCLUDE_BASTARDS: " << bvalue(INCLUDE_BASTARDS) << endl;

	out << "MIN_THREAD_READS: " << MIN_THREAD_READS << endl;
	out << "COV_THRESHOLD: "    << COV_THRESHOLD << endl;
	cerr.unsetf(ios::floatfield); // floatfield not set
	cerr.precision(5);
	out << "MIN_COV_RATIO: "    << MIN_COV_RATIO << endl;
	cerr.setf(ios::fixed,ios::floatfield);
	cerr.precision(1);
	out << "TIP_COV_THRESHOLD: "<< TIP_COV_THRESHOLD << endl;
	out << "DFS_LIMIT: "        << DFS_LIMIT << endl;
	out << "PATH_LIMIT: "       << PATH_LIMIT << endl;
	out << "MAX_INDEL_LEN: "    << MAX_INDEL_LEN << endl;
	out << "MAX_MISMATCH: "     << MAX_MISMATCH << endl;

	out << "SCAFFOLD_CONTIGS: " << bvalue(SCAFFOLD_CONTIGS) << endl;
	out << "INSERT_SIZE: "      << INSERT_SIZE << " +/- " << INSERT_STDEV << endl;

	out << "PRINT_ALL: "        << bvalue(PRINT_ALL) << endl;
	out << "PRINT_RAW: "        << bvalue(PRINT_RAW) << endl;

	out << "PRINT_DENOVO: "     << bvalue(PRINT_DENOVO) << endl;
	out << "PRINT_REFPATH: "    << bvalue(PRINT_REFPATH) << endl;
	out << "NODE_STRLEN: "      << NODE_STRLEN << endl;

	out << "QUAD_ASM: "         << bvalue(QUAD_ASM) << endl;
    out << "FASTQ_ASM: "        << bvalue(FASTQ_ASM) << endl;
}

// loadRef
//////////////////////////////////////////////////////////////

void Microassembler::createRef(const string & hdr, const string & s, int refstart, int refend) {
    Ref_t* ref = new Ref_t(minK);
    ref->hdr = hdr;
    ref->seq = s;
    ref->rawseq = s;
    ref->refchr = hdr;
    ref->refstart = refstart;
    ref->refend = refend;

    //cerr << "label:\t"    << ref->hdr << endl;
    //cerr << "refchr:\t"   << ref->refchr << endl;
    //cerr << "refstart:\t" << ref->refstart << endl;
    //cerr << "refend:\t"   << ref->refend << endl;
    
    reftable.insert(make_pair(hdr, ref));
}

// load Red Groups
//////////////////////////////////////////////////////////////

void Microassembler::loadRG(const string & filename, int member) {

	FILE * fp;
	
	switch (member) {
		case FATHER:
			fp = xfopen(PREFIX + "/father.rg", "r");
			break;
		case MOTHER:
			fp = xfopen(PREFIX + "/mother.rg", "r");
			break;
		case SELF:
			fp = xfopen(PREFIX + "/self.rg", "r");
			break;
		case SIBLING:
			fp = xfopen(PREFIX + "/sibling.rg", "r");
			break;
		default:
			fp = xfopen(PREFIX + "/" + filename, "r");
	}
	
	char rgbuffer[BUFFER_SIZE];
	string rg;

	while (fscanf(fp, "%s\n", rgbuffer) == 1) {
		//memcpy(a,s.c_str(),s.size());		
		switch (member) {
			case FATHER:
				RG_father.insert(rgbuffer);
				break;
			case MOTHER:
				RG_mother.insert(rgbuffer);
				break;
			case SELF:
				RG_self.insert(rgbuffer);
				break;
			case SIBLING:
				RG_sibling.insert(rgbuffer);
				break;
			default:
				readgroups.insert(rgbuffer);
		}
	}
	
	if(readgroups.empty()) { // insert empty symbol
		readgroups.insert("null");
	}

	xfclose(fp);
}


// processGraph
//////////////////////////////////////////////////////////////////////////

void Microassembler::processGraph(Graph_t & g, const string & refname, const string & prefix, int minkmer, int maxkmer)
{	
	if (!refname.empty())
	{
		graphCnt++;
		
		//VERBOSE = false;
		
		// skip region if no mapped reads
		if(g.countMappedReads()<=0) { return; }

		cerr << "== Processing " << graphCnt << ": " << refname 
			<< " numsequences: " << g.readid2info.size() 
			<< " mapped: " << g.countMappedReads()
			<< " bastards: " << g.countBastardReads()
			<< endl;

		cerr << "=====================================================" << endl;

		// Load the reference
		map<string, Ref_t *>::iterator ri = reftable.find(refname);
		if (ri == reftable.end())
		{
			cerr << "Can't find ref info for: " << refname << endl;
			exit(1);
		}
		Ref_t * ref = ri->second;
        
		
		bool rptInRef;
		bool rptInQry;
		bool cycleInGraph;

		// dinamic kmer mode
		for (int k=minkmer; k<=maxkmer; k++) {
			g.setK(k);
			ref->setK(k);
			
			rptInRef = false;
			rptInQry = false;
			cycleInGraph = false;

			// exit if the region has a repeat of size K
			if(isRepeat(ref->rawseq, k)) { 
				cerr << "Repeat in reference sequence for kmer " << k << endl;
				rptInRef = true;
				//return; 
				continue;
			} 

			// exit if the region has an almost perfect repeat of size K
			if(isAlmostRepeat(ref->rawseq, k, MAX_MISMATCH)) { 
				cerr << "Near-perfect repeat in reference sequence for kmer " << k << endl;
				rptInRef = true;
				//return; 
				continue;
			}
			
			//if no repeats in the reference build graph			
			g.buildgraph();
			
			double avgcov = ((double) g.totalreadbp_m) / ((double)ref->rawseq.length());			

			cerr << "reads: "   << g.readid2info.size()
                 << " reflen: "  << ref->rawseq.length()
                 << " readlen: " << g.totalreadbp_m
                 << " cov: "     << avgcov << endl;

			//printReads();
			g.printStats();
	
			string out_prefix = prefix + "."  + READSET + "." + refname;
	
			if (PRINT_ALL) { g.printDot(out_prefix + ".0.dot"); }

			// mark source and sink
			g.markRefEnds(ref);
			//g.markRefNodes();
			
			// if there is a cycle in the graph skip analysis
			if (g.hasCycle()) { g.clear(false); cycleInGraph = true; continue; }

			g.checkReadStarts();
	
			// Initial compression
			g.compress(); 
			g.printStats();
			if (PRINT_RAW) { g.printDot(out_prefix + ".1c.dot"); }

			// Remove low coverage
			g.removeLowCov();
			if (PRINT_ALL) { g.printDot(out_prefix + ".2l.dot"); }


			// Remove tips
			g.removeTips();
			if (PRINT_ALL) { g.printDot(out_prefix + ".3t.dot"); }

			// skip analysis if there is a cycle in the graph 
			if (g.hasCycle()) { g.clear(false); cycleInGraph = true; continue; }

			// skip analysis if there is a perfect or near-perfect repeat in the graph paths			
			if(g.hasRepeatsInGraphPaths()) { g.clear(false); rptInQry = true; continue; }
			
			// Thread reads
			// BUG: threding is off because creates problems if the the bubble is not covered (end-to-end) 
			// by the reads. This is particularly problematic for detecting denovo events
			//g.threadReads();
			//if (PRINT_ALL) { g.printDot(out_prefix + ".4thread.dot"); }
		
			// scaffold contigs
			if (SCAFFOLD_CONTIGS)
			{
				g.scaffoldContigs();
			}

			if (PRINT_DENOVO)
			{
				g.denovoNodes(out_prefix + ".denovo.fa", refname);
			}  

			if (PRINT_REFPATH)
			{
				g.markRefNodes();
				g.countRefPath(out_prefix + ".paths.fa", refname, false);
				//g.printFasta(prefix + "." + refname + ".nodes.fa");
			}

			if (PRINT_ALL) { g.printDot(out_prefix + ".final.dot"); }
			
			g.clear(true);
					
			break; // break loop if graph has been processed correctly
		}
		
		if(rptInRef) { cerr << " Found repeat in reference" << endl; }
		if(rptInQry) { cerr << " Found repeat in assembly" << endl; }	
		if(cycleInGraph) { cerr << " Found cycle in assembly" << endl; }	
		
		cerr << "FINISHED" << endl;
	}
}

// fastqAsm
//////////////////////////////////////////////////////////////////////////

void Microassembler::fastqAsm(Graph_t & g, const string & prefix)
{
	cerr << endl << endl;
	cerr << "fastqAsm" << endl;
	cerr << "=====================================================" << endl;

	string refname = "fastq";
	string out_prefix = prefix + "."  + READSET + "." + refname;

	g.printStats();
		
	if (PRINT_ALL) 
    {
        g.printDot(out_prefix + ".0c.dot");
        g.printFasta(out_prefix + ".0c.fa");
    }

	if (PRINT_ALL) 
    { 
        g.printDot(out_prefix + ".0.dot"); 
        g.printFasta(out_prefix + ".0.fa"); 
    }
    
	// Initial compression
	g.compress(); 
	g.printStats();

	if (PRINT_RAW) 
    {
        g.printDot(out_prefix + ".1c.dot");
        g.printFasta(out_prefix + ".1c.fa");
    }

	// Remove low coverage
	g.removeLowCov();
	if (PRINT_ALL) 
    { 
        g.printDot(out_prefix + ".2l.dot"); 
        g.printFasta(out_prefix + ".2l.fa"); 
    }

	// Remove tips
	g.removeTips();
	if (PRINT_ALL) 
    { 
        g.printDot(out_prefix + ".3t.dot"); 
        g.printFasta(out_prefix + ".3t.fa"); 
    }

	// Thread reads
	g.threadReads();
	if (PRINT_ALL) 
    { 
        g.printDot(out_prefix + ".4thread.dot"); 
        g.printFasta(out_prefix + ".4thread.fa"); 
    }

	// scaffold contigs
	if (SCAFFOLD_CONTIGS)
	{
		g.scaffoldContigs();
	}

	g.printDot(out_prefix + ".final.dot");
	g.printFasta(out_prefix + ".final.fa");
	g.printPairs(out_prefix + ".pairs.fa");
	g.clear(true);

	cerr << endl;
}

int Microassembler::run(int argc, char** argv)
{
	string USAGE = "Usage: Microassembler [options] -b bamfile -f ref.fa -r chr:start-end\n";

	if (argc == 1)
	{
		cerr << USAGE;
		exit(0);
	}

	cerr.setf(ios::fixed,ios::floatfield);
	cerr.precision(1);

	stringstream helptext;
	helptext << USAGE <<
		"\n"
		"   -f <reffile>  : multifasta file to reference regions\n"
		"   -b <bamfile>  : file of mapped reads\n"
        "   -r <region>   : target region in chr:start-end format e.g. chrX:0-10\n"
		"   -s <name>     : label for reads (default: " << READSET << ")\n"
		"   -p <prefix>   : use prefix (default: region)\n"
		"\n"
		"   -k <kmersize> : min kmersize (default: " << minK << ")\n"
		"   -K <kmersize> : max kmersize (default: " << maxK << ")\n"
		"   -q <qv>       : trim bases below qv at 5' and 3' (default: " << MIN_QV << ")\n"
		"   -Q <char>     : quality value range (default: " << (char) QV_RANGE << ")\n"
		"   -C <mq>       : minimum read mapping quality in Phred-scale (default: " << MIN_MAP_QUAL << ")\n"
		"   -l <tiplen>   : max tip length (default: " << MAX_TIP_LEN << ")\n"
		"   -t <reads>    : min number of reads to thread (default: " << MIN_THREAD_READS << ")\n"
		"   -c <cov>      : coverage threshold (default: " << COV_THRESHOLD << ")\n"
		"   -x <covratio> : minimum coverage ratio (default: " << MIN_COV_RATIO << ")\n"
		"   -d <tipcov>   : tip coverage threshold (default: " << TIP_COV_THRESHOLD << ")\n"
		"   -F <dfs>      : limit dfs search space (default: " << DFS_LIMIT << ")\n"
		"   -P <maxp>     : limit on number of paths to report (default: " << PATH_LIMIT << ")\n"
		"   -T <maxindel> : limit on size of detectable indel (default: " << MAX_INDEL_LEN << ")\n"
		"   -M <max-mismatch> : max number of mismatches for near-perfect repeats (default: " << MAX_MISMATCH << ")\n"
		"   -B            : include bastards\n"
		"\n"
		"   -D            : print de novo mutations (read map must be sorted and tagged with id)\n"
		"   -R            : print reference paths\n"
		"   -A            : print graph after every stage\n"
		"   -I            : don't print initial graph\n"
		"   -L <len>      : length of sequence to display at graph node (default: " << NODE_STRLEN << ")\n"
		"\n"
		"   -E            : fastq assembly, map file is in fq (experimental)\n"
		"   -v            : be verbose\n"
		"\n";

	bool errflg = false;
	int ch;

	optarg = NULL;

	while (!errflg && ((ch = getopt (argc, argv, "b:r:g:s:k:K:l:t:c:d:x:BDRAIhSL:T:M:vF:q:C:f:Q:P:p:E")) != EOF))
	{
		switch (ch)
		{
        case 'b': BAMFILE          = optarg;       break; 
        case 'f': REFFILE          = optarg;       break;
        case 'r': REGION           = optarg;       break;
        case 's': READSET          = optarg;       break;
        case 'p': PREFIX           = optarg;       break;

        case 'k': minK             = atoi(optarg); break;
        case 'K': maxK             = atoi(optarg); break;
        case 'l': MAX_TIP_LEN      = atoi(optarg); break;
        case 't': MIN_THREAD_READS = atoi(optarg); break;
        case 'c': COV_THRESHOLD    = atoi(optarg); break;
        case 'x': MIN_COV_RATIO    = atof(optarg); break;
        case 'd': TIP_COV_THRESHOLD= atoi(optarg); break;

        case 'B': INCLUDE_BASTARDS = 1;            break;

        case 'v': VERBOSE          = 1;            break;
        case 'D': PRINT_DENOVO     = 1;            break;
        case 'R': PRINT_REFPATH    = 1;            break;
        case 'A': PRINT_ALL        = 1;            break;
        case 'I': PRINT_RAW        = 0;            break;
        case 'L': NODE_STRLEN      = atoi(optarg); break;
        case 'F': DFS_LIMIT        = atoi(optarg); break;
        case 'P': PATH_LIMIT       = atoi(optarg); break;
        case 'T': MAX_INDEL_LEN    = atoi(optarg); break;
        case 'M': MAX_MISMATCH     = atoi(optarg); break;

        case 'S': QUAD_ASM         = 1;			   break;
        case 'E': FASTQ_ASM       = 1;            break;
		  
        case 'q': MIN_QV           = atoi(optarg); break;
        case 'C': MIN_MAP_QUAL     = atoi(optarg); break;
        case 'Q': QV_RANGE         = *optarg;      break;

        case 'h': errflg = 1;                      break;

        case '?':
			fprintf (stderr, "Unrecognized option -%c\n", optopt);

        default:
			errflg = true;
		}

		if (errflg)
		{
			cout << helptext.str();
			exit (EXIT_FAILURE);
		}
	}

    if (REGION.empty()) { cerr << "ERROR: A region is required" << endl; errflg++; }
    if (PREFIX.empty()) { PREFIX = REGION; }

    if (BAMFILE.empty()) { cerr << "ERROR: Must provide a bamfile (-b)" << endl; errflg++; }
    if (REFFILE.empty()) {
        cerr << "ERROR: Must provide a reffile (-f)" << endl; errflg++;
    } else {
        reference.open(REFFILE);
    }

	if (errflg) { exit(EXIT_FAILURE); }

	printConfiguration(cerr);

	// Process the reads
	FILE * fp = NULL;
	BamReader reader;
	//BamWriter writer;
	SamHeader header;
	RefVector references;
	if(!BAMFILE.empty()) { // attempt to open our BamMultiReader
		if ( !reader.Open(BAMFILE) ) {
			cerr << "Could not open input BAM files." << endl;
			return -1;
		}
		// retrieve 'metadata' from BAM files, these are required by BamWriter
		header = reader.GetHeader();
		references = reader.GetReferenceData();
		reader.LocateIndex(); // locate and load BAM index file

		/*
          string outputBamFilename = PREFIX + "/outputBam.bam"; 
          // attempt to open our BamWriter
          if ( !writer.Open(outputBamFilename, header, references) ) {
          cerr << "Could not open output BAM file" << endl;
          return -1;
          }
		*/
		
		/*
          vector<RefData>::iterator it;
          cout << "refvector contains:" << endl;
          for ( it=references.begin() ; it < references.end(); it++ ) {
          cout << (*it).RefName << " " << (*it).RefLength << endl;
          }
		*/

		//load the read group information
        readgroups.insert("null");
		
		//set<string>::iterator prova;
		//for ( prova=readgroups.begin() ; prova != readgroups.end(); prova++ ) {
		//	cout << (*prova) <<  endl;
		//}
	}

	Graph_t g;

	//set configuration parameters
	g.setK(minK);
	g.setVerbose(VERBOSE);
	g.setMinQual(MIN_QUAL);
	g.setIncludeBastards(INCLUDE_BASTARDS);
	g.setBufferSize(BUFFER_SIZE);
	g.setDFSLimit(DFS_LIMIT);
	g.setPathLimit(PATH_LIMIT);
	g.setCovThreshold(COV_THRESHOLD);
	g.setMinCovRatio(MIN_COV_RATIO);
	g.setTipCovThreshold(TIP_COV_THRESHOLD);
	g.setPrintDotReads(PRINT_DOT_READS);
	g.setNodeStrlen(NODE_STRLEN);
	g.setMaxTipLength(MAX_TIP_LEN);
	g.setMaxIndelLen(MAX_INDEL_LEN);
	g.setMinThreadReads(MIN_THREAD_READS);
	g.setScaffoldContigs(SCAFFOLD_CONTIGS);
	g.setInsertSize(INSERT_SIZE);
	g.setInsertStdev(INSERT_STDEV);
	g.setMaxMismatch(MAX_MISMATCH);

	string graphref = "";

	char code [BUFFER_SIZE];
	char chr  [BUFFER_SIZE];
	int  refstart;
	int  refend;
	char readname [BUFFER_SIZE];
	char seq1 [BUFFER_SIZE];
	char qv1  [BUFFER_SIZE];
	char seq2 [BUFFER_SIZE];
	char qv2  [BUFFER_SIZE];

	char refbuf[BUFFER_SIZE];

	int paircnt = 0;
	int graphcnt = 0;
	int readcnt = 0;

	int idx;

	if (!BAMFILE.empty()) {
		
		// for each reference location
		BamRegion region;
        string chr; int start, end;
        parseRegion(REGION, chr, start, end);
		region.LeftRefID = reader.GetReferenceID(chr); // atoi((refinfo->refchr).c_str());
        region.RightRefID = reader.GetReferenceID(chr); // atoi((refinfo->refchr).c_str());
        region.LeftPosition = start;
        region.RightPosition = end;
        cout << "region = " << chr << ":" << start << "-" << end << endl; 

        string refseq = reference.getSubSequence(chr, start, end-start);
        createRef(chr, refseq, start, end);

        // continue if the region has only Ns or prefect repeat of size maxK
        if(isNseq(refseq)) { cerr << "region is only N" << endl; return 0; } 
        if(isRepeat(refseq, maxK)) { cerr << "region is repeat of maxK=" <<maxK << endl; return 0; } 

	

        bool jump = reader.SetRegion(region);
        if(!jump) {
            cout << "Error: not able to jump successfully to the region's left boundary" << endl;
            return -1;
        }

        // iterate through all alignments
        //int num_PCR_duplicates = 0;
        BamAlignment al;
        string rg = "";
        int num_unmapped = 0;
        while ( reader.GetNextAlignment(al) ) { // get next alignment and populate the alignment's string data fields
				
            // skip alignments outside region
            //int alstart = al.Position;
            //int alend = al.GetEndPosition();
            //if( (alstart < region.LeftPosition) || (alend > region.RightPosition) ) { continue; }
				
            if ( (al.MapQuality >= MIN_MAP_QUAL) && !al.IsDuplicate() ) { // only keeping ones with high map quality and skip PCR duplicates
										
                int mate = 0;
                if(al.IsFirstMate()) { mate = 1; }
                if(al.IsSecondMate()) { mate = 2; }
				
                al.GetTag("RG", rg); // get the read group information for the read
                if(rg.empty()) { rg = "null"; }
					
                if ( (readgroups.find("null") != readgroups.end())  || (readgroups.find(rg) != readgroups.end()) ) { // select reads in the read group RG
							
                    //writer.SaveAlignment(al); // save alignment to output bam file
							
                    if (mate>1) { // mated pair
                        if( !(al.IsMapped()) ) { // unmapped read
                            g.addpaired(READSET, al.Name, al.QueryBases, al.Qualities, mate, Graph_t::CODE_BASTARD);
                            num_unmapped++; 
                        }
                        else { // mapped reads
                            g.addpaired(READSET, al.Name, al.QueryBases, al.Qualities, mate, Graph_t::CODE_MAPPED);								
                        }
                    }
                    else { // unpaired
                        g.addUnpaired(READSET, al.Name, al.QueryBases, al.Qualities, Graph_t::CODE_MAPPED);	
                    }
                    //cout << al.Name << endl;
                    readcnt++;

                    //void addMates(ReadId_t r1, ReadId_t r2)
                    //{
                    //g.readid2info[r1].mateid_m = r2;
                    //g.readid2info[r2].mateid_m = r1;
                    //}
                }

            }
            //else{ num_PCR_duplicates++; }
        }
        // close the reader & writer
        //cout << "Number of PCR duplicates: " << num_PCR_duplicates << endl;	
        //cout << "Number of unmapped reads: " << num_unmapped << endl;	
        processGraph(g, chr, PREFIX, minK, maxK);

        reader.Close();
        //writer.Close();
	}
	else
	{
		// single set of reads against a reference

		while (fscanf(fp, "%s\t%d\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s",
				code, &idx, chr, &refstart, &refend, readname, seq1, qv1, seq2, qv2) == 10)
		{
			paircnt++;
			readcnt+=2;

			// printf("c:%s s:%d e:%d n:%s s:%s q:%s s:%s q:%s\n", 
			//         chr, refstart, refend, readname, seq1, qv1, seq2, qv2);
			//
			snprintf(refbuf, BUFFER_SIZE, "%s:%d-%d", chr, refstart, refend);

			string refstr = refbuf;

			if (refstr != graphref)
			{
				processGraph(g, graphref, PREFIX, minK, maxK);
				graphref = refstr;
				graphcnt++;
			}
			g.addPair(READSET, readname, seq1, qv1, seq2, qv2, code[0]);
		}
		processGraph(g, graphref, PREFIX, minK, maxK);
	}

	cerr << "=======" << endl;

	cerr << "total reads: " << readcnt << " pairs: " << paircnt << " total graphs: " << graphcnt << " ref sequences: " << reftable.size() <<  endl;

	return 0;
}

// main
//////////////////////////////////////////////////////////////////////////


int main(int argc, char** argv)
{
	try {
		Microassembler* assembler = new Microassembler();
		assembler->run(argc, argv);
	}
	catch (int e) {
		cout << "An exception occurred. Exception Nr. " << e << endl;
	} 
	//catch(std::out_of_range& e) {
   	//	cerr << e.what( ) << endl;
 	//}
	//catch (...) { cout << "default exception"; }
}
