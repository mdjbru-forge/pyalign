### * Description

# Entry points for the command line scripts

DESCRIPTION_CONSENSUS = ("Determine consensus sequences from alignments in "
 "fasta format. The consensus sequences are written to stdout. Conservation "
 "profiles for each alignment file and average conservation score for each "
 "alignment can also be produced.")

### * Wishlist

# pyalign align clusters/* alignments # ungap sequences and run mafft
# pyalign consensus alignments/* consensus.fa
# other commands to produce conservation profiles?
# pyalign split alignments/* --outDir alignmentsSplitted --threshold 0.4
# pyalign align alignmentsSplitted/* # with no output dir, erase the original files

# To do
# pyquickmcl input.fa # perform a quick clustering with mcl to split a suspicious alignment

### * Set up

### ** Import

import sys
import os
import argparse
import hashlib
import collections
import shutil
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import pygenes as pygenes
import pyalign as pyalign

### * Parser and subparsers

def makeParser() :
    """Prepare the parser

    Returns:
        ArgumentParser: An argument parser

    """
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help = "")
    # Align
    sp_align = subparsers.add_parser("align",
                                     help = "Ungap and align sequences in fasta "
                                     "file(s) using mafft")
    sp_align.add_argument("input", metavar = "FASTA_FILE",
                          type = str, nargs = "+",
                          help = "Fasta file")
    sp_align.add_argument("-o", "--outDir", metavar = "DIR", type = str,
                          help = "Output directory (default: overwrite the "
                          "input file)")
    sp_align.set_defaults(action = "align")
    # Consensus
    sp_consensus = subparsers.add_parser("consensus",
                                         description = DESCRIPTION_CONSENSUS,
                                         help = "Generate consensus sequences "
                                         "and conservation information from "
                                         "fasta alignments")
    sp_consensus.add_argument("input", metavar = "FASTA_FILE", type = str,
                              nargs = "+",
                              help = "Alignment in fasta format")
    sp_consensus.add_argument("-p", "--profiles", metavar = "FILE", type = str,
                              help = "Conservation profiles output file")
    sp_consensus.add_argument("-c", "--conservation", metavar = "FILE", type = str,
                              help = "Average conservation output file")
    sp_consensus.set_defaults(action = "consensus")
    # # Split alignments (mcl)
    # sp_splitMcl = subparsers.add_parser("splitMcl",
    #                                  help = "Split alignments containing more "
    #                                  "than one family using mcl")
    # sp_splitMcl.add_argument("input", metavar = "FASTA_FILE", type = str,
    #                       nargs = "+",
    #                       help = "Alignment in fasta format")
    # sp_splitMcl.add_argument("conservation", metavar = "CONSERVATION", type = str,
    #                       help = "Conservation for each alignment file, output "
    #                       "from the -c option of pyalign consensus")
    # sp_splitMcl.add_argument("-t", "--threshold", metavar = "FLOAT",
    #                       type = float,
    #                       help = "Conservation threshold below which "
    #                       "an alignment is analyzed again (default: 0.8)",
    #                       default = 0.8)
    # sp_splitMcl.add_argument("-I", "--inflation", metavar = "FLOAT",
    #                       type = float,
    #                       help = "Inflation for the mcl algorithm (default: "
    #                       "10)",
    #                       default = 10)
    # sp_splitMcl.add_argument("-n", "--nThreads", metavar = "INT",
    #                       type = int,
    #                       help = "Number of cores to use for blastp "
    #                       "(default: 1)")
    # sp_splitMcl.add_argument("-o", "--outDir", metavar = "DIR", type = str,
    #                       default = ".",
    #                       help = "Output directory (default: current "
    #                       "directory)")
    # sp_splitMcl.add_argument("-k", "--keep", action = "store_true",
    #                       help = "Keep the original alignment when it is split "
    #                       "(default: remove the original file)")
    # sp_splitMcl.set_defaults(action = "splitMcl")
    # Split alignments (hierarchical clustering)
    sp_splitHclust = subparsers.add_parser("splitHclust",
                                          help = "Split alignments containing "
                                           "more than one family using "
                                           "hierarchical clustering.")
    sp_splitHclust.add_argument("input", metavar = "FASTA_FILE", type = str,
                                nargs = "+",
                                help = "Input fasta file. Sequences should have "
                                "the same length, gapped alignments are accepted.")
    sp_splitHclust.add_argument("conservation", metavar = "CONSERVATION", type = str,
                                help = "Conservation for each alignment file, output "
                                "from the -c option of pyalign consensus") 
    sp_splitHclust.add_argument("-t", "--threshold", metavar = "FLOAT",
                                type = float,
                                help = "Conservation threshold below which "
                                "an alignment is analyzed again (default: 0.8)",
                                default = 0.8)
    sp_splitHclust.add_argument("-k", "--keep", action = "store_true",
                                help = "Keep the original alignment when it is split "
                                "(default: remove the original file)")
    sp_splitHclust.add_argument("-d", "--dissim", metavar = "FLOAT", type = float,
                               default = 0.05,
                               help = "Maximum dissimilarity for merging and splitting"
                               "(between 0 and 1) (default: 0.05)")
    sp_splitHclust.add_argument("-o", "--outDir", metavar = "DIR", type = str,
                               help = "Output directory for splitted fasta "
                               "files (default: current directory)")
    sp_splitHclust.set_defaults(action = "splitHclust")
    # Validate alignments
    sp_validate = subparsers.add_parser("validate",
                                        help = "Validate alignments based on "
                                        "their conservation score and the "
                                        "number of sequences they contain")
    sp_validate.add_argument("input", metavar = "FASTA_ALN", type = str,
                             nargs = "+",
                             help = "Input alignments in fasta format. If "
                             "--outDir is not specified, input files which "
                             "do not pass the validation criteria will be "
                             "deleted.")
    sp_validate.add_argument("-s", "--seqcons", metavar = "FLOAT",
                             type = float, default = 0,
                             help = "Minimum conservation score for a sequence "
                             "to be kept in the alignment (default: 0)")
    sp_validate.add_argument("-n", "--n_seqs", metavar = "INT", type = int,
                             default = 1,
                             help = "Minimum number of sequences for "
                             "validation (after removing sequences with "
                             "--seqcons) (default: 1)")
    sp_validate.add_argument("-c", "--conservation", metavar = "FLOAT",
                             type = float, default = 0,
                             help = "Minimum conservation score for validation "
                             "(after removing sequences using --seqcons) "
                             "(default: 0)")
    sp_validate.add_argument("-o", "--outDir", metavar = "DIR", type = str,
                             help = "Output directory for splitted fasta "
                             "files (default: current directory)")
    sp_validate.set_defaults(action = "validate")
    # phaseNt
    sp_phaseNt = subparsers.add_parser("phaseNt",
                                        help = "Convert protein alignments to "
                                        "detailled nucleotide alignments")
    sp_phaseNt.add_argument("alnFiles", metavar = "FASTA_FILE", type = str,
                             nargs = "+",
                             help = "Alignment files in fasta format")
    sp_phaseNt.add_argument("geneTable", metavar = "GENE_TABLE", type = str,
                             help = "Gene table")
    sp_phaseNt.add_argument("-o", "--outDir", metavar = "DIR", type = str,
                             help = "Output directory for splitted fasta "
                             "files (default: current directory)")
    sp_phaseNt.set_defaults(action = "phaseNt")
    # Ungap alignments
    sp_ungap = subparsers.add_parser("ungap",
                                     help = "Ungap alignments and fasta files. "
                                     "To clean an alignment, one approach is "
                                     "first to remove gappy positions (-p) and "
                                     "then gappy sequences (-s)",
                                     description = "Ungap alignments and fasta files. "
                                     "To clean an alignment, one approach is "
                                     "first to remove gappy positions (-p) and "
                                     "then gappy sequences (-s)")
    sp_ungap.add_argument("input", metavar = "FASTA_FILE", type = str,
                          nargs = "+",
                          help = "Input sequences or alignments in fasta "
                          "format. If --outDir is not specified, input files "
                          "will be overwritten by their ungapped version.")
    sp_ungap.add_argument("-a", "--all", action = "store_true",
                          help = "Remove all gap characters, does not take "
                          "into account any alignment information")
    sp_ungap.add_argument("-p", "--position", metavar = "PROPORTION",
                          type = float,
                          help = "Remove gappy positions in alignment if gap "
                          "proportion is at least equal to PROPORTION")
    sp_ungap.add_argument("-s", "--seq", metavar = "PROPORTION",
                          type = float,
                          help = "Remove sequences with at least one gap "
                          "whose proportion is less than PROPORTION")
    sp_ungap.add_argument("-o", "--outDir", metavar = "DIR", type = str,
                          help = "Output directory for splitted fasta "
                          "files (default: current directory)")
    sp_ungap.add_argument("--syncNtDir", metavar = "DIR", type = str,
                          help = "Synchronize the ungapped alignments with "
                          "detailled nucleotide alignments in another "
                          "directory. The nucleotide alignment files will "
                          "be modified in place, not copied.")
    sp_ungap.set_defaults(action = "ungap")
    # Check origin
    sp_origin = subparsers.add_parser("origin",
                                      help = "Check the origin of the "
                                      "sequences in an alignment, based on a "
                                      "gene table file. For now only return "
                                      "the files for which one or several "
                                      "records have multiple entries")
    sp_origin.add_argument("geneTable", metavar = "GENE_TABLE", type = str,
                           help = "Gene table file")
    sp_origin.add_argument("alnFiles", metavar = "FASTA_FILE", type = str,
                           nargs = "+",
                           help = "Alignment files in fasta format")
    sp_origin.set_defaults(action = "origin")
    # Compile gene information
    sp_compile = subparsers.add_parser("compile",
                                       help = "Compile gene information for "
                                       "each cluster")
    sp_compile.add_argument("geneTable", metavar = "GENE_TABLE", type = str,
                            help = "Gene table file")
    sp_compile.add_argument("alnFiles", metavar = "FASTA_FILE", type = str,
                            nargs = "+",
                            help = "Alignment files in fasta format")
    sp_compile.add_argument("-i", "--info", metavar = "FIELD", type = str,
                            nargs = "+",
                            help = "Fields to extract")
    sp_compile.set_defaults(action = "compile")
    # Scan alignments to get composition information
    sp_scan = subparsers.add_parser("scan",
                                    help = "Scan alignment positions to get "
                                    "their composition")
    sp_scan.add_argument("alnFiles", metavar = "FASTA_FILE", type = str,
                         nargs = "+",
                         help = "Alignment files in fasta format")
    sp_scan.set_defaults(action = "scan")
    # SNP calling
    sp_snp = subparsers.add_parser("callSNP",
                                   help = "Call SNPs from detailled nucleotide "
                                   "alignments")
    sp_snp.add_argument("gene2record", metavar = "GENE2RECORD", type = str,
                        help = "Tabular file containing the mapping between "
                        "gene id and record id (e.g. produced by pygenes "
                        "extract)")
    sp_snp.add_argument("alnFiles", metavar = "ALNFILE", type = str,
                        nargs = "+",
                        help = "Detailled nucleotide alignment file(s)")
    sp_snp.set_defaults(action = "callSNP")
    # SNPtable
    sp_SNPtable = subparsers.add_parser("SNPtable",
                                        help = "Filter a SNP table")
    sp_SNPtable.add_argument("snpTable", metavar = "SNPTABLE", type = str,
                             help = "SNP table")
    sp_SNPtable.add_argument("-s", "--strand", action = "store_true",
                             help = "Drop the strand columns")
    sp_SNPtable.add_argument("-m", "--multipleEntries", action = "store_true",
                             help = "Drop alignments for which multiple "
                             "entries of a record were present")
    sp_SNPtable.set_defaults(action = "SNPtable")
    # Return
    return parser

### * Mains

### ** Main entry point (dispatch)

def main(args = None, stdout = None, stderr = None) :
    """Main entry point

    Args:
        args (namespace): Namespace with script arguments, parse the command 
          line arguments if None
        stdout (file): Writable stdout stream (if None, use `sys.stdout`)
        stderr (file): Writable stderr stream (if None, use `sys.stderr`)

    """
    if args is None :
        parser = makeParser()
        args = parser.parse_args()
    if stdout is None :
        stdout = sys.stdout
    if stderr is None :
        stderr = sys.stderr
    dispatch = dict()
    dispatch["align"] = main_align
    dispatch["consensus"] = main_consensus
    dispatch["splitHclust"] = main_splitHclust
    dispatch["validate"] = main_validate
    dispatch["phaseNt"] = main_phaseNt
    dispatch["ungap"] = main_ungap
    dispatch["origin"] = main_origin
    dispatch["compile"] = main_compile
    dispatch["scan"] = main_scan
    dispatch["callSNP"] = main_callSNP
    dispatch["SNPtable"] = main_SNPtable
    dispatch[args.action](args, stdout, stderr)

### ** Main align

def main_align(args, stdout, stderr) :
    if args.outDir is None :
        args.outDir = "."
    for fastaFile in args.input :
        pyalign.ungapFastaFile(fastaFile, fastaFile + ".ungap.tmp")
        pyalign.runMafft(fastaFile + ".ungap.tmp", fastaFile + ".tmp")
        out = os.path.join(args.outDir, os.path.basename(fastaFile))
        os.rename(fastaFile + ".tmp", out)
        os.remove(fastaFile + ".ungap.tmp")

### ** Main consensus

def main_consensus(args, stdout, stderr) :
    stderr.write("Note that positions with only gaps are removed before "
                 "calculations\n")
    consensus = dict()
    profiles = dict()
    for fastaFile in args.input :
        try :
            aln = pyalign.AlignIO.read(fastaFile, "fasta")
            aln = pyalign.ungapAln(aln)
            k = os.path.basename(fastaFile)
            assert not k in consensus
            consensus[k] = pyalign.makeConsensus(aln)
            stdout.write(">" + k + "\n")
            stdout.write(consensus[k] + "\n")
            profiles[k] = pyalign.conservationProfile(aln)
        except ValueError :
            msg = "Problem with " + fastaFile + "\n"
            stderr.write(msg)
    keys = list(consensus.keys())
    if args.profiles is not None :
        with open(args.profiles, "w") as fo :
            for k in keys :
                fo.write("\t".join([k] + [str(x) for x in profiles[k]]) + "\n")
    if args.conservation is not None :
        def mean(x) :
            return sum(x) * 1. / len(x)
        with open(args.conservation, "w") as fo :
            for k in keys :
                fo.write(k + "\t" + str(mean(profiles[k])) + "\n")

### ** Main splitMcl

# def main_splitMcl(args, stdout, stderr) :
#     # Load conservation file and determine files to realign
#     alnToSplit = set([])
#     with open(args.conservation, "r") as fi :
#         for l in fi :
#             if l.strip() != "" :
#                 fasta, cons = l.strip().split("\t")
#                 if float(cons) < args.threshold :
#                     alnToSplit.add(fasta)
#     # Go through the alignment files
#     for fastaFile in args.input :
#         if os.path.basename(fastaFile) in alnToSplit :
#             output = os.path.join(args.outDir, fastaFile)
#             pyalign.splitAlignment(fastaFile, args.outDir, args.inflation,
#                                    args.nThreads, not args.keep)

### ** Main splitHclust

def main_splitHclust(args, stdout, stderr) :
    # Load conservation file and determine files to realign
    alnToSplit = set([])
    with open(args.conservation, "r") as fi :
        for l in fi :
            if l.strip() != "" :
                fasta, cons = l.strip().split("\t")
                if float(cons) < args.threshold :
                    alnToSplit.add(fasta)
    # Go through the input files
    if args.outDir is None :
        args.outDir = "."
    for fastaFile in args.input :
        if os.path.basename(fastaFile) in alnToSplit :
            # Build mapping from peptide sequences to sequence names
            seqParser = SeqIO.parse(fastaFile, "fasta")
            seqRaw = [x for x in seqParser]
            seqs = dict()
            [seqs.update({x.description : str(x.seq)}) for x in seqRaw]
            pep2seqNames = collections.defaultdict(lambda : [])
            [pep2seqNames[v].append(k) for (k, v) in seqs.iteritems()]
            # Produce merged sequences
            uniqueSeqs = list(set(seqs.values()))
            mergedSeqs = pygenes.mergeSequences(uniqueSeqs,
                                                maxDistance = args.dissim)
            # Build mapping from merged sequences to original peptide sequences
            merged2pep = collections.defaultdict(lambda : [])
            [merged2pep[v].append(k) for (k, v) in mergedSeqs.iteritems()]
            # Output
            for (i, v) in enumerate(merged2pep.values()) :
                outFile = os.path.join(args.outDir,
                                       (os.path.basename(fastaFile) + ".split" +
                                        str(i) + ".fa"))
                with open(outFile, "w") as fo :
                    for originalPep in v :
                        for seqName in pep2seqNames[originalPep] :
                            fo.write(">" + seqName + "\n")
                            fo.write(originalPep + "\n")
            if args.outDir == "." and not args.keep :
                os.remove(fastaFile)
        else :
            if args.outDir != "." :
                shutil.copy(fastaFile,
                            os.path.join(args.outDir,
                                         os.path.basename(fastaFile)))

### ** Main phaseNt

def main_phaseNtOld(args, stdout, stderr) :
    if args.outDir is None :
        args.outDir = "."
    # Load gene table
    geneTable = pygenes.GeneTable()
    geneTable.loadTable(args.geneTable)
    # Go through the alignments
    for alnFile in args.alnFiles :
        stderr.write("Processing alignment " + os.path.basename(alnFile) +
                     "\n")
        try :
            aln = AlignIO.read(alnFile, "fasta")
            loaded = True
        except ValueError :
            loaded = False
        if loaded :
            alnNt = pyalign.phaseNtDetailled(aln, geneTable)
            outFile = os.path.join(args.outDir, os.path.basename(alnFile) + ".alnNt")
            alnNt.writeCompactFile(outFile)

def main_phaseNt(args, stdout, stderr) :
    if args.outDir is None :
        args.outDir = "."
    # Load gene table
    geneTable = pygenes.GeneTable()
    geneTable.loadTable(args.geneTable)
    # Go through the alignments
    for alnFile in args.alnFiles :
        stderr.write("Processing alignment " + os.path.basename(alnFile) +
                     "\n")
        outFile = os.path.join(args.outDir, os.path.basename(alnFile) + ".alnNt")
        pyalign.phaseNtDetailledFast(alnFile, geneTable, outFile)

### ** Main ungap

def main_ungap(args, stdout, stderr) :
    if args.all :
        assert args.position is None and args.seq is None and args.syncNtDir is None
        for inputFile in args.input :
            stderr.write("Processing file " + inputFile + "\n")
            outFile = os.path.join(args.outDir, inputFile)
            pyalign.ungapFastaFile(inputFile, outFile)
    if args.position is not None or args.seq is not None :
        if args.position is None :
            args.position = 1
        if args.seq is None :
            args.seq = 0
        for inputFile in args.input :
            stderr.write("Processing file " + inputFile + "\n")
            if args.outDir is None :
                outFile = inputFile
            else :
                outFile = os.path.join(args.outDir, os.path.basename(inputFile))
            try :
                inputAln = AlignIO.read(inputFile, "fasta")
                loaded = True
            except ValueError :
                loaded = False
                if args.outDir == "." :
                    os.remove(outFile)
            if loaded :
                ungapResult = pyalign.ungapAln(inputAln, args.position, args.seq)
                outputAln = ungapResult["ungappedAln"] 
                with open(outFile, "w") as fo :
                    for seq in outputAln :
                        fo.write(">" + seq.description + "\n")
                        fo.write(str(seq.seq) + "\n")
                if args.syncNtDir is not None :
                    ntFile = os.path.join(args.syncNtDir,
                                          os.path.basename(inputFile) + ".alnNt")
                    pyalign.ungapNtAlnFile(ntFile, ungapResult["removedSeq"],
                                           ungapResult["removedPos"], ntFile)
                        
### ** Main validate

def main_validate(args, stdout, stderr) :
    def mean(x) :
        return sum(x) * 1. / len(x)
    if args.outDir is None :
        args.outDir = "."
    for inputFile in args.input :
        try :
            stderr.write("Processing file " + inputFile + "\n")
            aln = AlignIO.read(inputFile, "fasta")
            seqScores = pyalign.sequenceConservation(aln)
            seqsKept = [seq for (seq, score) in zip(aln, seqScores) if score >= args.seqcons]
            cleanAln = MultipleSeqAlignment(seqsKept)
            alnCons = mean(pyalign.conservationProfile(cleanAln))
            outFile = os.path.join(args.outDir, inputFile)
            if (len(cleanAln) >= args.n_seqs) and (alnCons >= args.conservation) :
                with open(outFile, "w") as fo :
                    for seq in cleanAln :
                        fo.write(">" + seq.description + "\n")
                        fo.write(str(seq.seq) + "\n")
            else :
                if args.outDir == "." :
                    os.remove(inputFile)
        except :
            stderr.write("Problem with " + inputFile + "\n")
                
### ** Main origin

def main_origin(args, stdout, stderr) :
    stderr.write("Loading gene table\n")
    geneTable = pygenes.GeneTable()
    geneTable.loadTable(args.geneTable)
    for inputFile in args.alnFiles :
        aln = AlignIO.read(inputFile, "fasta")
        origins = collections.defaultdict(lambda : [])
        for seq in aln :
            origins[geneTable.geneId(seq.description).recordId].append(seq.description)
        multipleOrigins = [(x,y) for (x,y) in origins.items() if len(y) > 1]
        for (x,y) in multipleOrigins :
            stdout.write(inputFile + "\t" + str(x) + "\t" + str(len(y)) + "\t" +
                         ";".join(y) + "\n")

### ** Main compile

def main_compile(args, stdout, stderr) :
    stderr.write("Loading gene table\n")
    geneTable = pygenes.GeneTable()
    geneTable.loadTable(args.geneTable)
    for inputFile in args.alnFiles :
        try :
            aln = AlignIO.read(inputFile, "fasta")
            loaded = True
        except ValueError :
            loaded = False
        if loaded :
            for seq in aln :
                gene = geneTable.geneId(seq.description)._asdict()
                out = [gene[x] for x in args.info]
                stdout.write("\t".join([os.path.basename(inputFile)] + out) + "\n")

### ** Main scan

def main_scan(args, stdout, stderr) :
    characters = "-ACDEFGHIKLMNPQRSTUVWXY"
    stdout.write("cluster" + "\t" + "\t".join(characters) + "\n")
    n = str(len(args.alnFiles))
    for (i, inputFile) in enumerate(args.alnFiles) :
        stderr.write("Processing file " + str(i) + " out of " + n + "\n")
        try :
            aln = AlignIO.read(inputFile, "fasta")
            loaded = True
        except ValueError :
            loaded = False
        if loaded :
            for i in xrange(aln.get_alignment_length()) :
                compo = collections.Counter(aln[:, i])
                out = collections.defaultdict(lambda : 0)
                out.update(compo)
                stdout.write("\t".join([os.path.basename(inputFile)] +
                                       [str(out[x]) for x in characters]) + "\n")

### ** Main call SNP

def main_callSNP(args, stdout, stderr) :
    # Build the gene to record mapping
    gene2recordMapping = dict()
    with open(args.gene2record, "r") as fi :
        for line in fi :
            e = line.strip().split("\t")
            gene2recordMapping[e[0]]= e[1]
    records = set(gene2recordMapping.values())
    if False :
        # First pass through the alignment files to get the record names
        # Not needed if we assume all the records are in gene2recordMapping
        # already
        records = set([])
        for f in args.alnFiles :
            stderr.write("First pass (record names) - processing file " + f + "\n")
            with open(f, "r") as fi :
                for line in fi :
                    if line.startswith(">") :
                        records.add(gene2recordMapping[line.strip()[1:]])
        stderr.write(str(len(records)) + " records found\n")
    records = list(records)
    # Second pass to process each alignment
    headerBase = ["SNPid", "cluster", "clusterPos", "codonPos", "base", "alt",
                  "nBase", "nAlt", "nonSyn", "multipleRecordEntries"]
    headerGenotypes = ["geno_" + x for x in records]
    headerPositions = ["position_" + x for x in records]
    headerStrand = ["strand_" + x for x in records]
    headers = headerBase + headerGenotypes + headerPositions + headerStrand
    stdout.write("\t".join(headers) + "\n")
    n = str(len(args.alnFiles))
    for (i, f) in enumerate(args.alnFiles) :
        stderr.write("SNP calling - processing file (" + str(i) + "/" + n + ") " + f + "\n")
        SNPdata = pyalign.callSNP(f, gene2recordMapping)
        for SNP in SNPdata :
            o = [str(SNP[x]) for x in headerBase]
            o += [SNP["genotypes"].get(x, "NA") for x in records]
            o += [str(SNP["genomicPositions"].get(x, "NA")) for x in records]
            o += [str(SNP["genomicStrands"].get(x, "NA")) for x in records]
            stdout.write("\t".join(o) + "\n")

### ** Main SNP table

def main_SNPtable(args, stdout, stderr) :
    with open(args.snpTable, "r") as fi :
        headers = fi.next().strip().lstrip("#").split("\t")
        headersOut = headers + []
        if args.strand :
            headersOut = [x for x in headers if not x.startswith("strand_")]
        stdout.write("#" + "\t".join(headersOut) + "\n")
        if args.multipleEntries :
            for line in fi :
                e = dict(zip(headers, line.strip().split("\t")))
                if e["multipleRecordEntries"] == "0" :
                    o = [e[x] for x in headersOut]
                    stdout.write("\t".join(o) + "\n")
        else :
            for line in fi :
                e = dict(zip(headers, line.strip().split("\t")))
                o = [e[x] for x in headersOut]
                stdout.write("\t".join(o) + "\n")
