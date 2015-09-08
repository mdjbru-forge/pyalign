### * Description

# Tools to align and refine alignments of fasta sequences

### * Set up

### ** Import

import os
import sys
import shutil
import hashlib
import random
import collections
from Bio import AlignIO
from Bio import SeqIO
from Bio import Align
from Bio import SeqRecord
from Bio import Seq
from Bio import Alphabet
#import pyblastp as pyblastp
#import pycluster as pycluster
import pygenes as pygenes

# The pyblastp and pycluster import were used only for the split function using
# blastp and mcl (code still present in this file but commented out)

### ** Parameters

# from Wikipedia page for "amino acid"
AMINO_ACID_SORTED_MW = ["G", "A", "S", "P", "V", "T", "C", "I", "L", "N", "D",
                        "Q", "K", "E", "M", "H", "F", "R", "Y", "W", "-", "X"]

### * Functions

### ** runMafft(inputFile, out, n = 1)

def runMafft(inputFile, out, n = 1) :
    """Run mafft on a fasta file (perform alignment).
    The awk code to unwrap the fasta output is from Richard Finney:
    http://seqanswers.com/forums/showthread.php?t=27567

    Args:
        inputFile (str): Input fasta file name
        out (str): Output file name
        n (int): Number of threads to use with mafft

    """
    command = ""
    command += "mafft --thread " + str(n)
    command += " " + inputFile
    command += " | awk '{if (substr($0,1,1)==\">\"){if (p){printf \"\\n\";} print $0} else printf(\"%s\",$0);p++;} END {print \"\\n\"}' "
    command += "> " + out
    os.system(command)

### ** ungapFastaFile(inputFile, outputFile)

def ungapFastaFile(inputFile, outputFile) :
    """Ungap a fasta file (remove all "-" characters except in sequence names)

    Args:
        inputFile (str): Path to the input fasta file
        outputFile (str): Path to the output file
    """
    tmp = inputFile + ".ungap.tmp"
    with open(inputFile, "r") as fi :
        with open(tmp, "w") as fo :
            for l in fi :
                if l.strip().startswith(">") :
                    fo.write(l)
                else :
                    fo.write(l.replace("-", ""))
    shutil.move(tmp, outputFile)

### ** ungapAln(aln, gapPropDropPosition = 1, gapPropDropSequence = 0)

def ungapAln(aln, gapPropDropPosition = 1, gapPropDropSequence = 0) :
    """Ungap an alignment by removing positions for which only gaps are present
    
    Args:
        aln (Bio.Align.MultipleSeqAlignment object): An alignment object
        gapPropDropPosition (float): Float between 0 and 1. If a position has
          a proportion of gap at least equal to this value, the position is 
          removed.
        gapPropDropSequence (float): Float between 0 and 1. If a sequence has
          a gap for which gap proportion in this position is less than this
          value, the sequence is dropped.

    Return:
        dict: Dictionary {"ungappedAln" : Bio.Align.MultipleSeqAlignment 
          with gappy positions removed, "removedSeq" : list of removed sequence
          names, "removedPos" : list of removed positions starting at 0}
    """
    # https://www.biostars.org/p/90005/
    # http://biopython.org/DIST/docs/tutorial/Tutorial.html#chapter%3ABio.AlignIO
    # Find indices of gaps
    gaps = []
    seqsToDrop = set([])
    for i in xrange(aln.get_alignment_length()) :
        column = aln[:, i]
        if column.count("-") >= (len(column) * gapPropDropPosition) :
            gaps.append(i)
        if column.count("-") < (len(column) * gapPropDropSequence) :
            seqsToDrop.update([x for x in xrange(len(column)) if column[x] == "-"])
    if len(gaps) == 0 :
        out = aln[:, :]
    else :
        out = aln[:, 0:gaps[0]]
        if len(gaps) == 1 :
            out += aln[:, (gaps[0] + 1):]
        else :
            for i in xrange(len(gaps) - 1) :
                out += aln[:, (gaps[i] + 1):gaps[(i + 1)]]
            out += aln[:, (gaps[-1] + 1 ):]
    # Clear sequences to drop
    outAln = Align.MultipleSeqAlignment([out[x] for x in xrange(len(out)) if x not in seqsToDrop])
    out = dict()
    out["ungappedAln"] = outAln
    out["removedSeq"] = [x for (i, x) in enumerate([x.description for x in aln])
                         if i in seqsToDrop]
    out["removedPos"] = sorted(gaps)
    return out

### ** ungapNtAlnFile(ntFile, removedSeqs, removedPos, outFile)

def ungapNtAlnFile(ntFile, removedSeqs, removedPos, outFile) :
    """Synchronize the ungapping performed on a protein alignment fasta file
    with a detailled nucleotide alignment file. 
    
    Args:
        ntFile (str): Detailled nucleotide alignment file
        removedSeqs (list): List of removed sequence identifiers in the protein
          alignment
        removedPos (list): List of removed positions (starting at 0) in the
          protein alignment
        outFile (str): Output file

    """
    out = ""
    with open(ntFile, "r") as fi :
        geneId = ""
        curOut = ""
        ngaps = 0
        for line in fi :
            if line.startswith(">") :
                if geneId not in removedSeqs :
                    out += curOut
                geneId = line.strip()[1:]
                curOut = ">" + geneId + "\n"
                ngaps = 0
            elif line.strip() != "" :
                content = line.strip().split("\t")
                if int(content[0]) not in removedPos :
                    pepAlnPos = str(int(content[0]) - ngaps)
                    curOut += "\t".join([pepAlnPos] + content[1:]) + "\n"
                else :
                    ngaps += 1
        if geneId not in removedSeqs :
                    out += curOut
    with open(outFile, "w") as fo :
        fo.write(out)

### ** positionConservation(positionAlignment)

def positionConservation(positionAlignment) :
    """Calculate a measure of conservation in one sequence of amino-acids (for
    example, one position (i.e. column) in a multiple alignment).

    Args:
        positionAlignment (iterable): Amino acids sequence at one position of 
          the alignment. Can be a string, list or set.

    Returns:
        float: Conservation score, calculated as the sum of the squared 
          frequencies for each amino acid, which should be equivalent to 
          the probability of picking up two identical amino acid for this 
          position when choosing two sequences in the alignment randomly

    """
    positionFrequencies = iterableFrequency(positionAlignment)
    return(sum([x**2 for x in positionFrequencies.values()]))

### ** conservationProfile(aln)

def conservationProfile(aln) :
    """Calculate a conservation profile for an alignment. For each position in 
    the alignment, a conservation value is calculated with
    :func:`positionConservation`.

    Args:
        aln (``MultipleSeqAlignment``): A MultipleSeqAlignment object from 
          Biopython

    Returns
        list: A list of numerical values corresponding to the conservation 
          indices at each position, calculated with 
          :func:`positionConservation`

    """
    return([positionConservation(aln[:, x]) for x in range(aln.get_alignment_length())])

### ** sequenceConservation(aln)

def sequenceConservation(aln) :
    """Calculate a conservation score for each sequence in an alignment.

    Args:
        aln (``MultipleSeqAlignment``): A MultipleSeqAlignment object from 
          Biopython

    Returns:
        list: Average conservation scores for each sequence of the alignment. 
          The conservation score for each amino acid along a sequence is the
          frequency of this amino acid at this position in the whole alignment.

    """
    compositionPerPosition = [iterableFrequency(aln[:, x])
                              for x in range(aln.get_alignment_length())]
    scores = []
    for sequence in aln :
        score = 0
        aaSeq = str(sequence.seq)
        for (i, p) in enumerate(aaSeq) :
            score += compositionPerPosition[i][p]
        score = score * 1. / aln.get_alignment_length()
        scores.append(score)
    return(scores)

### ** checkDuals(aln, locSepInName, locField)

def checkDuals(aln, locSepInName, locField, stderr = sys.stderr) :
    """Check if an alignment contains several sequences of same origin

    Args:
        aln (``MultipleSeqAlignment``): A MultipleSeqAlignment object from 
          Biopython
        locSepInName (str): Separator used to split each sequence name in the
          alignment
        locField (int): Index corresponding to sequence origin after sequence
          name is split by locSepInName

    Returns:
        bool: False if all sequences have different origins, True if not

    """
    origins = []
    for seq in aln :
        name = seq.description
        origins.append(name.split(locSepInName)[locField])
    counts = set(iterableComposition(origins).values())
    if max(counts) > 1 :
        stderr.write(repr(iterableComposition(iterableComposition(origins).values())) + "\n")
    return max(counts) > 1

### ** makeConsensus(aln, raiseErrorIfDraw)

def makeConsensus(aln, raiseErrorIfDraw = False, stderr = None) :
    """Make a consensus sequence from an alignment
    
    Args:
        aln (``MultipleSeqAlignment``): A MultipleSeqAlignment object from 
          Biopython
        raiseErrorIfDraw (bool): if True, raise an error when there are two or 
          more amino acids with the best count. If False, draws are resolved by 
          taking the lightest amino acid using the AMINO_ACID_SORTED_MW list 
          (made using information from Wikipedia). Note that glutamine and 
          lysine have the same MW, but glutamine is sorted first (neutral vs 
          positive side chain).
        stderr (file): Stderr stream to write warnings when draws are found.
          If None, sys.stderr is used.

    Returns:
        str: Sequence (string) made of the most frequent character for each 
          column in the alignment. If there are two characters or more with 
          the highest frequency, raise an error if asked to.
    """
    if stderr is None :
        stderr = sys.stderr
    countsPerPosition = [iterableComposition(aln[:, x])
                         for x in range(aln.get_alignment_length())]
    consensus = ""
    for (i, counts) in enumerate(countsPerPosition) :
        sortedCounts = sorted(counts.items(),
                              key = lambda x : x[1],
                              reverse = True)
        mostFrequentChar = sortedCounts[0][0]
        mostFrequentCount = sortedCounts[0][1]
        totalCount = sum([x[1] for x in sortedCounts])
        if (len(sortedCounts) > 1) :
            if mostFrequentCount == sortedCounts[1][1] and raiseErrorIfDraw :
                raise Exception("Two characters with same frequency")
            if mostFrequentCount == sortedCounts[1][1] :
                stderr.write("Draw in amino-acid frequencies in position " +
                             str(i + 1) + " (freq = " +
                             str(mostFrequentCount * 1./totalCount) + ")\n")
            candidates = [x[0] for x in sortedCounts if x[1] == mostFrequentCount]
            sortedCandidates = sorted(candidates, key = lambda x : AMINO_ACID_SORTED_MW.index(x))
            mostFrequentChar = sortedCandidates[0]
        consensus += mostFrequentChar
    assert len(consensus) == aln.get_alignment_length()
    return consensus

### ** randomTag(n)

def randomTag(n) :
    """Generate a random tag of size n.

    Args:
        n (int): Length of the tag

    Returns:
        str: A tag of length n, which each character in [0-9][a-f]
    """
    o = ""
    allowed = "0123456789abcdef"
    for x in range(n) :
        o += random.choice(allowed)
    return o

### ** iterableComposition(i)

def iterableComposition(i) :
    """Determine the composition (elements and their counts) of an iterable

    Args:
        i (iterable): Iterable

    Returns:
        dict: Dictionary mapping element identities and their respective 
          counts

    """
    o = dict()
    for e in i :
        o[e] = o.get(e, 0)
        o[e] += 1
    return o

### ** iterableFrequency(i)

def iterableFrequency(i) :
    """Determine the composition (elements and their frequency) of an iterable

    Args:
        i (iterable): Iterable

    Returns:
        dict: Dictionary mapping element identities and their respective 
          frequencies

    """
    o = iterableComposition(i)
    total = sum(o.values())
    for k in o.keys() :
        o[k] = o[k] * 1. / total
    return o

# ### ** splitAlignment(inputFile, outDir = ".", inflation = 1.4, n = 1, remove = True)

# def splitAlignment(inputFile, outDir = ".", inflation = 1.4, n = 1,
#                    remove = True) :
#     """Split a suspicious alignment in fasta format into sub-alignments.
#     First perform a self-blastp, run mcl, extract clusters and realign

#     Args:
#         inputFile (str): Path to the input alignment (fasta file)
#         outDir (str): Path to the output directory
#         inflation (float): Inflation parameter for mcl
#         n (int): Number of parallel processes for blastp
#         remove (bool): Remove the original alignment file
#     """
#     toRemove = []
#     if remove :
#         toRemove.append(inputFile)
#     # Ungap
#     ungapFasta(inputFile, inputFile + ".ungap.tmp")
#     toRemove.append(inputFile + ".ungap.tmp")
#     # Make unique sequences
#     originalInput = SeqIO.parse(inputFile + ".ungap.tmp", "fasta")
#     mapping = dict()
#     hashes = dict()
#     for seq in originalInput :
#         name = seq.description
#         sequence = str(seq.seq)
#         h = hashlib.md5()
#         h.update(sequence)
#         hStr = h.hexdigest()
#         hashes[hStr] = sequence
#         mapping[hStr] = mapping.get(hStr, [])
#         mapping[hStr].append(name)
#     with open(inputFile + ".unique.tmp", "w") as fo :
#         for h, seq in hashes.items() :
#             fo.write(">" + h + "\n")
#             fo.write(seq + "\n")
#     toRemove.append(inputFile + ".unique.tmp")
#     # Run blastp
#     pyblastp.makeBlastDb(inputFile + ".unique.tmp", "prot", inputFile + ".db.tmp")
#     pyblastp.runBlastp(query = inputFile + ".unique.tmp",
#                        db = inputFile + ".db.tmp",
#                        evalMax = 10,
#                        task = "blastp",
#                        out = inputFile + ".blastp.tmp",
#                        max_target_seqs = 10000, cores = n)
#     pyblastp.removeBlastpDb(inputFile + ".db.tmp")
#     toRemove.append(inputFile + ".blastp.tmp")
#     # Run mcl
#     pycluster.prepareBlastpABC(inputFile + ".blastp.tmp",
#                                inputFile + ".blastp.ABC.tmp")
#     pycluster.runMcxload(inputFile + ".blastp.ABC.tmp",
#                          inputFile + ".mci.tmp",
#                          inputFile + ".tab.tmp")
#     pycluster.runMcl(inputFile + ".mci.tmp", inflation,
#                      inputFile + ".tab.tmp", inputFile + ".clusters")
#     toRemove.append(inputFile + ".blastp.ABC.tmp")
#     toRemove.append(inputFile + ".mci.tmp")
#     toRemove.append(inputFile + ".tab.tmp")
#     toRemove.append(inputFile + ".clusters")
#     # Extract clusters
#     clusterFiles = []
#     with open(inputFile + ".clusters", "r") as fi :
#         i = 0
#         for l in fi :
#             if l.strip() != "" :
#                 out = (os.path.join(outDir, os.path.basename(inputFile)) +
#                        ".cl" + str(i) + ".fa")
#                 clusterFiles.append(out)
#                 with open(out, "w") as fo :
#                     for h in l.strip().split("\t") :
#                         for seq in mapping[h] :
#                             fo.write(">" + seq + "\n")
#                             fo.write(hashes[h] + "\n")
#                 i += 1
#     # Align clusters
#     for cluster in clusterFiles :
#         ungapFasta(cluster, cluster + ".ungap.tmp")
#         runMafft(cluster + ".ungap.tmp", cluster + ".tmp")
#         toRemove.append(cluster + ".ungap.tmp")
#         os.rename(cluster + ".tmp", cluster)
#     # Clean tmp files
#     for f in toRemove :
#         os.remove(f)

### ** sequenceConservation(aln)

def sequenceConservation(aln) :
    """Calculate a conservation score for each sequence in an alignment.

    Args:
        aln (``MultipleSeqAlignment``): A MultipleSeqAlignment object from 
          Biopython

    Returns:
        list: Average conservation scores for each sequence of the alignment. 
          The conservation score for each amino acid along a sequence is the
          frequency of this amino acid at this position in the whole alignment.

    """
    compositionPerPosition = [iterableFrequency(aln[:, x])
                              for x in range(aln.get_alignment_length())]
    scores = []
    for sequence in aln :
        score = 0
        aaSeq = str(sequence.seq)
        for (i, p) in enumerate(aaSeq) :
            score += compositionPerPosition[i][p]
        score = score * 1. / aln.get_alignment_length()
        scores.append(score)
    return(scores)

### ** phaseNtDetailled(aln, geneTable)

def phaseNtDetailled(aln, geneTable) :
    """Convert a protein alignment to a detailled nucleotide alignment table

    Args:
        aln (Bio.Align.MultipleSeqAlignment): Alignment with geneId as names
        geneTable (GeneTable from pygenes): Gene table containing the 
          information for the genes whose geneId is in the alignment

    Returns:
        AlnPosTable from pygenes: Detailled alignment data

    """
    alnNt = pygenes.AlnPosTable()
    for seq in aln :
        geneId = seq.description
        pepAln = str(seq.seq)
        geneEntry = geneTable.geneId(geneId)
        alnNt.addAlnRow(geneId, pepAln, geneEntry, ignoreFuzzy = True)
    return alnNt

### ** phaseNtDetailledFast(inputAlnFastaFile, geneTable, outputFile)

def phaseNtDetailledFast(inputAlnFastaFile, geneTable, outputFile) :
    """Convert an input protein alignment to a detailled nucleotide
    alignment file

    Args:
        inputAlnFastaFile (str): Alignment file
        geneTable (GeneTable from pygenes): Gene table containing the 
          information for the genes whose geneId is in the alignment
        outputFile (str): Output file name

    Returns:
        bool: False if no alignment was found, True if it was found
    """
    try :
        aln = AlignIO.read(inputAlnFastaFile, "fasta")
    except ValueError :
        return False
    fo = open(outputFile, "w")
    for row in aln :
        geneId = row.description
        pepAln = str(row.seq)
        geneEntry = geneTable.geneId(geneId)
        ntSeq = geneEntry.codingSeq
        recordPos = pygenes.locStr2int(geneEntry.location, ignoreFuzzy = True)
        assert len(recordPos) == len(ntSeq)
        fo.write(">" + geneId + "\n")
        gaps = 0
        for (i, aa) in enumerate(pepAln) :
            line = [str(i), aa]
            if aa == "-" :
                gaps += 1
            else :
                j = (i - gaps) * 3
                line += [ntSeq[j:(j+3)]]
                line += [",".join([str(x[0]) for x in recordPos[j:(j+3)]])]
                line += ["".join([x[1] for x in recordPos[j:(j+3)]])]
            fo.write("\t".join(line) + "\n")
    fo.close()
    return True
    

### ** prot2nucOld(aln, geneTable)

def prot2nucOld(aln, geneTable) :
    """Convert a protein alignment to a nucleotide alignment

    Args:
        aln (Bio.Align.MultipleSeqAlignment): Alignment with geneId as names
        geneTable (GeneTable from pygenes): Gene table containing the 
          information for the genes whose geneId is in the alignment

    Returns:
        aln (Bio.Align.MultipleSeqAlignment): Nucleotide alignment

    """
    seqReady = []
    for seq in aln :
        gene = geneTable.geneId(seq.description)
        protSeq = str(seq.seq).replace("-", "")
        # Get the truncated peptide seq if needed (due to alignment processing)
        start = gene.peptideSeq.find(protSeq)
        assert start != (-1)
        assert gene.peptideSeq.find(protSeq, start + 1) == (-1)
        nucSeq = gene.codingSeq[(start * 3) : ((start + len(protSeq)) * 3)]
        # Phase gaps
        ntSeqOut = ""
        k = 0
        for l in str(protSeq) :
            if l == "-" :
                ntSeqOut += "---"
            else :
                ntSeqOut += nucSeq[k:k+3]
                k += 3
        seqReady.append(SeqRecord.SeqRecord(Seq.Seq(ntSeqOut,
                                            alphabet = Alphabet.generic_dna)))
    return(Align.MultipleSeqAlignment(seqReady))
                                  

### ** loadAlnFromDetailledNtAln(alnNtFile)

def loadAlnFromDetailledNtAln(alnNtFile) :
    """Load the sequence data from a detailled nucleotide alignment file into
    one protein alignment and one nucleotide alignment

    Args:
        alnFile (str): Path to an alignment file

    Returns:
        tuple: A tuple (nucleotide alignment, protein alignment,
          genomicPositionsList, strandList)

    """
    with open(alnNtFile, "r") as fi :
        protSeqs = []
        nucSeqs = []
        genomicPosList = []
        strandList = []
        protSeq = ""
        nucSeq = ""
        genomicPos = []
        strand = []
        currentGeneId = ""
        for line in fi :
            if line.startswith(">") :
                if currentGeneId != "" :
                    protSeqs.append(SeqRecord.SeqRecord(Seq.Seq(protSeq,
                                                  alphabet = Alphabet.generic_protein),
                                    id = currentGeneId,
                                    description = currentGeneId))
                    nucSeqs.append(SeqRecord.SeqRecord(Seq.Seq(nucSeq,
                                                 alphabet = Alphabet.generic_dna),
                                   id = currentGeneId,
                                   description = currentGeneId))
                    genomicPosList.append(genomicPos)
                    strandList.append(strand)
                currentGeneId = line.strip()[1:]
                protSeq = ""
                nucSeq = ""
                genomicPos = []
                strand = []
            else :
                elements = line.strip().split("\t")
                if len(elements) > 1 :
                    if elements[1] == "-" :
                        protSeq += "-"
                        nucSeq += "---"
                        genomicPos += ["NA", "NA", "NA"]
                        strand += ["NA", "NA", "NA"]
                    else :
                        protSeq += elements[1]
                        nucSeq += elements[2]
                        genomicPos += elements[3].split(",")
                        strand += list(elements[4])
        protSeqs.append(SeqRecord.SeqRecord(Seq.Seq(protSeq,
                                      alphabet = Alphabet.generic_protein),
                        id = currentGeneId,
                        description = currentGeneId))
        nucSeqs.append(SeqRecord.SeqRecord(Seq.Seq(nucSeq,
                                     alphabet = Alphabet.generic_dna),
                       id = currentGeneId,
                       description = currentGeneId))
        genomicPosList.append(genomicPos)
        strandList.append(strand)
    # Build the MultipleSeqAlignment objects (prot and nt)
    protAln = Align.MultipleSeqAlignment(protSeqs)
    nucAln = Align.MultipleSeqAlignment(nucSeqs)
    return (nucAln, protAln, genomicPosList, strandList)

### ** callSNP(alnNtFile, geneToRecordMapping)

def callSNP(alnNtFile, geneToRecordMapping) :
    """Call SNP in a detailled nucleotide alignment file and return the 
    SNPdata as a list.
    NOTE: If a record has several genes involved in the alignment, its genotype
    string is a concatenation of the genotypes.

    Args:
        alnNtFile (str): Path to an alignment file
        geneToRecordMapping (dict): Mapping between gene id and record id

    Returns:
        list: List of SNP data dictionaries
    
    """
    clusterName = os.path.basename(alnNtFile)
    SNPdata = list()
    # Load the alignments
    (nucAln, protAln, genomicPos, strands) = loadAlnFromDetailledNtAln(alnNtFile)
    assert len(nucAln) == len(genomicPos)
    assert len(nucAln) == len(strands)
    # Go through the nt positions to detect SNPs
    for i in xrange(nucAln.get_alignment_length()) :
        nucPosComp = list(collections.Counter(nucAln[:,i]).items())
        multipleOccurrenceOneRecord = False
        if len(nucPosComp) > 1 :
            SNP = dict()
            nucPosComp = sorted(nucPosComp, key = lambda x: x[1], reverse = True)
            SNP["clusterPos"] = int(i / 3)
            SNP["codonPos"] = i % 3
            SNP["base"] = nucPosComp[0][0]
            SNP["nBase"] = nucPosComp[0][1]
            SNP["alt"] = "".join([x[0] for x in nucPosComp[1:]])
            SNP["nAlt"] = sum([x[1] for x in nucPosComp[1:]])
            protPosComp = list(collections.Counter(
                protAln[:, SNP["clusterPos"]]).items())
            SNP["nonSyn"] = int(len(protPosComp) > 1)
            SNP["cluster"] = clusterName
            SNP["SNPid"] = "-".join([clusterName, str(SNP["clusterPos"]),
                                     str(SNP["codonPos"])])
            genotypes = dict()
            genomicPositions = dict()
            genomicStrands = dict()
            for j in range(len(nucAln)) :
                recordId = geneToRecordMapping[nucAln[j].description]
                if recordId in genotypes.keys() :
                    multipleOccurrenceOneRecord = True
                    genotypes[recordId] += "/" + nucAln[j, i]
                    genomicPositions[recordId] += "/" + str(genomicPos[j][i])
                    genomicStrands[recordId] += "/" + strands[j][i]
                else :
                    genotypes[recordId] = nucAln[j, i]
                    genomicPositions[recordId] = str(genomicPos[j][i])
                    genomicStrands[recordId] = strands[j][i]
            SNP["genotypes"] = dict(genotypes)
            SNP["genomicPositions"] = genomicPositions
            SNP["genomicStrands"] = genomicStrands
            SNP["multipleRecordEntries"] = int(multipleOccurrenceOneRecord)
            SNPdata.append(SNP)
    return SNPdata
            
            
### ** mapSequenceToAln(directory)

def mapSequenceToAln(directory) :
    """Build a dictionary mapping sequence name to alignment names for all
    alignment files in a directory

    Args:
        directory (str): Path to the directory to explore

    Returns:
        dict: Mapping (geneId, alnFilename)
    """
    o = dict()
    l = os.listdir(directory)
    for f in l :
        aln = AlignIO.read(f, "fasta")
        for g in aln :
            o[g.description] = f
    return o

### ** splitGeneTable(geneTableFile, mapSeqAln, outDir)

def splitGeneTable(geneTableFile, mapSeqAln, outDir) :
    """Split a large gene table file into smaller files corresponding to the gene
    information for individual alignments

    Args:
        geneTableFile (str): Gene table filename
        mapSeqAln (dict): Mapping between geneId and alignment, output from
          mapSequenceAln
        outDir (str): Path to the output directory

    """
    alnStarted = dict()
    with open(geneTableFile, "r") as fi :
        headerLine = fi.next()
        headers = headerLine.strip().split("\t")
        geneIdIndex = headers.index("geneId")
        for l in fi :
            if l.strip() != "" :
                geneId = l.strip().split("\t")[geneIdIndex]
                aln = mapSeqAln.get(geneId, False)
                if aln :
                    if not alnStarted.get(aln, False) :
                        alnStarted[aln] = True
                        fo = open(os.path.join(outDir, aln + ".geneTable"), "w")
                        fo.write(headerLine)
                        fo.write(l)
                        fo.close()
                    else :
                        fo = open(os.path.join(outDir, aln + ".geneTable"), "a")
                        fo.write(l)
                        fo.close()
