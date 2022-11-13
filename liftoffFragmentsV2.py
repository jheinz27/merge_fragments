from pyfaidx import Fasta
from Bio import pairwise2
import argparse
import gffutils

parser = argparse.ArgumentParser(description='Check for copies in Liftoff fragments')
parser.add_argument('-i', type=float, required=False, default=1.00, help='Proportion threshold for copy identification')
parser.add_argument('-fi', metavar='<input.gff>', required=True, help='desired input file')
parser.add_argument('-fa', metavar='<target.fasta>', required=True, help='Fasta of target')
parser.add_argument('-min', type=int, metavar='<min_intron_length_bp>', required=False, default=80,help='min intron bp length')
parser.add_argument('-max', type=int, metavar='<max_intron_length_bp>', required=False, default=500000,help='max intron bp length')
parser.add_argument('-o', metavar='<output.gff>', required=True, help='desired output file')

args = parser.parse_args()
threshold = args.i
file = args.fi
fastaFile = args.fa
minIntronLength = args.min
maxIntronLength = args.max
outFile = args.o

target = Fasta(fastaFile)
inputDB = gffutils.create_db(file, file + '.db', merge_strategy="create_unique", force=True, verbose=True)


def removeCopies(genes):
    seqs = []
    hasCopy = [0] * len(genes)
    for tempGene in genes:
        gene = inputDB[tempGene]

        #j = 0
        chr = gene.seqid
        #will have to check how this holds up if no transcripts present
        transcripts = inputDB.children(gene.id,level = 1)
        if transcripts:
            for t in transcripts:
                start = t.start
                end = t.end
               # if j == 0:
                break

            seq = str(target[chr][start:end])
            seqs.append(seq)
        else:
            first = True
            mergedseq = ''
            for c in inputDB.children(gene.id, level = 2):
                typ = ''
                if first:
                    first = False
                    typ = c.featuretype
                if c.featuretype == typ:
                    start = c.start
                    end = c.end
                    curseq = str(target[chr][start:end])
                    mergedseq = mergedseq + curseq
            seqs.append(mergedseq)

    for i in range(len(seqs)):
        tempSeq = seqs[i]
        j = i + 1
        while j < len(seqs):
            toCompareSeq = seqs[j]
            alignments = pairwise2.align.globalxx(tempSeq, toCompareSeq)
            a1 = alignments[0][0]
            a2 = alignments[0][1]
            matchCount = 0
            for base in range(len(a1)):
                if a1[base] == a2[base]:
                    matchCount += 1
            identityProportion = matchCount / (base + 1)
            if identityProportion >= threshold:
                hasCopy[i] = 1
                hasCopy[j] = 1
            j += 1
    nonCopyGenes = []
    copyGenes = []
    for i in range(len(genes)):
        if hasCopy[i] == 0:
            nonCopyGenes.append(genes[i])
        else:
            copyGenes.append(genes[i])
    allInfo = [nonCopyGenes, copyGenes]
    return allInfo


def get_feature_order(gene_db):
    feature_types = list(gene_db.featuretypes())
    print(feature_types)
    index = 0
    feature_order = {}
    if 'exon' in feature_types:
        feature_order['exon'] = index
        index += 1
    if 'CDS' in feature_types:
        feature_order['CDS'] = index
        index += 1
    for feature_type in feature_types:
        if feature_type not in feature_order:
            feature_order[feature_type] = index
            index += 1
    return feature_order

def checkSplits(genes):
    prev = inputDB[genes[0]]
    splits = []
    for i in range(1, len(genes)):
        cur = inputDB[genes[i]]
        prevEnd = prev.end
        curStart = cur.start
        diff = curStart - prevEnd
        if (0 < diff <= minIntronLength):
            splits.append(1)
        elif (diff < maxIntronLength):
            splits.append(2)
        else:
            splits.append(0)
        prev = cur
    return splits


# lots of repetitive code here!!!
def checkChr(genes):
    allChrs = []
    singularChr = []
    for i in range(len(genes)):
        cur = inputDB[genes[i]]
        chr = cur.seqid
        allChrs.append(chr)
    ids = [genes[0]]
    sepIds = []
    prevChr = allChrs[0]
    for i in range(1, len(allChrs)):
        chr = allChrs[i]
        if chr == prevChr:
            ids.append(genes[i])
        else:
            if len(ids) > 1:
                sepIds.append(ids)
            ids = [genes[i]]
        prevChr = chr
    if len(ids) > 1:
        sepIds.append(ids)
    else:
        singularChr.append(ids)
    return [sepIds, singularChr]


idsTest = []
final_parent_list = list(inputDB.features_of_type(featuretype="gene"))
final_parent_list.sort(key=lambda x: (x.seqid, x.start))

featureOrder = get_feature_order(inputDB)

#for key in featureOrder:
    #value = featureOrder[key]
    #print(str(key) + "  " + str(value))

#final_features = final_parent_list.sort(key=lambda x: featureOrder[x.featuretype])

with open(outFile, 'w') as out:

    prevPlusName = ''
    prevNegName = ''
    plusStrandName = ''
    negStrandName=''
    prevPlusGeneID = ''
    prevNegGeneID = ''
    plusInterested = []
    negInterested = []
    plusFirst = 1
    negFirst = 1
    toBePrinted= []
    for gene in final_parent_list:
        geneID = gene.id
        strand = gene.strand
        if strand == '+':
            prevPlusName = plusStrandName
            a = gene.attributes['Name']
            plusStrandName = a[0]
            if plusStrandName == prevPlusName:
                if plusFirst == 1:
                    plusInterested.extend([prevPlusGeneID, geneID])
                    plusFirst = 0
                else:
                    plusInterested.append(geneID)
            else:
                interestedGeneIDs = plusInterested
                plusFirst = 1
                plusInterested = []

            prevPlusGeneID = geneID

        else:
            prevNegName = negStrandName
            a = gene.attributes['Name']
            negStrandName = a[0]
            if negStrandName == prevNegName:
                if negFirst == 1:
                    negInterested.extend([prevNegGeneID, geneID])
                    negFirst = 0
                else:
                    negInterested.append(geneID)

            else:
                interestedGeneIDs = negInterested
                negFirst = 1
                negInterested = []

            prevNegGeneID = geneID

        if interestedGeneIDs:
            all= checkChr(interestedGeneIDs)

            sameChrGenes = all[0]
            possible = all[1]
            if len(possible) > 0:
                toBePrinted.extend(possible)
            print("same Chr Genes   " + str(sameChrGenes))
            for interestedGeneIDsSameChr in sameChrGenes:
                allInfo = removeCopies(interestedGeneIDsSameChr)
                nonCopies = allInfo[0]
                toBePrinted.extend(allInfo[1])
                if len(nonCopies) > 1:
                    print("after copies    " + str(nonCopies))
                    splitTypes = checkSplits(nonCopies)
                    print(splitTypes)
                    allTempGenes = []
                    prev = inputDB[nonCopies[0]]
                    temp = prev
                    prevLev1Children = list(inputDB.children(prev.id,level = 1))
                    prevLev2Children = []
                    for transcript in prevLev1Children:
                        children =  list(inputDB.children(transcript))
                        prevLev2Children.append(children)

                    tempLev1Children = prevLev1Children
                    tempLev2Children = prevLev2Children
                    yes = False

                    for i in range(1, len(nonCopies)):
                        cur = inputDB[nonCopies[i]]
                        splitType = splitTypes[i - 1]
                        curEnd = cur.end
                        curLev1Children = list(inputDB.children(cur.id,level = 1))
                        curLev2Children = []
                        for transcript in curLev1Children:
                            children = list(inputDB.children(transcript))
                            curLev2Children.append(children)
                        if len(curLev2Children) == 0:
                            curLev2Children = curLev1Children
                            curLev1Children = []

                        if splitType == 1:
                            temp.end = cur.end
                  #          if curLev1Children:
                            tempLev1Children[-1].end = curLev1Children[0].end

                            #not going to holdup for multiexon last transcripts!!
                            for i in range(len(tempLev2Children[-1])):
                                tempLev2Children[-1][i].end = curLev1Children[0].end


                        elif splitType == 2:
                            #print(str(tempLev1Children))
                            temp.end = cur.end
                            tempLev1Children[-1].end = curLev1Children[0].end
                            childrenOfMerged = []
                            if prevLev2Children[0]:
                                print(prevLev2Children)
                                originalParent = prevLev2Children[0][0].attributes['Parent']
                                for child in curLev2Children[0]:
                                    child.attributes['Parent'] = originalParent
                                    tempLev2Children[-1].append(child)


                            for i in range(1,len(curLev1Children)):
                                tempLev1Children.append(curLev1Children[i])
                                tempLev2Children.append(curLev2Children[i])

                        #going to have a problem if there is a 0 and more stuff behind it.
                        #thinking I can change the getSplit types function to return separate lists and gene ids of the ones that were cut out.
                        elif splitType == 0:
                            a=1


                    out.write(str(temp) + '\n')
                    for i in range(len(tempLev1Children)):
                        out.write(str(tempLev1Children[i]) + '\n')
                        for j in range(len(tempLev2Children[i])):
                            out.write(str(tempLev2Children[i][j])+ '\n')

                    prevLev2Children = curLev2Children
                    prevLev1Children = curLev1Children

            interestedGeneIDs = []
        else:

            out.write(str(gene)+'\n')
            curChildrenList = list(inputDB.children(geneID, level = 1))
            for child in curChildrenList:

                out.write(str(child) + '\n')
                grandChildren = list(inputDB.children(child))
                if grandChildren:

                    grandChildren.sort(key=lambda x: featureOrder[x.featuretype])
                    for gc in grandChildren:
                        out.write(str(gc) + '\n')

        if len(toBePrinted) > 0:
            current = toBePrinted
            prev = []
            while type(current) is 'list':
                prev = current
                current = current[0]
            print(prev)
            for geneID in prev:
                print(geneID)
                gene = inputDB[geneID]
                out.write(str(gene) + '\n')
                curChildrenList = list(inputDB.children(geneID, level=1))
                for child in curChildrenList:
                    out.write(str(child) + '\n')
                    grandChildren = list(inputDB.children(child))
                    if grandChildren:
                        grandChildren.sort(key=lambda x: featureOrder[x.featuretype])
                        for gc in grandChildren:
                            out.write(str(gc) + '\n')

            toBePrinted = []
