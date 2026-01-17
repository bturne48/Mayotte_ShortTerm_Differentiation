from __future__ import division
import sys
import gzip


# It is assumed that the variants have been mapped to a haplotype 
class StructVariantData:
    # Object of type file that will be run through the class
     def __init__(self, myFile):
          self.myFile = myFile
        

     def readSVFile(self):
        # Input: Read in structural variant file of format:
        # 1) Mutation origin chromosome 
        # 2) Mutation origin starting loci
        # 3) Mutation orgin ending loci
        # 4) Mutation destination chromosomes
        # 5) Mutation destination starting loci
        # 6) Mutation destination ending loci
        # 7) Mutation haplotype
        # 8) Strain of organism
        # Output: A list of all variants with strains added at the end
        
        # * START: Reading in file * #
        # TODO: Make it so you can choose stuff to ignore in try statments ex; pick chroms
        # Read file into memory, might crash on larger files
          with open(self.myFile) as inFile:
            SVList = [line.rstrip() for line in inFile]
          # Make sure file is in correct format
          filteredSV = []
          for x in range(0, len(SVList)):
               splitList = SVList[x].split()
               # Cast columns to expected type, throw warning if wrong
               try:
                    tempString = str(splitList[0]) +'\t'+ str(int(splitList[1])) +'\t'+ str(int(splitList[2])) +'\t'+ str(splitList[3]) \
                                   +'\t'+ str(int(splitList[4])) +'\t'+ str(int(splitList[5])) +'\t'+ str(splitList[6]) +'\t'+ str(splitList[7])
               except ValueError:
                    print('WARNING ValueError in readSVFile:  Removing line ' + SVList[x])
                    continue
               else:
                    # brandon this used to be tempstring, please remember that you changed it
                    filteredSV.append(SVList[x])
          return filteredSV
          # *END: Reading in file * #

     def readHapFile(self):
     # INPUT: List of variants produced from clusterSV function and haplotype file of format:
     # 1) Haplotype chromosome - str
     # 2) Haplotype starting - int
     # 3) Haplotype ending loci - int 
     # 4) Haplotype - str 
     # 5) Haplotype strain - str
     # * START: Reading in file * #
     
     # TODO: Make it so you can choose stuff to ignore in try statments ex; pick chroms
     # Read file into memory, might crash on larger files
          with open(self.myFile) as inFile:
               SVList = [line.rstrip() for line in inFile]
          # Make sure file is in correct format
          filteredSV = []
          for x in range(0, len(SVList)):
               splitList = SVList[x].split()
               # filter headers out
               if splitList[4].startswith("Dsan") or splitList[4].startswith("Dyak"):
                    splitList[4] = splitList[4][5:]
               filteredSV.append(SVList[x])     
               # Cast columns to expected type, throw warning if wrong
          #   try:
          #       tempString = str(splitList[0]) +'\t'+ str(int(splitList[1])) +'\t'+ str(int(splitList[2])) +'\t'+ str(splitList[3]) +'\t'+ str(splitList[4].split("_haps.txt")[0])                
          #   except ValueError:
          #       print('WARNING ValueError in readHapFile:  Removing line ' + SVList[x])
          #       continue
          #   else:
          #     filteredSV.append(tempString)
          # *END: Reading in file * #
          return filteredSV
        

class StructVariants:
     # Must use an object from StructVariantData class
     def __init__(self, structVariantData):
          self.structVariantData = structVariantData
          

     def clusterSV(self, globalChromList):
          # * START: Create bins for SVs * #
          chromDict, chromList = {}, globalChromList
          for x in range(0, len(self.structVariantData)):
               splitSV = self.structVariantData[x].split()
               # needs to be in the chrom list for us to use in anaylsis
               if splitSV[0] not in chromList or splitSV[3] not in chromList:
                    continue
               # reformat the lines so you dont have to worry about checking for the second breakpoint matching
               if (chromList.index(splitSV[0]) > chromList.index(splitSV[3])):
                    # Ex) 2L will always come before 2R on the line
                    splitSV = splitSV[3:6] + splitSV[:3] + splitSV[-2:]
               # within chromosome
               elif chromList.index(splitSV[0]) == chromList.index(splitSV[3]):
                    if int(splitSV[1]) > int(splitSV[4]):
                         splitSV = splitSV[3:6] + splitSV[:3] + splitSV[-2:]
               chrom = splitSV[0]
               originStartLoci = splitSV[1]
               # Create an empty dict for new chromosomes to cluster on
               # Since I'm clustering based on the first origin chrom you don't need to check the second
               if chrom not in chromDict:
                    chromDict[chrom] = {}
               # Create bins based on the mutation origin loci
               # This is a nested dict inside the dict organized by mutation origin chrom
               targetBinDict = chromDict.get(chrom)
               myBin = int(originStartLoci) - int(originStartLoci)%325
               # If a bin key has already been created then append the mutation to a list, otherwise create a new bin  
               if myBin in targetBinDict:
                    prev = targetBinDict.get(myBin)
                    targetBinDict[myBin] = prev + splitSV
               else:
                    targetBinDict[myBin] = splitSV
          # * END: Create bins for SVs * #
        
    	     # * START: Function to check for SVs needing to be clustered * #
          # Function to check if variants inside bins should be clustered
          def multMutantBin(myList):
               # if the mutant only has one match at this locus run this bit
               # returns the variant and an empty list to signify there are no other mutants to scan in this bin
               # I dont feel like rereferencing currentBin cause vi is stupid, deal with this redundancy or perish
               currentBin = myList
               strainList = []
               leftovers = []
               if len(currentBin) == 8:
                    ulLoc = currentBin[0]
                    ulStart = int(currentBin[1])
                    ulEnd = int(currentBin[2])
                    ulLoc2 = currentBin[3]
                    ulStart2 = int(currentBin[4])
                    ulEnd2 = int(currentBin[5])
                    haplo = currentBin[6]
                    strain = currentBin[7]
                    variant = currentBin
                    return [variant, [], ulLoc, ulStart, ulEnd, ulLoc2, ulStart2, ulEnd2, [haplo +'\t'+ strain]]

               else:
                    # run this bit if there are multiple variants in the same bin 
                    # call this in a loop and run it on itself until the leftovers are exhausted and save variants per return
                    strainList = []
                    leftovers = []
                    variant = []
                    # Choose 1 variant at a time to check against
                    # Needs more explanation
                    ulLoc = currentBin[0]
                    ulStart = int(currentBin[1])
                    ulEnd = int(currentBin[2])
                    ulLoc2 = currentBin[3]
                    ulStart2 = int(currentBin[4])
                    ulEnd2 = int(currentBin[5])
                    haplo = currentBin[6]
                    strain = currentBin[7]
                    strainList.append(haplo +'\t'+ strain)
                    # Create references for these variables needed later
                    ulMin = ulStart
                    ulMax = ulEnd
                    ulMin2 = ulStart2
                    ulMax2 = ulEnd2
                    variant = currentBin[0:8]
                    # Loop through the rest of the variants to see if they should be clustered
                    for y in range(8, len(currentBin), 8):
                            nl_ulLoc = currentBin[y]
                            nl_ulStart = int(currentBin[y+1])
                            nl_ulEnd = int(currentBin[y+2])
                            nl_ulLoc2 = currentBin[y+3]
                            nl_ulStart2 = int(currentBin[y+4])
                            nl_ulEnd2 = int(currentBin[y+5])
                            nl_haplo = currentBin[y+6]
                            nl_strain = currentBin[y+7]
                            # If the chroms for origin and dest are the same, and the start and stop loci are within 325bp
                            if ulLoc == nl_ulLoc and ulLoc2 == nl_ulLoc2 and abs(ulStart - nl_ulStart) <=325 and abs(ulEnd - nl_ulEnd) <= 325 \
                                        and abs(ulStart2 - nl_ulStart2) <=325 and abs(ulEnd2 - nl_ulEnd2) <= 325:
                                    # Check if min and max loci need to be updated
                                    if nl_ulStart <= ulMin:
                                            ulMin = nl_ulStart
                                    if nl_ulEnd >= ulMax:
                                            ulMax = nl_ulEnd
                                    if nl_ulStart2 <= ulMin2:
                                            ulMin2 = nl_ulStart2
                                    if nl_ulEnd2 >= ulMax2:
                                            ulMax2 = nl_ulEnd2
                                    # *ASK* If this strain has already been added to the list then we need to skip as we already have an entry 
                                    # Should it be counted twice? 
                                    # ****FLAWED????? REMOVE *JUST* BECAUSE OF MATCHING STRAINS????
                                    if nl_strain in strainList:
                                            continue
                                    else:
                                            strainList.append(nl_haplo +'\t'+ nl_strain)
                                    variant = variant + currentBin[y:y+8]
                            else:
                                    leftovers = leftovers + currentBin[y:y+8]
                    return [variant, leftovers, ulLoc, ulMin, ulMax, ulLoc2, ulMin2, ulMax2, strainList]
          # * END: Function to check for SVs needing to be clustered * #

          # * START: Use the multMutBin function to cluster the variants * #
          # ClustVariants will hold the clusters with the max and min loci and chrom postions as key, and full mutation data as values *fix*
          clustVariants = {}
          for key in chromDict:
               currentChrom = chromDict.get(key)
               # I'm a fucking idiot and hate my past code so i'm reformatting it
               for mybin in currentChrom:
                    def fixBin(dumbBin):
                         stupidBin = currentChrom.get(dumbBin)
                         goodBin = []
                         for i in range(0, len(stupidBin),8):
                              goodBin.append('\t'.join(stupidBin[i:i+8]))
                         currentChrom[dumbBin] = goodBin
                    # It's... It's code that undoes my fixed code..... yup..
                    def unfixBin(goodBin, retBin):
                         badBin = []
                         for i in range(0, len(goodBin)):
                              badBin = badBin + goodBin[i].split('\t')
                         currentChrom[retBin] = badBin
                    fixBin(mybin)
               # Back to big brain code
               # * START: Merge variants in next door bins * # 
               for mybin in currentChrom:
                    # Use this function to see if a specific clustered variant already has that strain associated it when comparing 2 bins
                    def strainVarAlreadyInBin1(bin1, canidateForMerge):
                         for x in range(0, len(bin1)):
                              totalBin = bin1[x].split('\t')
                              if totalBin[0]==canidateForMerge[0] and abs(int(totalBin[1]) - int(canidateForMerge[1]))<=325 and abs(int(totalBin[2]) - int(canidateForMerge[2]))<=325 and\
                                   totalBin[3]==canidateForMerge[3] and abs(int(totalBin[4]) - int(canidateForMerge[4]))<=325 and abs(int(totalBin[5]) - int(canidateForMerge[5]))<=325 and\
                                   totalBin[7]==canidateForMerge[7]:
                                   return True
                              else:
                                   return False
                    # Through merging the bin may now be empty, if so skip
                    if currentChrom.get(mybin) == []:
                         continue
                    # Check if neigboring bins have variants that should cluster
                    # Feed this the values of the currentChromDict.get(mybin) and currentChromDict.get(mybin+/-325)
                    def mergeBin(firstBin, secondBin):
                         len1 = len(firstBin)
                         len2 = len(secondBin)
                         bin1 = firstBin
                         bin2 = secondBin
                         correctedBin1 = []
                         correctedBin2 = []
                         for x in range(0, len(bin1)):
                              splitBin1 = bin1[x].split('\t')
                              # We are merging on bin1, so anything in bin1 stays
                              if bin1 in correctedBin1:
                                   continue
                              else:
                                   correctedBin1.append(bin1[x])
                              for y in range(0, len(bin2)):
                                   splitBin2 = bin2[y].split('\t')
                                   # Variants to cluster in different bins that are from different strains
                                   if splitBin1[0]==splitBin2[0] and abs(int(splitBin1[1]) - int(splitBin2[1]))<=325 and abs(int(splitBin1[2]) - int(splitBin2[2]))<=325 and\
                                   splitBin1[3]==splitBin2[3] and abs(int(splitBin1[4]) - int(splitBin2[4]))<=325 and abs(int(splitBin1[5]) - int(splitBin2[5]))<=325 and\
                                   splitBin1[7]!=splitBin2[7]:
                                        # Annoyingly, you have to check if this variant+strain maybe slightly different and is already in the bin
                                        if strainVarAlreadyInBin1(bin1, splitBin2)==True:
                                             continue
                                             # print('\n')
                                        # If not append to bin 1
                                        else:
                                             if bin2[y] in correctedBin1:
                                                  continue
                                             else:
                                                  if bin2[y] in correctedBin2:
                                                       correctedBin2.remove(bin2[y])
                                                  correctedBin1.append(bin2[y])
                                   # Variants that dont match when parsed through should stay in their bin
                                   else:
                                        if bin2[y] in correctedBin2 or bin2[y] in correctedBin1:
                                             continue
                                        else:
                                             correctedBin2.append(bin2[y])
                         len3 = len(correctedBin1)
                         len4 = len(correctedBin2)
                         # if len1+len2 != len3+len4:
                         #      print('Your code sucks')
                         # Update the bin with the new merged bin
                         return [correctedBin1, correctedBin2]
                    if mybin-325 in currentChrom:
                         postMerge = mergeBin(currentChrom.get(mybin),currentChrom.get(mybin-325))
                         currentChrom[mybin] = postMerge[0]
                         currentChrom[mybin-325] = postMerge[1]
                    if mybin+325 in currentChrom:
                         postMerge = mergeBin(currentChrom.get(mybin),currentChrom.get(mybin+325))
                         currentChrom[mybin] = postMerge[0]
                         currentChrom[mybin+325] = postMerge[1]
                    # * END: Merge variants in next door bins * # 
               for mybin in currentChrom:
                    # Missing bins from merging code
                    if currentChrom.get(mybin) == []:
                         continue
                    # BIG.BRAIN.CODING
                    unfixBin(currentChrom.get(mybin), mybin)
                    # * START: Run function to cluster variants on all bins * #
                    # For each bin in each chromosome     
                    roundN = multMutantBin(currentChrom.get(mybin))
                    # If there are no remaining varaints left in a bin then you just need to call the multMutBin once
                    if roundN[1] == []:
                         clustVariants[roundN[2]+'\t'+str(roundN[3])+'\t'+str(roundN[4])+'\t'+roundN[5]+'\t'+str(roundN[6])+'\t'+str(roundN[7])] = '\t'.join(roundN[8])
                    # If there are still unclustered variants then keep calling the function until there are none remaining	
                    else:
                         while roundN[1] != []:
                              clustVariants[roundN[2]+'\t'+str(roundN[3])+'\t'+str(roundN[4])+'\t'+roundN[5]+'\t'+str(roundN[6])+'\t'+str(roundN[7])] = '\t'.join(roundN[8])
                              roundN = multMutantBin(roundN[1])
                         clustVariants[roundN[2]+'\t'+str(roundN[3])+'\t'+str(roundN[4])+'\t'+roundN[5]+'\t'+str(roundN[6])+'\t'+str(roundN[7])] = '\t'.join(roundN[8])
                    # * END: Run function to cluster variants on all bins * #
          return clustVariants
        # # * END: Use the multMutBin function to cluster the variants * #


     def newCluster(self):
          chromDict, chromList = {'2L':{}, '2R':{}, '3L':{}, '3R':{}, 'X':{}}, ['2L', '2R', '3L', '3R', 'X']
          # manually create bins to make sure there are no errors later
          for chrom in chromDict:
               for x in range(0, 30000000, 325):
                    chromDict[chrom][x] = []
          
          # read in sv data to the bins
          for x in range(0, len(self.structVariantData)):
               splitSV = self.structVariantData[x].split('\t')       
               # reformat the lines so you dont have to worry about checking for the second breakpoint matching
               if (chromList.index(splitSV[0]) > chromList.index(splitSV[3])):
                    # Ex) 2L will always come before 2R on the line
                    splitSV = splitSV[3:6] + splitSV[:3] + splitSV[-2:]
               # within chromosome
               elif chromList.index(splitSV[0]) == chromList.index(splitSV[3]):
                    if int(splitSV[1]) > int(splitSV[4]):
                         splitSV = splitSV[3:6] + splitSV[:3] + splitSV[-2:]
               # binning vars
               binChrom1, binPos1 = splitSV[0], int(splitSV[1])
               binNum1 = binPos1 - binPos1%325
               targetChrom1 = chromDict[binChrom1]
               # add line to bin
               targetChrom1[binNum1].append('\t'.join(splitSV))

          # check the bins for BOTH breakpoints to see if there is overlap at the start/stop pos
          newClusterDict, chromList = {'2L':{}, '2R':{}, '3L':{}, '3R':{}, 'X':{}}, ['2L', '2R', '3L', '3R', 'X']
          # manually create bins to make sure there are no errors later
          for chrom in newClusterDict:
               for x in range(0, 30000000, 325):
                    newClusterDict[chrom][x] = []

          # since both breakpoints are put into the dictionary above, need to make sure mut calls are only used once
          mutCallUsedList = []
          
          # loop over nested dict
          for chrom in chromDict:
               for bins in chromDict[chrom]:
                    # if the bin has no mutations continue 
                    if chromDict[chrom][bins] == []:
                         continue
                    
                    # mut calls in the bin for this loop
                    combineMutListLoop = chromDict[chrom][bins]
                    for mut in combineMutListLoop:
                         if mut not in mutCallUsedList:
                              mutCallUsedList.append(mut)

                    # use this to get flanking bins up and down stream
                    def adjacentBins(ogMutList, ogBin, asc=True):
                         # vars needed for first loop 
                         prevMutCallListLength = 1
                         binIter = ogBin

                         while prevMutCallListLength != len(ogMutList):
                              newMutList = []
                              # update the lenght of the previous loop to see if any more muts get added for this loop
                              prevMutCallListLength = len(ogMutList)
                              # check if we are ascendign or descending
                              if asc==True:
                                   binIter = binIter+325
                              else:
                                   binIter = binIter-325
                              # check the adjacent bins
                              if binIter in chromDict[chrom]:
                                   # check if this mut call should be added by looking at muts in existing bin and if it matches with bins nextdoor
                                   for mut in ogMutList:     
                                        splitMut = mut.split()
                                        chrom1, start1, end1 = splitMut[0], int(splitMut[1]), int(splitMut[2])
                                        chrom2, start2, end2 = splitMut[3], int(splitMut[4]), int(splitMut[5])
                                        for newMuts in chromDict[chrom][binIter]:
                                             # if you've already used this mut call
                                             if newMuts in mutCallUsedList:
                                                 continue
                                             splitNewMut = newMuts.split()
                                             new_chrom1, new_start1, new_end1 = splitNewMut[0], int(splitNewMut[1]), int(splitNewMut[2])
                                             new_chrom2, new_start2, new_end2 = splitNewMut[3], int(splitNewMut[4]), int(splitNewMut[5])
                                             # check for overlap bewteen bin nextdoor and current bin
                                             if chrom1 == new_chrom1 and abs(start1-new_start1)<=325 and abs(start2-new_start2)<=325 and\
                                                  chrom2 == new_chrom2 and abs(end1-new_end1)<=325 and abs(end2-new_end2)<=325:
                                                  newMutList.append(newMuts)
                                                  # since you have merged this mut call, you need to remove it so you dont dup it
                                                  chromDict[chrom][binIter].remove(newMuts)
                                                  mutCallUsedList.append(newMuts)

                              # add new muts to og list
                              ogMutList += newMutList
                              
                         # add all mut calls to left and right
                         return ogMutList

                    # append all muts in the adjacent bins to the right (ascending bins)
                    combineMutListUpper = adjacentBins(combineMutListLoop, bins)
                    # append all muts in the adjacent bins to the left (descending)
                    combineMutListLower = adjacentBins(combineMutListUpper, bins, False)

                    # append all merged bins back to a new dict
                    for mut in combineMutListLower:
                         if mut not in newClusterDict[chrom][bins]:
                              newClusterDict[chrom][bins].append(mut)
               # NOTE TEMP
               #break

          # clustering of stuff in bins
          for chrom in newClusterDict:
               for bins in newClusterDict[chrom]:
                    # dont iterate on empty bins
                    if newClusterDict[chrom][bins] == []:
                         continue

                    # intial list
                    startMutList = newClusterDict[chrom][bins]

                    strainList = []
                    # cluster stuff inside bins, first breakpoint is already clustered no matter what, it's just seperating out different second bps
                    def clusterUntilDone(mutList):
                         sortDict = {}
                         # results of using this as a search param
                         splitStartMut = mutList[0].split()
                         #print('init mut:' + mutList[0])
                         chrom1, start1, end1 = splitStartMut[0], int(splitStartMut[1]), int(splitStartMut[2])
                         chrom2, start2, end2 = splitStartMut[3], int(splitStartMut[4]), int(splitStartMut[5])
                         hap, strain = splitStartMut[-2], splitStartMut[-1]
                         
                         # sort any muts that have the same second bp chrom to check for final clustering
                         sortDict = {start2: mutList[0].split()}
                         for x in range(1, len(mutList)):
                              n_splitStartMut = mutList[x].split()
                              n_chrom1, n_start1, n_end1 = n_splitStartMut[0], int(n_splitStartMut[1]), int(n_splitStartMut[2])
                              n_chrom2, n_start2, n_end2 = n_splitStartMut[3], int(n_splitStartMut[4]), int(n_splitStartMut[5])
                              n_strain = n_splitStartMut[-1]
                              # only look at the ones with correct chroms to save iteration time
                              if chrom2==n_chrom2:
                                   # make sure you are not overwriting any keys
                                   if n_start2 not in sortDict:
                                        sortDict[n_start2] = n_splitStartMut
                                   else:
                                        sortDict[n_start2] += n_splitStartMut
                         
                         # do clustering for all muts that match the inital mutation (first index)
                         sortKeys = sortDict.keys()
                         sortKeys.sort()
                         # need to make min and max as to properly cluster everything 
                         minStart1, maxEnd1, minStart2, maxEnd2 = start1, end1, start2, end2
                         for key in sortKeys:
                              sortedMut = sortDict[key]
                              # list of muts to remove once done because they were clustered
                              #print(key, sortedMut)
                              for x in range(0, len(sortedMut), 8):
                                   sortStart1, sortEnd1, sortStart2, sortEnd2 = int(sortedMut[x+1]), int(sortedMut[x+2]), int(sortedMut[x+4]), int(sortedMut[x+5]),
                                   n_hap, n_strain =  sortedMut[x+6], sortedMut[x+7]
                                   # check for 325+/- overlap of muts with first index mut for the second breakpoint
                                   if abs(sortStart2 - minStart2)<=325 and abs(sortEnd2 - maxEnd2)<=325:
                                        #print('match mut:' + '\t'.join(sortedMut[x:x+8]))
                                        # update min and max for second bp
                                        minStart2, maxEnd2 = sortStart2, sortEnd2
                                        # update first bp min max
                                        if sortStart1 < minStart1:
                                             minStart1 = sortStart1
                                        if sortEnd1 > maxEnd1:
                                             maxEnd1 = sortEnd1
                                        # append hap and strain and to remove list
                                        if n_strain not in strainList:
                                             strainList.append(n_hap)
                                             strainList.append(n_strain)
                                        # even if we didnt add to the strain list, we dont need that mut call anymore
                                        #print('rem:\t'+'\t'.join(sortedMut[x:x+8]))
                                        startMutList.remove('\t'.join(sortedMut[x:x+8]))


                         # append the init mut strain and hap to a list and add to removal
                         if strain not in strainList:
                              strainList.append(hap)
                              strainList.append(strain)


                         #print('cluster')
                         print(chrom1 +'\t'+ str(minStart1)+'\t'+ str(maxEnd1)+'\t'+ chrom2+'\t'+ str(minStart2)+'\t'+ str(maxEnd2) +'\t'+ '\t'.join(strainList))

                    # run clustering until there is only 1 thing or nothing left in the bin
                    prevLen = 0
                    while len(startMutList)!=prevLen and len(startMutList)>1:
                         prevLen = len(startMutList)
                         #print('start list:')
                         #print(startMutList)
                         clusterUntilDone(startMutList)
                         #print('after list:')
                         #print(startMutList)
                    #print('next mut\n')

                    # # if there was one thing left we need to print that as a singleton cluster
                    if len(startMutList) != 0:
                         for mut in startMutList:
                              print(mut)
                    #print('next\n')


     def localVsGenomicCov(self, clusteredVars, bamPath):
          # import/modules
          import os, sys, subprocess, re

          # BAM location
          bamDir = os.listdir(bamPath)
          bamDir = [bam for bam in bamDir if bam.endswith('.bam')]

          # loop over variants to find coverage increase to determine which part of variant is donor region
          counter = 1
          resultsList = []
          for variant in clusteredVars:
               #print(variant)
               counter += 1
               splitVar, varName = variant.split(), '\t'.join(variant.split()[:6])
               var1, var2, strains  = splitVar[:3], splitVar[3:6], splitVar[7::2]
               # find average coverage of region +/- 1kb, my var should be [chrom,startpos,endpos]
               def calcAvgCov(myVar, bpThreshold):
                    # loop over bam files in scratch and match them to strain names with regex
                    covDict = {}
                    for strain in strains:
                         for bamStrain in bamDir:
                              
                              test_strain, test_bamStrain =  strain, bamStrain.split('.')[1]
                              # if strain starts with cy/ny it doesnt have this strain and lane trash notation to filter out
                              # if (bamStrain.lower()).startswith('cy') or (bamStrain.lower()).startswith('ny'):
                              #      test_strain, test_bamStrain = ''.join(strain.split("_")).lower(), ''.join(bamStrain.split("_")).lower()  
                              # else:
                              #      # split by lane so you can remove that and the strain, to properly match strain names with regex 
                              #      test_strain, test_bamStrain = ''.join(strain.split("_")).lower(), ''.join(bamStrain.split('_L00')[0].split('_')[:-1]).lower()  
                              # check for name match with regex
                              if re.match(r'\b' + test_strain + r'\b', test_bamStrain) and bamStrain.endswith(".bam"):
                                   # uses samtools to pull the coverage of the region
                                   covOut1 = subprocess.check_output('samtools depth -r ' +myVar[0]+ ':' +str(max(1,int(myVar[1])-bpThreshold))+ '-' +str(int(myVar[1]))+ " " +bamPath+ bamStrain\
                                                                    +"| awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }'", shell=True)
                                   covOut2 = subprocess.check_output('samtools depth -r ' +myVar[0]+ ':' +str(max(1, int(myVar[2])))+ '-' +str(int(myVar[2])+bpThreshold)+ " " +bamPath+ bamStrain\
                                                                    +"| awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }'", shell=True)
                                   # if there was no region to create coverage set it to 0.00001 to prevent divide by 0 error 
                                   if covOut1 == '':
                                        covOut1 = '0.00001'
                                   if covOut2 == '':
                                        covOut2 = '0.00001'
                                   # create a dict with strain name and avg coverage to compare variants with mulitple strains
                                   covDict[strain] = [varName, {myVar[0] +' '+ str(max(1, int(myVar[1])-bpThreshold)) +'-'+ myVar[1]:covOut1.rstrip(), \
                                                                 myVar[0] +' '+ myVar[2] +'-'+ str(max(1,int(myVar[2])+bpThreshold)):covOut2.rstrip()}, '\t'.join(myVar)]
                                   # you already found the match so you can go to the next strain in the variants strain list
                                   break
                              else: # testing whether there are bad bam names
                                   pass
                    if covDict=={}: # debug 
                         pass
                         #print(strains)
                    return covDict
               # store the results of each variant and both of their sides into a list of lists, each element of parent list is a variant, child is the sides
               resultsList.append([calcAvgCov(var1, 100), calcAvgCov(var2, 100)])
               #if counter == 100:
                #    break
          return resultsList


     def findDonorRegion(self, localVsGenomicCovOutput, covTable, filterTE=''):
          # * TEMP * Check how TE results of stuff 
          teDict = {}
          path = '/users/bturne48/RogersLab/RemoteEdit/ParseIlluminaData/Data/OutputFromScripts/'
          outNone = open(path+'outNoneSides_100.txt', 'w')
          outOne = open(path+'outOneSides_100.txt', 'w')
          outTwo = open(path+'outTwoSides_100.txt', 'w')
          with open(filterTE) as inTE:
               for line in inTE:
                    teVarName = '\t'.join(line.split()[-6:])
                    # put it in a dictionary that specifies either one sided or two sided
                    if len(line.split())==9:
                         teSide = '\t'.join(line.split()[:3])
                         teDict[teVarName] = [1, teSide]
                    else:
                         teSide = '\t'.join(line.split()[:6])
                         teDict[teVarName] = [2, teSide]

          # compare the average coverage of both sides to the genomic coverage of the strain
          with open(covTable) as genomicCovTable:
               genCovTable = {}
               for line in genomicCovTable:
                    if line != '' and len(line.split())>2:
                         # if strain starts with cy/ny it doesnt have this strain and lane trash notation to filter out
                         if (line.lower()).startswith('cy') or (line.lower()).startswith('ny'): 
                              genCovTable[line.split()[0].lower()] = line.split()[1:]
                         else:
                              genCovTable[''.join(line.split('_L00')[0].split('_')[:-1]).lower()] = line.split()[1:]

          def genomicVsLocalCov(avgCovDict):
               outList, foldList  = (), ()
               # loop over dict that has table info and the dict that has the local cov
               for strainKey in avgCovDict:
                    # row already in .lower() for mat
                    for row in genCovTable:
                         # if strains match from genomic cov table and local coverage
                         if row.rstrip() == ''.join(strainKey.split('_')).lower().rstrip():
                              # track full name of the variant in question
                              varSide, varName = avgCovDict[strainKey][2], avgCovDict[strainKey][0]
                              # the list of chrom coverage from the table
                              splitRow = genCovTable[row]
                              # this is the format of the input
                              chrom2Rcov, chrom2Lcov, chrom3Rcov, chrom3Lcov, chromXcov = float(splitRow[0]), float(splitRow[1]), float(splitRow[2]), float(splitRow[3]), float(splitRow[4])
                              covList = {'2R':chrom2Rcov, '2L':chrom2Lcov, '3R':chrom3Rcov, '3L':chrom3Lcov, 'X':chromXcov}
                              for side in avgCovDict[strainKey][1]:
                                   # which chrom from the table are we pulling in relation to my variant
                                   targetChrom = avgCovDict[strainKey][0].split('\t')[0]
                                   # pull the correct chrom from the table using the parent loops chromsosome as the .get() search term
                                   tableCov, localCov = float(covList[targetChrom]), float(avgCovDict[strainKey][1].get(side))
                                   # check to see if local cov is 1.25x < x < 1.75 greater than table cov
                                   #print('Variant Side: ' + side, 'Strain: '+row, 'Local Cov: '+str(localCov), 'Table Cov: '+ str(tableCov), 'Fold Change: ' + str(localCov/tableCov))
                                   #print('\n')
                                   if localCov/tableCov >= 1.25:# and localCov/tableCov <= 1.75:
                                        outList += (1,)
                                        foldList += (localCov/tableCov, )
                                   else:
                                        outList += (0,)
                                        foldList += (localCov/tableCov, )
               return [outList, foldList, varSide, varName]
          
          # parse over the output of localVsGenomic, list of variants and each element is a dictionary of both sides coverage for each strain
          gList = []
          for variant in localVsGenomicCovOutput:
               sidesList = []
               for sideDict in variant:
                    if sideDict=={}: # temp debug CAR_1566, missing that strains info
                         continue 
                    sidesList.append(genomicVsLocalCov(sideDict))
               if sidesList==[]: # temp debug car
                    continue
               gList.append(sidesList)
          
          # # for each variant check if it has TE stuff
          byeList = []
          for variant in gList:
               g1, g2= variant[0], variant[1]
               varName = g1[2] +'\t'+ g2[2]
               #print(g1)
               #print(g2[2])
               #print('\n')
               #if you are looking at 0-side TE, duplication events that are non-TE
               if varName not in teDict:
                    if sum(g1[0])>sum(g2[0]):
                         outNone.write(g1[2]+'\t'+varName+'\tnonambig\n')
                    elif sum(g1[0])<sum(g2[0]):
                         outNone.write(g2[2]+'\t'+varName+'\tnonambig\n')
                    else:
                         #print('0', str(sum(g1[0])), str(sum(g2[0])), g1[1], g2[1], varName)
                         outNone.write(varName+'\tambig\n')
                    # check if they aare TEs in coverage

               # if you are looking at 1-side TE
               if varName in teDict and teDict[varName][0]==1:
                    if sum(g1[0])>sum(g2[0]):
                         outOne.write(g1[2]+'\t'+varName+'\tnonambig\n')
                    elif sum(g1[0])<sum(g2[0]):
                         outOne.write(g2[2]+'\t'+varName+'\tnonambig\n')
                    else:
                         #print('0', str(sum(g1[0])), str(sum(g2[0])), g1[1], g2[1], varName)
                         outOne.write(varName+'\tambig\n')    

               # if you are looking at 2-side TE, ectopic recomb? 
               if varName in teDict and teDict[varName][0]==2:
                    if sum(g1[0])!=0 and sum(g2[0])!=0:
                         #print('2', str(sum(g1[0])), str(sum(g2[0])), g1[1], g2[1], teDict[varName][1], varName)
                         outTwo.write(varName+'\tnonambig\n')
                    else:
                         #print('2', str(sum(g1[0])), str(sum(g2[0])), g1[1], g2[1], teDict[varName][1], varName)
                         outTwo.write(varName+'\tambig\n')
               #print('\n')
          return byeList
               

     def heteroCorrection(self, localVsGenomicCovOutput, covTable, chromListSnake):
          # compare the average coverage of both sides to the genomic coverage of the strain
          with open(covTable) as genomicCovTable:
               genCovTable = {}
               for line in genomicCovTable:
                    if line.startswith('Chrom'):
                         continue
                    genCovTable[line.split()[0]] = line.split()[1:]
                    # THIS ONLY WORKS FOR DSAN STUFF
                    # if line != '' and len(line.split())>2:
                    #      # if strain starts with cy/ny it doesnt have this strain and lane trash notation to filter out
                    #      if (line.lower()).startswith('cy') or (line.lower()).startswith('ny'): 
                    #           genCovTable[line.split()[0].lower()] = line.split()[1:]
                    #      else:
                    #           genCovTable[''.join(line.split('_L00')[0].split('_')[:-1]).lower()] = line.split()[1:]
                    #           #print(''.join(line.split('_L00')[0].split('_')[:-1]).lower())
          chromList = chromListSnake

          def genomicVsLocalCov(avgCovDict):
               foldList, varName  = [], []
               #print(avgCovDict)
               outList, outStrain = [], ''
               varSide = ''
               # loop over dict that has table info and the dict that has the local cov
               for strainKey in avgCovDict:
                    #outList = []
                    for row in genCovTable:
                         # if strains match from genomic cov table and local coverage
                         if row.lower() == strainKey.lower():
                              outStrain = strainKey
                              # track full name of the variant in question
                              varSide = avgCovDict[strainKey][2]
                              # the list of chrom coverage from the table
                              splitRow = genCovTable[row]
                              #print(splitRow)
                              
                              covList = {}
                              #print(len(chromList), len(splitRow))
                              for x in range(0, len(chromList)):
                                   covList.update({chromList[x]:float(splitRow[x])})
                              #print(covList)
                              # this is the format of the input
                              # ONLY WORKS FOR DSAN
                              #chrom2Rcov, chrom2Lcov, chrom3Rcov, chrom3Lcov, chromXcov = float(splitRow[0]), float(splitRow[1]), float(splitRow[2]), float(splitRow[3]), float(splitRow[4])
                              #covList = {'2R':chrom2Rcov, '2L':chrom2Lcov, '3R':chrom3Rcov, '3L':chrom3Lcov, 'X':chromXcov}
                              
                              # per each breakpoint
                              for side in avgCovDict[strainKey][1]:
                                   #print(row)
                                   #print(side)
                                   localChrom = varSide.split()[0]
                                   #print(localChrom)
                                   #print(avgCovDict[strainKey][1])

                                   # skip if missing info at sites
                                   if avgCovDict[strainKey][1].get(side) == b'':
                                        continue

                                   # pull the correct chrom from the table using the parent loops chromsosome as the .get() search term
                                   tableCov, localCov = float(covList.get(localChrom)), float(avgCovDict[strainKey][1].get(side))
                                   #print('Variant Side: ' + side, 'Strain: '+row, 'Local Cov: '+str(localCov), 'Table Cov: '+ str(tableCov), 'Fold Change: ' + str(localCov/tableCov))
                                   # hetero bounds 
                                   if localCov/tableCov >= 1.3 and localCov/tableCov<=1.7:# and localCov/tableCov <= 1.75:
                                        #print('hetero')
                                        #print(avgCovDict[strainKey][0] +'\t'+ side +'\t'+ row +'\t'+ str(localCov) +'\t'+ str(tableCov) +'\t'+ str(localCov/tableCov) +'\t'+ 'hetero')
                                        outList.append('hetero')
                                        foldList.append(str(localCov/tableCov))
                                        continue
                                   # inbred bounds
                                   if localCov/tableCov >= 1.8 and localCov/tableCov<=2.2:
                                        #print('homo')
                                        #print(avgCovDict[strainKey][0] +'\t'+ side +'\t'+ row +'\t'+ str(localCov) +'\t'+ str(tableCov) +'\t'+ str(localCov/tableCov) +'\t'+ 'homo')
                                        outList.append('homo')
                                        foldList.append(str(localCov/tableCov))
                                        continue
                                   # te bounds
                                   if localCov/tableCov >= 3:
                                        #print('te')
                                        #print(avgCovDict[strainKey][0] +'\t'+ side +'\t'+ row +'\t'+ str(localCov) +'\t'+ str(tableCov) +'\t'+ str(localCov/tableCov) +'\t'+ 'te')
                                        outList.append('te')
                                        foldList.append(str(localCov/tableCov))
                                        continue
                                   if localCov/tableCov <= 0.25:
                                        #print('del/inv')
                                        #print(avgCovDict[strainKey][0] +'\t'+ side +'\t'+ row +'\t'+ str(localCov) +'\t'+ str(tableCov) +'\t'+ str(localCov/tableCov) +'\t'+ 'del/inv')
                                        outList.append('del/inv')
                                        foldList.append(str(localCov/tableCov))
                                        continue
                                   #print('ambig')
                                   #print(avgCovDict[strainKey][0] +'\t'+ side +'\t'+ row +'\t'+ str(localCov) +'\t'+ str(tableCov) +'\t'+ str(localCov/tableCov) +'\t'+ 'ambig')
                                   outList.append('ambig')
                                   foldList.append(str(localCov/tableCov))
                              #print(avgCovDict[strainKey][0] +'\t'+ '\t'.join(outList) +'\t'+ row)
                              #print('\n')
               return [outList, foldList, varSide, outStrain]

          # parse over the output of localVsGenomic, list of variants and each element is a dictionary of both sides coverage for each strain
          gList = []
          for variant in localVsGenomicCovOutput:
               sidesList = []
               for sideDict in variant:
                    sidesList.append(genomicVsLocalCov(sideDict))
               gList.append(sidesList)
          
          # # for everyvariant in list, check strains coverage stuff
          outList = []
          for variant in gList:
               print(variant)
               for side in variant:
                    # assign a correction if needed
                    predictedTypeList = side[0]
                    if ('homo' in predictedTypeList and 'ambig' in predictedTypeList) and 'te' not in predictedTypeList:
                         #print(variant[0][2] +'\t'+ variant[1][2] +'\t'+ side[3])
                         #print(predictedTypeList)
                         outList.append(variant[0][2] +'\t'+ variant[1][2] +'\t'+ side[3])
                         #print('\n')
               #print('\n')
          return outList


     def checkSampledHaps(self, clusteredVars, structVariantData_Hap):
          import re

          variantList = []
          sampledList = []
          
          # Load in clustered variants from previous class object
          for variant in clusteredVars:
               variantList.append(variant +'\t'+ clusteredVars.get(variant))

          # Compare the the variants to the haplotypes to look for sampled haplotypes
          mutationList = []
          for x in range(0,len(variantList)):
               # A list of variants, their location, and all strains/haplotypes associated with it
               splitVariants = variantList[x].split()
               ulLoc = splitVariants[0]
               ulStart = splitVariants[1]
               ulEnd = splitVariants[2]
               ulLoc2 = splitVariants[3]
               ulStart2 = splitVariants[4]
               ulEnd2 = splitVariants[5]
               obsvStrains = splitVariants[6:]

               # Create a list for the variant and the strains that were sampled to check if all haplotypes get listed, fixes the dsan inbred issues seen with map.py
               inStrainList = []
               sampStrainList = []
               for q in range(1,len(obsvStrains), 2):
                    #print((''.join(obsvStrains[q].split('_')).lower()))
                    inStrainList.append(''.join(obsvStrains[q].split('_')).lower())

               # ouput lists
               outList = [ulLoc, ulStart, ulEnd, ulLoc2, ulStart2, ulEnd2]
               # map to haplotypes and look for "new" and/or longest region
               tempMutationList = []
               newMutationList = []
               # fix i hate myself
               tempList = []
               dsanStrains = []
               for y in range(0, len(structVariantData_Hap)):
                    splitHaplo = structVariantData_Hap[y].split()
                    haploLoc = splitHaplo[0]
                    haploStart = splitHaplo[1]
                    haploEnd = splitHaplo[2]
                    haploType = splitHaplo[3]
                    haploStrain = splitHaplo[-1]
                    # search for dsan strains 

                    if re.search('^(?:(?!CY|cy|NY|ny|Grey|grey|Oran|oran).)+$', haploStrain):
                         if haploStrain not in dsanStrains:
                              dsanStrains.append(haploStrain)
                    # fix so you can compare names
                    if haploStrain.startswith("Dsan") or haploStrain.startswith("Dyak"):
                         haploStrain = haploStrain[5:]
                    
                    # map to haplotypes
                    if (haploLoc==ulLoc and int(haploStart)<=int(ulStart) and int(haploEnd)>=int(ulEnd)) or \
                          (haploLoc==ulLoc2 and int(haploStart)<=int(ulStart2) and int(haploEnd)>=int(ulEnd2)):               
                         if haploStrain in outList:
                              continue
                         else:
                              # take out underscores and rejoin to test strain names are already in list
                              sampStrainList.append(''.join(haploStrain.split('_')).lower())
                              outList.append(haploType)	
                              outList.append(haploStrain)
                              tempList.append(haploType)
                              tempList.append(haploStrain)

               # for strain in dsanStrains:
               #      found = []
               #      for ele in outList:
               #           if re.search(ele, strain):
               #                found.append(strain)
               #      if found == []:
               #           outList.append('Inbred')
               #           outList.append(strain[5:])
               #           #print('not found: ' + strain)




                         
                     

               # check to see if their were observed dsan strains that have no haplotypes to match to are inbred, map.py inbred dsan issue (missing inbred in hap calls)
               for z in range(0, len(inStrainList)):
                    if inStrainList[z] not in sampStrainList:
                         #print(sampStrainList)
                         #print(inStrainList[z])
                         #print('\n')
                         outList.append('Inbred')
                         outList.append(inStrainList[z])
                    

               # output
               # print(splitVariants)
               # print('\t'.join(outList))
               # print('\n')
               sampledList.append('\t'.join(outList))
          return sampledList
    

     def variantsToFasta(self, inVars, refGenome, outName, flankAmount):
          # output
          outFile = open(outName,'w')
          with gzip.open(refGenome, 'rt') as inFile:
               lines = [line.rstrip() for line in inFile]
          # Using the chromosomes as dict keys
          newDict = {}
          print(len(lines))
          for x in range(0,len(lines)):
               print(lines[x])
               # If you are starting a new chromosome
               if lines[x].startswith(">"):
                    mySplit = lines[x].split()
                    # Add new chroms to dict
                    newSplit = mySplit[0][1:]
                    increment = 0
                    loopDict = {}
                    y = x+1
                    while not lines[y].startswith('>'):
                         loopVal = lines[y]
                         loopDict[increment] = loopVal
                         increment = increment + 80
                         y = y+1
                         # Reaches end of file
                         if y == len(lines)-1:
                              break
                    newDict[newSplit] = loopDict
               
          # max key val for each chrom in dict
          maxBinValPerChrom = {}
          for chrom in newDict:
               chromKeys = newDict[chrom].keys()
               chromKeys.sort()
               maxBinValPerChrom[chrom] = chromKeys[-1]

          for x in range(0, len(inVars)):
               splitLine = inVars[x].split()
               loc1, loc2 = splitLine[0], splitLine[3]
               start1, end1, start2, end2 = int(splitLine[1]), int(splitLine[2]), int(splitLine[4]), int(splitLine[5])
               # The first 3 elements are the location that the fasta referes too, the next 6 refer to the full variant
               name1 = '_'.join(splitLine[:3]) +'|'+ '_'.join(splitLine[:6])
               name2 = '_'.join(splitLine[3:6]) +'|'+ '_'.join(splitLine[:6])

               # regions to capture +/- the flank amount
               startDict1 = ((int(splitLine[1])-1)-flankAmount) - ((int(splitLine[1])-1)-flankAmount) % 80
               endDict1 = ((int(splitLine[2])-1)+flankAmount) - ((int(splitLine[2])-1)+flankAmount) % 80
               posStart1 = ((int(splitLine[1])-1)-flankAmount) % 80
               posEnd1 = ((int(splitLine[2])-1)+flankAmount) % 80
               startDict2 = ((int(splitLine[4])-1)-flankAmount) - ((int(splitLine[4])-1)-flankAmount) % 80
               endDict2 = ((int(splitLine[5])-1)+flankAmount) - ((int(splitLine[5])-1)+flankAmount) % 80
               posStart2 = ((int(splitLine[4])-1)-flankAmount) % 80
               posEnd2 = ((int(splitLine[5])-1)+flankAmount) % 80

               # check both breakpoints for their fasta info
               dicts = {startDict1:[endDict1, name1, loc1, posStart1, posEnd1], startDict2:[endDict2, name2, loc2, posStart2, posEnd2]}

               def checkDicts(dictList):
                    for key in dicts:
                         # loopVars
                         currentStartDict = key
                         currentEndDict = dicts[key][0]
                         currentName = dicts[key][1]
                         currentLoc = dicts[key][2]
                         currentStartPos = dicts[key][3]
                         currentEndPos = dicts[key][4]

                         outFile.write('> ' + currentName+'\n')
                         # if you only need to look at one bin
                         if currentStartDict==currentEndDict:
                              outFile.write(newDict.get(currentLoc)[currentStartDict][currentStartPos:currentEndPos] +'\n')
                         else:
                              # append whole bin to a list until you get to last one
                              blah = []
                              loopCounter = currentEndDict

                              # current chromosomes max value
                              maxVal = maxBinValPerChrom[currentLoc]
                                                           
                              # append whole bins unless on the you need to splice first and last bins 
                              while currentStartDict <= loopCounter:
                                   # end bin 
                                   if currentEndDict==loopCounter:
                                        #print('end')
                                        blah.append(newDict.get(currentLoc)[min(maxVal, currentEndDict)][:currentEndPos])
                                   # start bin
                                   elif currentStartDict==loopCounter:
                                        #print('start')
                                        blah.append(newDict.get(currentLoc)[currentStartDict][currentStartPos:])
                                   # whole bin
                                   else:
                                        #print('whole')
                                        blah.append(newDict.get(currentLoc)[min(maxVal, loopCounter)])
                                   # for next bin on next loop
                                   loopCounter = loopCounter - 80	
                                   # if no bins left because at start of chrom
                                   if loopCounter == -80:
                                        break
                              # append all bin data collected, reverse the list because you added them in backwards????
                              blah.reverse()
                              for y in range(0, len(blah)):
                                   outFile.write(blah[y]+'\n')
               # check fasta for both sides of var
               checkDicts(dicts)
          return


     def subsetFasta(self, inputFile, outputFile, flank):
          import os, sys, subprocess

          # ref genome 
          refGenome = '/scratch/bturne48/BlastDB/dyakRef.fasta'

          # open output
          outFile = open(outputFile, 'w')

          # load into list
          with open(inputFile) as inFile:
               lines = [line.rstrip() for line in inFile]

          # loop over variants
          for var in lines:
               splitLine = var.split()
               chrom1, start1, end1 = splitLine[0], splitLine[1], splitLine[2]
               chrom2, start2, end2 = splitLine[3], splitLine[4], splitLine[5]
               # write seqs to file
               def writeSeq(chrom, start, end):
                    x = subprocess.check_output(['seqkit', 'subseq', '--chr', chrom ,'-r', str(max(1, int(start)-flank))+':'+str(int(end)+flank), refGenome])
                    # format output
                    x = x.rstrip()
                    x = x.split('\n')
                    x[0] = '>'+'_'.join([chrom, start, end]) +'|'+ '_'.join(splitLine[:6])
                    x = '\n'.join(x)
                    outFile.write(x+'\n')
               # for both breakpoints
               writeSeq(chrom1, start1, end1)
               writeSeq(chrom2, start2, end2)
    

     def polarizeRearrangments(self, blastResultsPath, outFileName):
          import os, subprocess
          # dictionary to track variants that were fed to BLAST and match them up
          outFile = open(outFileName, 'w')
          blastReslts= open(blastResultsPath)
          varDict = {}
          for line in blastReslts:
               splitHit = line.split()
               splitVarName = splitHit[0].split('|')
               varSide, fullVarName = splitVarName[0], splitVarName[1]
               hitLoc, hitStart, hitEnd, alnLen, pctIden = splitHit[1].split('-')[0], int(splitHit[8]), int(splitHit[9]), int(splitHit[3]), float(splitHit[2])
               
               # filter based on the alignment length and sequence identity
               if int(alnLen)<=151 or float(pctIden)<95:
                    continue

               # use the fullVarName to create a dict to match both sides of variants together
               if fullVarName not in varDict:
                    #print(hitLoc +'\t'+ str(min(int(hitStart),int(hitEnd))) +'\t'+ str(max(int(hitStart),int(hitEnd))) +'\t'+ str(alnLen) +'\t'+ str(pctIden) +'\t'+ varSide)
                    varDict[fullVarName] = {varSide:[hitLoc +'\t'+ str(min(int(hitStart),int(hitEnd))) +'\t'+ str(max(int(hitStart),int(hitEnd))) +'\t'+ str(alnLen) +'\t'+ str(pctIden) +'\t'+ varSide]}
               else:
                    # if this side of the variant doesn't exist for the full variant dict
                    if varSide not in varDict[fullVarName]:
                         #print(hitLoc +'\t'+ str(min(int(hitStart),int(hitEnd))) +'\t'+ str(max(int(hitStart),int(hitEnd))) +'\t'+ str(alnLen) +'\t'+ str(pctIden) +'\t'+ varSide)
                         varDict[fullVarName].update({varSide:[hitLoc +'\t'+ str(min(int(hitStart),int(hitEnd))) +'\t'+ str(max(int(hitStart),int(hitEnd))) +'\t'+ str(alnLen) +'\t'+ str(pctIden) +'\t'+ varSide]})
                    # if this half of the variant already exits filter out similar results to save memory
                    else:    
                         varDict[fullVarName][varSide].append(hitLoc +'\t'+ str(min(int(hitStart),int(hitEnd))) +'\t'+ str(max(int(hitStart),int(hitEnd))) +'\t'+ str(alnLen) +'\t'+ str(pctIden) +'\t'+ varSide)
                         #print(hitLoc +'\t'+ str(min(int(hitStart),int(hitEnd))) +'\t'+ str(max(int(hitStart),int(hitEnd))) +'\t'+ str(alnLen) +'\t'+ str(pctIden) +'\t'+ varSide)
                         #print('\n')

          # for each full variant in dictionary
          polarList = []
          halfPolarList = {}
          for key in varDict:
               #print(key, varDict[key])
               # check each key/value to see if it found both sides of the variant
               if len(varDict[key]) >= 2:
                    polarizeVar = False
                    # if both sides found, create list of positions to see if the sides come within 1kb of each other
                    for side1 in varDict[key]:
                         # positions list for side 1
                         for posList1 in varDict[key][side1]:
                              # if you already found a match for the variant and will polarize, time to move to the next variant
                              if polarizeVar == True:
                                   break
                              splitPos1 = posList1.split('\t')
                              sideLoc1, startPos1, endPos1, alnLen1, pctIden1 = splitPos1[0], int(splitPos1[1]), int(splitPos1[2]), int(splitPos1[3]), float(splitPos1[4])
                              # put this half of the variant into a dict to track confi rate
                              halfPolarList[key]=[key, splitPos1[5]]
                              for side2 in varDict[key]:
                                   # make sure you're looking at different sides, idk why i have to set it up like this, probs cause im dumb
                                   if side1!=side2:
                                        # positions list for side 2
                                        for posList2 in varDict[key][side2]:
                                             splitPos2 = posList2.split('\t')
                                             sideLoc2, startPos2, endPos2, alnLen2, pctIden2 = splitPos2[0], int(splitPos2[1]), int(splitPos2[2]), int(splitPos2[3]), float(splitPos2[4])
                                             
                                             # put other half of the variant into a dict to track confi rate
                                             if splitPos2[5] not in halfPolarList[key]:
                                                  halfPolarList[key] = halfPolarList[key] + [splitPos2[5]]

                                             # make sure % overlap is ok
                                             def getOverlap(a, b):
                                                  return max(0, min(a[1], b[1]) - max(a[0], b[0]))
                                             if getOverlap([startPos1, endPos1], [startPos2, endPos2]) != 0:
                                                  case = 0
                                                  # is it overlap case 1 or 2
                                                  if endPos1 > startPos2:
                                                       case = 1
                                                  elif endPos2 > startPos1:
                                                       case = 2
                                                  # use case to calc overlap length
                                                  if case == 1:
                                                       overlapLength = endPos1 - startPos2
                                                  elif case == 2:
                                                       overlapLength = endPos2 - startPos1
                                                  # og seq length
                                                  splitVarSide = varSide.split("-")
                                                  varSideStart, varSideEnd = int(splitVarSide[1]), int(splitVarSide[2])
                                                  ogSeqLength = varSideEnd - varSideStart
                                                  # snp?
                                                  if ogSeqLength == 0:
                                                       continue
                                                  # if more than 10% overlap of original sequence length
                                                  if float(overlapLength)/float(ogSeqLength) > 0.1:
                                                       # print(posList1)
                                                       # print(posList2)
                                                       # print(overlapLength)
                                                       continue

                                             # check to see if within 1kb of eachother in BLAST result and mapped to same chrom
                                             if sideLoc1==sideLoc2 and (abs(min(startPos1, endPos1)-min(startPos2, endPos2)) <=1000 or abs(max(startPos1, endPos1)-max(startPos2, endPos2))<=1000):
                                                  #print(posList1)
                                                  #print(posList2)
                                                  #print('\n')
                                                  #if getOverlap([startPos1, endPos1], [startPos2, endPos2]) != 0:
                                                       #print(posList1)
                                                       #print(posList2)
                                                       #print('overlap:'+str(getOverlap([startPos1, endPos1], [startPos2, endPos2])))
                                                  # print(min(startPos1, endPos1), min(startPos2, endPos2),max(startPos1, endPos1), max(startPos2, endPos2))
                                                  # print(abs(min(startPos1, endPos1)-min(startPos2, endPos2)), abs(max(startPos1, endPos1)-max(startPos2, endPos2)<=1000))
                                                  #print('\n')
                                                  polarizeVar = True
                                                  polarList.append('\t'.join(key.split('-'))+'\t'+str(alnLen1)+'\t'+str(pctIden1)+'\t'+str(alnLen2)+'\t'+str(pctIden2))   
                                                  #print('\t'.join(key.split('_'))+'\t'+str(alnLen1)+'\t'+str(pctIden1)+'\t'+str(alnLen2)+'\t'+str(pctIden2)+'\t'+'\n')
                                                  outFile.write('\t'.join(key.split('-'))+'\t'+str(alnLen1)+'\t'+str(pctIden1)+'\t'+str(alnLen2)+'\t'+str(pctIden2)+'\t'+'\n')
                                                  #print('\n')
                                                  break  


          return polarList
          # poo = []
          # for x in halfPolarList:
          #      if x not in polarList:
          #           poo.append
          # #print(len(poo))
          # return polarList        


     def calcSFS(self, clusteredVars, sampledHaps, targetPop, polarVars='', shouldProject=False, maxno=0, repbaseFilterTE=False, repbaseBlastResults=[], makeRepbaseSFS=False, hetCorr=''):
          # Input:
          # 1) Clustered variant object
          # 2) Sampled haplotypes object or a file that has the list
          # 3) Target population: Leave blank if no subpop desired
          # 4) Will you be polarizing the mutations using blast results or a list of variants
          # 5) List of variants or list of blast results used to polarize
          # 6) Project down: If true uses the hypergeometric distribution to project SFS to a smaller population size. 
          # 7) Max number in population: Number desired if shouldProject is true
          #    If projection is true then the function will return a list of the proportion of the population per maxno
          #    EX) maxno=15; n=1:0.4, n=2:0.3, ...... , n=15:0.03   
          #    
          #    IF projection is false then the function ouputs a list with all variants and there frequency 

          from collections import Counter
          from scipy import special
          import re
          import os
          #import math

          # * START: Make sure input parameters are correct types * #

          # * END: Make sure input parameters are correct types * #

          # * START: Innit projection code * # 
          # Set up to use the geometric function later
          def choose(n,k):
               return special.binom(n, k)
          # Some smart shit that rebekahs code does 
          SFSVec=[0]*(maxno+1)
          # * END: Innit projection code * # 

          # # read in het correction stuff
          hetCorrDict = {}
          if hetCorr != '':
               with open(hetCorr) as hetCorrectionFile:
                    for line in hetCorrectionFile:
                         name = '\t'.join(line.split()[:6])
                         if name not in hetCorrDict:
                              hetCorrDict[name] = line.split()[-1]

          
          # * START: Read in previously created sampled strains or a list from sv class * # 
          sampledDict = {}
          if type(sampledHaps) == str:
               with open(sampledHaps) as inFile:
                    sampled = [line.rstrip() for line in inFile]
          else:
               sampled = sampledHaps

               
          for x in range(0,len(sampled)):
               splitSamp = sampled[x].split('\t')
               sampledDict['\t'.join(splitSamp[0:6])] = '\t'.join(splitSamp[6:])
          # * END: Read in previously created sampled strains or a list from sv class * # 

          # * START: Fix None Value * #
          # If there is no value for a haplotype when you use the Counter function later on, I need to fix a None value
          def fixNone(val):
               if val==None:
                    return 0
               else:
                    return val
          # * END: Fix None Value * #

          # *START* Polarize Variants
          if polarVars != '':
               with open(polarVars) as inPolar:
                    polarVarList = ['\t'.join(line.rstrip().split('\t')[:6]) for line in inPolar]
          else:
               polarVarList = []
          
          # * START: Process TEs from Repbase' * #
          repbaseFilterDict = {}
          if repbaseFilterTE==True and type(repbaseBlastResults)==str:
               # loop over blast results to find matches with variants
               inTE = open(repbaseBlastResults)
               for line in inTE:
                    #print(line)
                    splitRepLine = line.split('\t')
                    splitLoc = '\t'.join((splitRepLine[0].split('-')))
                    # Store which half of the variant is referenced in the blast result
                    refVar = splitLoc.split('|')[0]
                    # split away the actual result of the blast
                    splitLoc = splitLoc.split('\t')
                    splitLoc = '\t'.join([splitLoc[2].split('|')[1]] + splitLoc[3:8])
                    # check if both sides of variant are a blast hit
                    if splitLoc in repbaseFilterDict:
                         alreadyRef = repbaseFilterDict.get(splitLoc)
                         # If this half of the variant is already stored
                         if refVar in alreadyRef:
                              continue
                         else:
                              repbaseFilterDict[splitLoc].append(refVar)
                    else:
                         # If the variant is new and not in the filter dictionary
                         repbaseFilterDict[splitLoc] = [refVar]
          ## Check number of blast results that hit 0,1,2 sides of a variant, do this by dividing by 3
          for dab in repbaseFilterDict:
               #print(dab)
               print(dab+'\t'+'\t'.join(repbaseFilterDict.get(dab)))
               #print('\n')
          # # #* END: Process TEs from Repbase' * #
        
          # If subpop target desired, Regex search what strains you want
          if type(targetPop) == str:


               outList = []
               # I can do samp in sampled because checkSampled method uses clustered vars as a parameter
               for samp in sampled:
                    # new for tes
                    splitSamp = samp.split('\t')
                    #splitSamp[0], splitSamp[3] = splitSamp[0].split('_')[0], splitSamp[3].split('_')[0]
                    #splitSamp[0], splitSamp[3] = splitSamp[0].split('_')[0], splitSamp[3].split('_')[0]                    
                    teName = '\t'.join(splitSamp[:6])    
                    # back to normal
                    #splitSamp = samp.split('\t')
                    varName = '\t'.join(splitSamp[:6])  
                    #print(varName) 
                    #print(clusteredVars.get(varName))              
                    splitObsv = clusteredVars.get(varName).split('\t')
                    # print(teName)
                    # print(varName)
                    # print('doo')
                    #print('\t'.join(splitObsv) +'\t'+str(len(splitObsv)))
                    #print(repbaseFilterDict)
                    # Filtering TEs out 
                    if teName in repbaseFilterDict and makeRepbaseSFS==False and repbaseFilterTE==True:
                         continue
                    # Looking for just TEs
                    if makeRepbaseSFS==True and repbaseFilterTE==True:
                         if teName in repbaseFilterDict:
                              # print(teName)
                              # print(varName)
                              # print('\n')
                              splitObsv = clusteredVars.get(varName).split('\t')
                         else:    
                              continue

                    # Check to see if there were any observed in target strain, if not skip
                    hapListObsv = []
                    doneList = []
                    for x in range(0, len(splitObsv),2):
                         if re.search(targetPop, splitObsv[x+1]):
                              # add in het correction stuff here, change from inbred to pred_homo which will be 2/x instead of 1/x
                              if varName in hetCorrDict and hetCorrDict[varName] == ''.join(splitObsv[x+1].split('_')).lower() and splitObsv[x]=='hetero':
                                   #print(varName)
                                   hapListObsv.append('Pred_homo')
                              if splitObsv[x]+splitObsv[x+1] not in doneList:
                                   hapListObsv.append(splitObsv[x])
                                   doneList.append(splitObsv[x]+splitObsv[x+1])
                    # Count sampled haps based on strain 
                    hapListSamp = []
                    for x in range(6, len(splitSamp),2):
                         #print(splitSamp)
                         if re.search(targetPop, splitSamp[x+1]):
                              hapListSamp.append(splitSamp[x])
                      
                    print(hapListObsv)
                    print(hapListSamp)

                    # Count observed and sampled haplotypes
                    countHapsObsv, countHapsSamp = Counter(hapListObsv), Counter(hapListSamp)
                    inbObsv, hetObsv = fixNone(countHapsObsv.get('Inbred')), fixNone(countHapsObsv.get('hetero'))
                    # het correction
                    predHomo = fixNone(countHapsObsv.get('Pred_homo'))
                    inbSamp, hetSamp = fixNone(countHapsSamp.get('Inbred')), fixNone(countHapsSamp.get('hetero'))
                    
                    numerator, denominator = float(inbObsv) + (float(hetObsv)) + (float(predHomo)*2), float(inbSamp) + (2*float(hetSamp))

                    print(numerator, denominator, inbObsv, hetObsv)

                    # #NOTE: for when making an SFS without the 0's
                    # if numerator == 0:
                    #     continue

                    # allele freq for the subpop
                    if denominator==0:
                         sfsSubPop = 0
                    else:
                         sfsSubPop = float(numerator/denominator)

                    #print('sfssub '+str(sfsSubPop))

                    #Polarize mutations if needed
                    if '\t'.join(splitSamp[:6]) in polarVarList:
                         numerator, sfsSubPop = float(denominator - numerator), float(1 - sfsSubPop)

                    #print('sfssub2 '+str(sfsSubPop))


                    if sfsSubPop > 1:
                         print(numerator, sfsSubPop)
                         print('ERROR')
                         print('\n')
                         continue

                    #print('sfssub3 '+str(sfsSubPop))

                    # * START: Should Project * #
                    # #If you want to project the SFS down to a smaller pop size
                    if shouldProject == True:
                         # Update sfs value if denominator is greater than specifed projection sample size
                         #print('denmom'+str(denominator))
                         if denominator >=maxno:
                              for i in range(0, maxno+1):
                                   SFSVec[i] = SFSVec[i] + choose(numerator, i) * choose(denominator-numerator, maxno-i)/ (choose(denominator, maxno))
                              #print(SFSVec)
                    # * END: Should Project * #
                    else:
                         #outList.append(varName +'\t'+ str(numerator) +'\t'+ str(denominator) +'\t'+ str(sfsSubPop))
                         # NOTE: if filtering out small sample size rearrangements
                         if denominator < 15:
                              outList.append(varName +'\t'+ str(numerator) +'\t'+ str(denominator) +'\t'+ 'NA')
                         else:
                              outList.append(varName +'\t'+ str(numerator) +'\t'+ str(denominator) +'\t'+ str(sfsSubPop))

          # * END: SFS Calc and Projection
          if shouldProject == True:
               projList = []
               #print(SFSVec)
               for item in SFSVec:
                    projList.append(str(item/sum(SFSVec)))
                    #print(item/sum(SFSVec))
               return projList
          else:
               return outList


     def matchTE(self, clusteredVars, blastResults, outFile):
          writeToOut = open(outFile, 'w')
          # * START: Process TEs from Repbase' * #
          repbaseFilterDict = {}
          # loop over blast results to find matches with variants
          inTE = open(blastResults)
          for line in inTE:
               splitRepLine = line.split('\t')
               splitLoc = '\t'.join((splitRepLine[0].split('-')))
               # Store which half of the variant is referenced in the blast result
               refVar = splitLoc.split('|')[0]
               # split away the actual result of the blast
               splitLoc = splitLoc.split('\t')
               #print(splitLoc)
               splitLoc = '\t'.join([splitLoc[2].split('|')[1]] + splitLoc[3:8])
               # if float(splitRepLine[2])<95 or int(splitRepLine[3])<151:
               #      continue
               # check if both sides of variant are a blast hit
               if splitLoc in repbaseFilterDict:
                    alreadyRef = repbaseFilterDict.get(splitLoc)
                    # If this half of the variant is already stored
                    if refVar in alreadyRef:
                         continue
                    else:
                         repbaseFilterDict[splitLoc].append(refVar)
               else:
                    # If the variant is new and not in the filter dictionary
                    repbaseFilterDict[splitLoc] = [refVar]

               
          ## Check number of blast results that hit 0,1,2 sides of a variant, do this by dividing by 3
          for dab in repbaseFilterDict:
               writeToOut.write(dab +'\t'+ '\t'.join(repbaseFilterDict.get(dab)) +'\n')
          #* END: Process TEs from Repbase' * #
          return 


     def calcFST(self, clusteredVars, sampledHaps, subPopRegexList, subPopNamesList, compareStrains):
          import re
          from collections import Counter
     
          # * START: Read in previously created sampled strains or a list from sv class * # 
          sampledDict = {}
          if type(sampledHaps) == str:
               with open(sampledHaps) as inFile:
                    sampled = [line.rstrip() for line in inFile]
          else:
               sampled = sampledHaps
          for x in range(0,len(sampled)):
               splitSamp = sampled[x].split('\t')
               sampledDict['\t'.join(splitSamp[0:6])] = '\t'.join(splitSamp[6:])
          # * END: Read in previously created sampled strains or a list from sv class * # 

          # * START: Fix None Value * #
          # If there is no value for a haplotype when you use the Counter function later on, I need to fix a None value
          def fixNone(val):
               if val==None:
                    return 0
               else:
                    return val
          # * END: Fix None Value * #
          
          # init output dict
          varDict = {}
          for subPop in subPopNamesList:
               varDict[subPop] = []

          # read in het correction stuff
          # hetCorrDict = {}
          # with open('/users/bturne48/RogersLab/RemoteEdit/ParseIlluminaData/Data/OutputFromScripts/hetCorrFinal.txt') as hetCorrectionFile:
          #      for line in hetCorrectionFile:
          #           #print(name)
          #           name = '\t'.join(line.split()[:6])
          #           if name not in hetCorrDict:
          #                hetCorrDict[name] = line.split()[-1]
          
          # calc fst for each observed variant
          for samp in sampled:
                    splitSamp = samp.split('\t')
                    varName = '\t'.join(splitSamp[:6])
                    splitObsv = clusteredVars.get(varName).split('\t') 
                    # Check to see if there were any observed in current target population, if not skip and try next
                    # reset strain dict for each variant
                    strainsDict = {}
                    subPopListIndex = 0
                    for subPop in subPopRegexList:
                         targetPop = subPop
                         hapListObsv = []
                         strains = []
                         # use regex list to search subpops
                         for x in range(0, len(splitObsv),2):
                              if re.search(targetPop, splitObsv[x+1]):
                                   # add in het correction stuff here
                                   if varName in hetCorrDict and hetCorrDict[varName] == splitObsv[x+1]:
                                        strains.append(splitObsv[x+1])
                                        hapListObsv.append('Pred_homo')
                                   else:
                                        strains.append(splitObsv[x+1])
                                        hapListObsv.append(splitObsv[x])

                         # Count sampled haps based on strain 
                         hapListSamp = []
                         for x in range(7, len(splitSamp),2):
                              if re.search(targetPop, splitSamp[x]):
                                   hapListSamp.append(splitSamp[x-1])
                              
                         # allele freqs
                         countHapsObsv = Counter(hapListObsv)
                         countHapsSamp = Counter(hapListSamp)
                         inbObsv = fixNone(countHapsObsv.get('Inbred'))
                         hetObsv = fixNone(countHapsObsv.get('hetero'))
                         inbSamp = fixNone(countHapsSamp.get('Inbred'))
                         hetSamp = fixNone(countHapsSamp.get('hetero'))
                         predHomo = fixNone(countHapsObsv.get('Pred_homo'))

                         # variables for fst
                         #print(countHapsObsv)
                         subPopObserved, subPopSampled = float(inbObsv) + float(hetObsv) + float(predHomo), (float(inbSamp) + 2*float(hetSamp))
                         #print(subPop, subPopObserved, subPopSampled)
                         strainsDict[subPopNamesList[subPopListIndex]] = [subPopObserved, subPopSampled]

                         # move to next name in subpop list for output
                         subPopListIndex += 1

                    #if strainsDict.get(compareStrains[0]) != None and strainsDict.get(compareStrains[1]) != None:
                    #     continue

                    # if there wasnt mutation data for a population that we want to comapre
                    if strainsDict.get(compareStrains[0]) == None or strainsDict.get(compareStrains[1]) == None:
                         # if there wasnt mutation data for either we continue to the next variant
                         if strainsDict.get(compareStrains[0]) == None and strainsDict.get(compareStrains[1]) == None:
                              continue
                         # if there was mutation data for one of them we want to be able to calculate it later
                         if strainsDict.get(compareStrains[0]) == None:
                              strainsDict[compareStrains[0]] = [0,0]
                         if strainsDict.get(compareStrains[1]) == None:
                              strainsDict[compareStrains[1]] = [0,0]


                              
                    # FST calculation here
                    pop1 = strainsDict.get(compareStrains[0])
                    pop2 = strainsDict.get(compareStrains[1])
                    # if pop1 has no sampled (denom) data fix the division by 0 error
                    if pop1[1] == 0:
                         p_1, q_1 = 0, 1
                         c_1, c_2 = 0, 1
                    else:
                         p_1, q_1= pop1[0]/pop1[1], 1 - (pop1[0]/pop1[1])
                         c_1, c_2 = pop1[1]/(pop1[1]+pop2[1]), pop2[1]/(pop1[1]+pop2[1])
                    # same error, pop2
                    if pop2[1] == 0:
                         p_2, q_2 = 0, 1
                         c_1, c_2 = 0, 1
                    else:
                         p_2, q_2 = pop2[0]/pop2[1], 1 - (pop2[0]/pop2[1])
                         c_1, c_2 = pop1[1]/(pop1[1]+pop2[1]), pop2[1]/(pop1[1]+pop2[1])
                    # same error but for both
                    if pop1[1]==0 and pop2[1]==0:
                         p_bar, q_bar = 0, 1                  
                    else:
                         p_bar, q_bar = (pop1[0] + pop2[0])/(pop1[1]+pop2[1]), 1 - ((pop1[0] + pop2[0])/(pop1[1]+pop2[1]))
                    # print(splitObsv)
                    # print(splitSamp)
                    # print("P_bar: \n" + str(p_bar))
                    # print(compareStrains[0]+": ")
                    # print(p_1, q_1, pop1, c_1)
                    # print(compareStrains[1]+": ")
                    # print(p_2, q_2, pop2, c_2)
                    #print('\n')
                    # if p_bar is 0 or 1 you can't divide by 0 so set to 0
                    if p_bar==0.0 or p_bar==1.0:
                         fst = 0
                    else:
                         fst = ((p_bar * q_bar) - ((c_1 * p_1 * q_1) + (c_2 * p_2 * q_2)))/(p_bar * q_bar)
                    print(varName +'\t'+ str(fst))# +'\t'+ '\t'.join(strains))
                    #print('\n')
                    # print('FST: ' + str(fst))
                    # print('\n')
                    #varDict[]

          return varDict


     def popgenMap(self, outFileName, mutFile={}, yakubaInput=True, greyInput=True, objInput=False):
          # Parse whether it's an object input from previous stuff or if it's a file
          if objInput==True:
               # for now it's a dict weird thing that i need to fix
               variants = mutFile
               if yakubaInput==False:
                    variants = mutFile.get("Dsan")
                    with open("/users/bturne48/RogersLab/RemoteEdit/ParseIlluminaData/Data/MiscInput/Pi.sant") as robFilesDsan:
                         dsanTajD = [line.rstrip() for line in robFilesDsan]
                    # create output file
                    outRegSant = open(outFileName,'w')
                    dsanWindowDict = {'2L':[], '2R':[],'3L':[],'3R':[],'X':[]}
                    for x in range(0, len(dsanTajD)):
                         splitLine = dsanTajD[x].split()
                         frameLoc = splitLine[0]
                         for chrom in dsanWindowDict:
                              if frameLoc == chrom:
                                   currentList = dsanWindowDict.get(chrom)
                                   currentList.append(dsanTajD[x])
                    
               else:
                    if greyInput == True:
                         variants = mutFile.get("CY/NY/Grey")
                         with open("/users/bturne48/RogersLab/RemoteEdit/ParseIlluminaData/Data/MiscInput/Pi.grey") as robFilesDsan:
                              dyakTajD = [line.rstrip() for line in robFilesDsan]                      
                         outRegYak = open(outFileName,'w')
                         dyakWindowDict = {'2L':[], '2R':[],'3L':[],'3R':[],'X':[]}
                         for x in range(0, len(dyakTajD)):
                              splitLine = dyakTajD[x].split()
                              frameLoc = splitLine[0]
                              for chrom in dyakWindowDict:
                                   if frameLoc == chrom:
                                        currentList = dyakWindowDict.get(chrom)
                                        currentList.append(dyakTajD[x])
                    else:
                         variants = mutFile.get("Oran") 
                         with open("/users/bturne48/RogersLab/RemoteEdit/ParseIlluminaData/Data/MiscInput/Pi.yak") as robFilesDsan:
                              dyakTajD = [line.rstrip() for line in robFilesDsan]                      
                         outRegYak = open(outFileName,'w')
                         dyakWindowDict = {'2L':[], '2R':[],'3L':[],'3R':[],'X':[]}
                         for x in range(0, len(dyakTajD)):
                              splitLine = dyakTajD[x].split()
                              frameLoc = splitLine[0]
                              for chrom in dyakWindowDict:
                                   if frameLoc == chrom:
                                        currentList = dyakWindowDict.get(chrom)
                                        currentList.append(dyakTajD[x])
                    

          # if muts from a file
          else:
               # if input mut is dsan file 
               if yakubaInput==False:
                    # open seperated variants
                    with open('/nobackup/rogers_research/Brandon/PipelineResults/MyOgCluster/dsan_grey_deltap_plot.txt') as mutFile:
                         variants = [line.rstrip() for line in mutFile]
                    # open correct popgen file
                    with open("/users/bturne48/RogersLab/RemoteEdit/ParseIlluminaData/Data/MiscInput/Pi.sant") as robFilesDsan:
                         dsanTajD = [line.rstrip() for line in robFilesDsan]
                    # create output file
                    outRegSant = open(outFileName,'w')
                    # read in 10kb windows to dict
                    dsanWindowDict = {'2L':[], '2R':[],'3L':[],'3R':[],'X':[]}
                    for x in range(0, len(dsanTajD)):
                         splitLine = dsanTajD[x].split()
                         frameLoc = splitLine[0]
                         for chrom in dsanWindowDict:
                              if frameLoc == chrom:
                                   currentList = dsanWindowDict.get(chrom)
                                   currentList.append(dsanTajD[x])
               # if input mut is dyak file
               else:
                    if greyInput==True:
                         with open('/nobackup/rogers_research/Brandon/PipelineResults/MyOgCluster/dsan_grey_deltap_plot.txt') as mutFile:
                              variants = [line.rstrip() for line in mutFile]
                         with open("/users/bturne48/RogersLab/RemoteEdit/ParseIlluminaData/Data/MiscInput/Pi.grey") as robFilesDsan:
                              dyakTajD = [line.rstrip() for line in robFilesDsan]
                         outRegYak = open(outFileName,'w')
                         dyakWindowDict = {'2L':[], '2R':[],'3L':[],'3R':[],'X':[]}
                         for x in range(0, len(dyakTajD)):
                              splitLine = dyakTajD[x].split()
                              frameLoc = splitLine[0]
                              for chrom in dyakWindowDict:
                                   if frameLoc == chrom:
                                        currentList = dyakWindowDict.get(chrom)
                                        currentList.append(dyakTajD[x])
                    else:
                         with open('/nobackup/rogers_research/Brandon/PipelineResults/MyOgCluster/dsan_oran_deltap_plot.txt') as mutFile:
                              variants = [line.rstrip() for line in mutFile]
                         with open("/users/bturne48/RogersLab/RemoteEdit/ParseIlluminaData/Data/MiscInput/Pi.yak") as robFilesDsan:
                              dyakTajD = [line.rstrip() for line in robFilesDsan]
                         outRegYak = open(outFileName,'w')
                         dyakWindowDict = {'2L':[], '2R':[],'3L':[],'3R':[],'X':[]}
                         for x in range(0, len(dyakTajD)):
                              splitLine = dyakTajD[x].split()
                              frameLoc = splitLine[0]
                              for chrom in dyakWindowDict:
                                   if frameLoc == chrom:
                                        currentList = dyakWindowDict.get(chrom)
                                        currentList.append(dyakTajD[x])

          # parse the 10kb windows and associate the variant with the window closest to the midpoint of the variant/window
          for x in range(0, len(variants)):
               splitVar = variants[x].split()
               varLoc1 = splitVar[0]
               varStart1 = int(splitVar[1])
               varEnd1 = int(splitVar[2])
               mid1 = (varEnd1+varStart1)/2.0
               varLoc2 = splitVar[3]
               varStart2 = int(splitVar[4])
               varEnd2 = int(splitVar[5])
               mid2 = (varEnd2+varStart2)/2.0
               myList1_Dsan = []
               myList2_Dsan = []
               myList1_Dyak = []
               myList2_Dyak = []

               #print(varLoc1, mid1, varLoc2, mid2)

               if yakubaInput==False:
                    for key in dsanWindowDict:
                         if varLoc1 == key:
                              chromList = dsanWindowDict.get(key)
                              if varStart1 <= 5500:
                                   for y in range(0,6):
                                        splitWindow = chromList[y].split()
                                        window = int(splitWindow[1])
                                        if abs(varStart1-window) <= 500:
                                             myList1_Dsan.append(varLoc1 +' '+ str(varStart1) +' '+ str(varEnd1) + ' ' + ' '.join((chromList[y].split()[2:])))
                                             break
                              else:
                                   for y in range(0,len(chromList)):
                                        splitWindow = chromList[y].split()
                                        windowMid = int(splitWindow[1])+500
                                        if abs(windowMid - mid1) <= 500:
                                             myList1_Dsan.append(varLoc1 +' '+ str(varStart1) +' '+ str(varEnd1) + ' ' + ' '.join((chromList[y].split()[2:])))
                                             break
                         if varLoc2 == key:
                              chromList = dsanWindowDict.get(key)
                              if varStart2 <= 5500:
                                   for y in range(0,6):
                                        splitWindow = chromList[y].split()
                                        window = int(splitWindow[1])
                                        if abs(varStart2-window) <= 500:
                                             myList2_Dsan.append(varLoc2 +' '+ str(varStart2) +' '+ str(varEnd2) + ' ' + ' '.join((chromList[y].split()[2:])))
                                             break
                              else:
                                   for y in range(0,len(chromList)):
                                        splitWindow = chromList[y].split()

                                        windowMid = int(splitWindow[1])+500
                                        if abs(windowMid - mid2) <= 500:
                                             #print(splitWindow)
                                             myList2_Dsan.append(varLoc2 +' '+ str(varStart2) +' '+ str(varEnd2) + ' ' + ' '.join((chromList[y].split()[2:])))
                                             break
               else:
                    for key in dyakWindowDict:
                         if varLoc1 == key:
                              chromList = dyakWindowDict.get(key)
                              if varStart1 <= 5500:
                                   for y in range(0,6):
                                        splitWindow = chromList[y].split()
                                        window = int(splitWindow[1])
                                        if abs(varStart1-window) <= 500:
                                             myList1_Dyak.append(varLoc1 +' '+ str(varStart1) +' '+ str(varEnd1) + ' ' +' '.join((chromList[y].split()[2:])))
                                             break
                              else:
                                   for y in range(0,len(chromList)):
                                        splitWindow = chromList[y].split()
                                        windowMid = int(splitWindow[1])+500
                                        if abs(windowMid - mid1) <= 500:
                                             myList1_Dyak.append(varLoc1 +' '+ str(varStart1) +' '+ str(varEnd1) + ' ' + ' '.join((chromList[y].split()[2:])))
                                             break
                         if varLoc2 == key:
                              chromList = dyakWindowDict.get(key)
                              if varStart2 <= 5500:
                                   for y in range(0,6):
                                        splitWindow = chromList[y].split()
                                        window = int(splitWindow[1])
                                        if abs(varStart2-window) <= 500:
                                             myList2_Dyak.append(varLoc2 +' '+ str(varStart2) +' '+ str(varEnd2) + ' '+ ' '.join((chromList[y].split()[2:])))
                                             break
                              else:
                                   for y in range(0,len(chromList)):
                                        splitWindow = chromList[y].split()
                                        windowMid = int(splitWindow[1])+500
                                        if abs(windowMid - mid2) <= 500:
                                             myList2_Dyak.append(varLoc2 +' '+ str(varStart2) +' '+ str(varEnd2) +' '+ ' '.join((chromList[y].split()[2:])))
                                             break

               # used to produce file that has missing lines removed	
               if yakubaInput==False:
                    #print(' '.join(myList1_Dsan) +' '+ ' '.join(splitVar))
                    #print(' '.join(myList2_Dsan) +' '+ ' '.join(splitVar))
                    outRegSant.write(' '.join(myList1_Dsan) +' '+ ' '.join(splitVar)+'\n')
                    outRegSant.write(' '.join(myList2_Dsan) +' '+ ' '.join(splitVar)+'\n')
               else:
                    # checks to see if their are regions missing
                    if myList1_Dyak != []:
                         #print(' '.join(myList1_Dyak) +' '+ ' '.join(splitVar))
                         outRegYak.write(' '.join(myList1_Dyak) +' '+ ' '.join(splitVar)+'\n')
                    if myList2_Dyak != []:
                         outRegYak.write(' '.join(myList2_Dyak) +' '+ ' '.join(splitVar)+'\n')
                         #print(' '.join(myList2_Dyak) +' '+ ' '.join(splitVar))
               #print('\n')


     def statSigVars(self, statFile, statIndexTarget, pval=0.05, mathstuff=False, detail=False):
          ### * Returns a list of lists with lower bound and upper bound * ### 

          import sys, os, subprocess, numpy as np
          # created a sorted list based off your desired stat 
          popgenList = subprocess.check_output(["sort", "-k"+str(statIndexTarget+1), "-n", statFile])
          # list format and take out the trailing newline
          popgenList = popgenList.split('\n')[:-1]

          # divide pval by 2
          pval = pval/2.0

          if mathstuff==False:
               # calculate the tails for each side 
               lowerBoundIndex, uppperBoundIndex = int(len(popgenList) * pval), int(len(popgenList) * (1-pval))
               lowerSigLines, upperSigLines = popgenList[:lowerBoundIndex], popgenList[uppperBoundIndex:]
               return[lowerSigLines, upperSigLines]
          else:
               # create stat list with input param
               statList = []
               for line in popgenList:
                    statList.append(float(line.split()[statIndexTarget])) # maybe add a try/value error statment here
                    #print(float(line.split()[statIndexTarget]))
               # take the s.d of variants list and calculate anything more than 2 sd away
               def mean(data):
                    """Return the sample arithmetic mean of data."""
                    n = len(data)
                    if n < 1:
                         raise ValueError('mean requires at least one data point')
                    return sum(data)/n # in Python 2 use sum(data)/float(n)

               def _ss(data):
                    """Return sum of square deviations of sequence data."""
                    c = mean(data)
                    ss = sum((x-c)**2 for x in data)
                    return ss               

               def stddev(data, ddof=0):
                    """Calculates the population standard deviation
                    by default; specify ddof=1 to compute the sample
                    standard deviation."""
                    n = len(data)
                    if n < 2:
                         raise ValueError('variance requires at least two data points')
                    ss = _ss(data)
                    pvar = ss/(n-ddof)
                    return pvar**0.5
               # currently takes 2sd away from mean
               lowSdBound, highSdBound = mean(statList) - 2*stddev(statList), mean(statList) + 2*stddev(statList)
               
               if detail == True:
                    print("Summary Stats for column " + str(statIndexTarget+1) +":\nMean: " + str(mean(statList)) + "\tStdev: " + str(stddev(statList)))
               for line in popgenList:
                    if float(line.split()[statIndexTarget]) < lowSdBound or float(line.split()[statIndexTarget]) > highSdBound:
                         print(line)
                         continue
     
     
     def intergenicRegions(self, geneFile, snpFile, chromList, outFile):
          # parse blast snp results?
          
          # read in files
          with open(geneFile) as inFile:
               genes = [line.rstrip() for line in inFile]
          
          # filter out alt splicing to most exon juncts
          geDict, goodChroms = {}, chromList
          for gene in genes:
               ge, chrom = gene.split()[0].split('-')[0], gene.split()[2].split('=')[1]
               # skip the shit if it isnt a chrom we want
               if chrom not in goodChroms:
                    continue
               # save alt splice lines to a dict
               if ge not in geDict:
                    geDict[ge] = [gene]
               else:
                    geDict[ge].append(gene)

          finalGeDict = {}
          for geId in geDict:
               # count number of exon juncts, if there are more replace the value
               for line in geDict[geId]:
                    if geId not in finalGeDict:
                         # save how many exon juncts there are -1 because there is an extra '|'
                         finalGeDict[geId] = [line, len(line.split()[5].split('|'))-1]
                    # if there are the most exon juncts
                    else:
                         if len(line.split()[5].split('|'))-1 > finalGeDict[geId][1]:
                              finalGeDict[geId] = [line, len(line.split()[5].split('|'))-1]

          # reformat so code below can use
          genesFilter = []
          for key in finalGeDict:
               genesFilter.append(finalGeDict[key][0])

          # create a sorted list so you can find intron regions
          normPosList = {'2L':[], '2R':[], '3L':[], '3R':[], 'X':[]}
          for gene in genesFilter:
               splitLine = gene.split()
               geNum, chrom, strand, exonJuncts = splitLine[0], splitLine[2].split('=')[1], int(splitLine[4]), splitLine[5]
               targetDict = normPosList[chrom]
               fullGeneRange = [exonJuncts.split("|")[0].split('..')[0], exonJuncts.split("|")[-2].split('..')[1]]
               
               # print genes?
               #print(geNum +'\t'+ chrom +'\t'+ '\t'.join(fullGeneRange))
               # if reverse strand
               if strand == -1:
                    # string split the last exon junct for the reverse strand
                    targetJunctStart = int(exonJuncts.split("|")[-2].split('..')[0])
                    regionStart, regionStop = targetJunctStart-8, targetJunctStart-30
                    # save target region 8-30bp
                    #print(chrom+'\t'+str(regionStart)+'\t'+ str(regionStop))
                    targetDict += range(regionStart, regionStop)
               # if normal strand
               else:
                    targetJunctEnd = int(exonJuncts.split("|")[0].split('..')[1])
                    regionStart, regionStop = targetJunctEnd+8, targetJunctEnd+30
                    # save target region 8-30bp
                    #print(geNum, chrom, targetJunctStart, targetJunctEnd)
                    #print(chrom+'\t'+str(regionStart)+'\t'+ str(regionStop))

                    targetDict += range(regionStart, regionStop)
          #return
          # total = 0
          # for key in normPosList:
          #      total += len(normPosList[key])
          # print(total)
          # read in file and if is not in the snp pos list skip the line
          
          out = open(outFile, 'w')
          with open(snpFile) as inFile2:
               # skip first line
               next(inFile2)
               for line in inFile2:
                    chrom, pos, alleleFreq, freqFrac = line.split()[0].split('-')[0], int(line.split()[0].split('-')[1]), line.split()[-3], line.split()[-2]
                    # not the chrom we want
                    if chrom not in normPosList:
                         continue
                    else:
                         # not in the pos list
                         if pos not in normPosList[chrom]:
                              continue
                    out.write(chrom+'\t'+str(pos)+'\t'+alleleFreq +'\t'+freqFrac)
               

     def snpAnscState(self, blastResults, merge=False):
          with open(blastResults) as inFile:
               lines = [line.rstrip() for line in inFile]

          # if you havent cleaned the blast results file and merged multiple hits that are nearby
          if merge ==True:
               skipLines = 0
               for x in range(0, len(lines)):
                    # if there was a merge downstream we need to skip the line that was merged
                    if skipLines != 0:
                         skipLines -= 1
                         continue
                    splitLine = lines[x].split()
                    sampleStart, sampleEnd = int(splitLine[6]), int(splitLine[7])
                    mergePos = [sampleStart, sampleEnd]
                    for y in range(x+1, len(lines)):
                         splitLine2 = lines[y].split()
                         sampleStart2, sampleEnd2 = int(splitLine2[6]), int(splitLine2[7])
                         # check to see if looking at same sample region as parent loop
                         if splitLine[0] == splitLine2[0]:
                              skipLines += 1
                              # check to make sure the sample positions are within 100bp of each other 
                              if abs(sampleEnd-sampleStart2) < 100 or abs(sampleEnd2-sampleStart) < 100:
                                   mergePos[0], mergePos[1] = min(mergePos[0], sampleStart2), max(mergePos[1], sampleEnd2)
                                   #print(splitLine)
                                   #print(splitLine2)
                                   #print('\n')
                                   
                         else:
                              #print('\t'.join(splitLine[:6]) +'\t'+ str(mergePos[0]) +'\t'+ str(mergePos[1]) +'\t'+ '\t'.join(splitLine[8:]))
                              break



def bootstrap(inFile1, sampleSize, replicates, stat):
     # imports
     import random
     import numpy as np

     # # Read in files 
     with open(inFile1) as inFly:
        	windowList = [line.rstrip() for line in inFly]
	# Create lists to hold stats
     flyTajD = []
     # Allows user to choose which stat to look at. Tajimas D = D, Pi, Watersons Theta = @
     if stat == 'D':
          stat = 6
     if stat == 'Pi':
          stat = 2
     if stat == 'W':
          stat = 3
     # Format the stats values from file
     for line in windowList:
          splitLine = line.split()
          flyTajD.append(float(splitLine[stat]))
     flyTajD.sort()
     # Designate a value in the stats list that will act as a threshold for the bootstrap
     lowTail = flyTajD[int(len(flyTajD) * 0.05) - 1]
     # Uses sample size to calculate what amount of values is the 5% lower tail threshold 
     lowEle = int(sampleSize * 0.05)
     # Count the number of replicates that are greater than the 
     totalCount = 0.0
     for x in range(0, replicates):
          bs = random.sample(flyTajD, sampleSize) # Bootstrap using input file and indicated sample size
          bs.sort()
          # Count the number of times that a boostrap value is less than the threshold value
          count = 0.0
          for y in range(0, len(bs)):
               if float(bs[y])<float(lowTail): 
                    count = count + 1.0
          # if there are more than the expected number of boostrap values than update counter
          print(count, lowEle)
          if count > lowEle:
               totalCount = totalCount + 1.0
     print(totalCount/replicates)
#bootstrap('Pi.yak', 35660, 1000, 'D')

# # #Read in mut and hap data 
# svData = StructVariantData('/nobackup/rogers_research/Brandon/RemoteEdit/ComplexVarAlign/MayotteHMM/workflow/results/SvCalls/mapped_ul_hap.txt')
# svMut = svData.readSVFile()
# svData = StructVariantData('/nobackup/rogers_research/Brandon/RemoteEdit/ComplexVarAlign/MayotteHMM/workflow/results/SvCalls/allHaps.txt')
# svHap = svData.readHapFile()

# # # Create cluster object for multiple steps
# sv = StructVariants(svMut)

# ## clustering 
# # #sv.newCluster()
# cluster = sv.clusterSV()
# allClust = []
# for dab in cluster:
#      allClust.append(dab +'\t'+ cluster.get(dab))
#      print(dab +'\t'+ cluster.get(dab))


#synths read in cluster file instead of object
# finClustDict, allClust = {}, []
# with open('/users/bturne48/RogersLab/RemoteEdit/ParseIlluminaData/Scripts/post_cluster_myOg.txt') as inFile:
#      for line in inFile:
#           line = line.rstrip()
#           #print(line)
#           splitLine = line.split()
#           allClust.append(line)
#           finClustDict['\t'.join(splitLine[:6])] = '\t'.join(splitLine[6:])


# # find donor regions for variants
# outputLoc = '/users/bturne48/RogersLab/RemoteEdit/TempAsc/'
# x = sv.localVsGenomicCov(allClust, '/scratch/bturne48/NewBAM/')
# sv.heteroCorrection(x, '/users/bturne48/RogersLab/RemoteEdit/ParseIlluminaData/Data/OutputFromScripts/sort_genomicCovAllStrains.txt')
#highFreq = sv.findDonorRegion(x, outputLoc+'sort_genomicCovAllStrains.txt', outputLoc+'matchTE.txt')


# fasta info
# with open('/users/bturne48/RogersLab/RemoteEdit/ParseIlluminaData/StructVarResults/Reb_Old/rebCluster_post_old.txt') as inWeirdCov:
#      case = [line.rstrip() for line in inWeirdCov]   
#sv.variantsToFasta(allClust, '/scratch/bturne48/BlastDB/dyakRef.fasta', '/nobackup/rogers_research/Brandon/PipelineResults/MyNewCluster/toBlast.txt', 1000)


# polarize variants based on BLAST results
#polarList = sv.polarizeRearrangments('/nobackup/rogers_research/Brandon/PipelineResults/MyOgCluster/blastPolar.txt', 'polarizedVariantsList.txt')


# #Create sampled object 
# sampled = sv.checkSampledHaps(finDict, svHap)
# for dab in sampled:
#     print(dab)


# # # intron stuff
#sv.intergenicRegions('/users/bturne48/RogersLab/RemoteEdit/ParseIlluminaData/Data/MiscInput/ExonCoords.r1.05.txt',
#                   '/nobackup/rogers_research/rob/inbred/vcf-megatable/allelefrequency/pos-grey-cols-freqall.txt')


# Create sfs object with projection and polarization
#def calcSFS(self, clusteredVars, sampledHaps, targetPop, polarVars='', shouldProject=False, maxno=0, repbaseFilterTE=False, repbaseBlastResults=[], makeRepbaseSFS=False)
#sv.calcSFS(cluster, sys.argv[1], sys.argv[2], sys.argv[3], bool(sys.argv[4]), int(sys.argv[5]), bool(sys.argv[6]), sys.argv[7], bool(sys.argv[8]))
#allPop = sv.calcSFS(cluster, 'sampRes.txt', '', 'Polarize/polarOut.txt', True, 22, True, 'Blast/blastResultsRepbase.txt', True)
# oranPop = sv.calcSFS(cluster, '/users/bturne48/RogersLab/RemoteEdit/ParseIlluminaData/StructVarResults/Reseq_2sup/sampRes_post.txt',
#                     '(Oran|oran)', '',
#                      False, 22, False, '/scratch/bturne48/ParseIllumina/BigFileOutputs/blastResultsRepbase.txt', False)
# greyPop = sv.calcSFS(finClustDict, '/users/bturne48/RogersLab/RemoteEdit/ParseIlluminaData/StructVarResults/Reseq_2sup/sampRes_post.txt',
#                                '(CY|cy|NY|ny|Grey|grey)', '',
#                                   True, 19, False, '/scratch/bturne48/ParseIllumina/BigFileOutputs/blastResultsRepbase.txt', False)
# dsanPop = sv.calcSFS(finClustDict, '/users/bturne48/RogersLab/RemoteEdit/ParseIlluminaData/Scripts/sampRes_myOg_post.txt',
#                   '^(?:(?!CY|cy|NY|ny|Grey|grey|Oran|oran).)+$', '',
#                   False, 19, False, '/scratch/bturne48/ParseIllumina/BigFileOutputs/blastResultsRepbase.txt', False)

# synths
# oranPop = sv.calcSFS(finClustDict, '/users/bturne48/RogersLab/RemoteEdit/ParseIlluminaData/Data/OutputFromScripts/synth_sampRes_approx.txt',
#                     '(Oran|oran)', '/users/bturne48/RogersLab/RemoteEdit/ParseIlluminaData/Data/OutputFromScripts/polarizedVariantsList_synth.txt',
#                      False, 15, False, '/scratch/bturne48/ParseIllumina/BigFileOutputs/blastResultsRepbase_synth.txt', False)
# greyPop = sv.calcSFS(finClustDict, '/users/bturne48/RogersLab/RemoteEdit/ParseIlluminaData/Data/OutputFromScripts/synth_sampRes_approx.txt',
#                                '(CY|cy|NY|ny|Grey|grey)', '/users/bturne48/RogersLab/RemoteEdit/ParseIlluminaData/Data/OutputFromScripts/polarizedVariantsList_synth.txt',
#                                   False, 15, False, '/scratch/bturne48/ParseIllumina/BigFileOutputs/blastResultsRepbase_synth.txt', False)
# dsanPop = sv.calcSFS(finClustDict, '/users/bturne48/RogersLab/RemoteEdit/ParseIlluminaData/Data/OutputFromScripts/synth_sampRes_approx.txt',
#                   '^(?:(?!CY|cy|NY|ny|Grey|grey|Oran|oran).)+$', '/users/bturne48/RogersLab/RemoteEdit/ParseIlluminaData/Data/OutputFromScripts/polarizedVariantsList_synth.txt',
#                    False, 15, False, '/scratch/bturne48/ParseIllumina/BigFileOutputs/blastResultsRepbase_synth.txt', False)


#fst calcs
# svFST = sv.calcFST(finClustDict, '/users/bturne48/RogersLab/RemoteEdit/ParseIlluminaData/Data/OutputFromScripts/sampRes.txt',
#                      ['(Oran|oran)','(CY|cy|NY|ny|Grey|grey)','^(?:(?!CY|cy|NY|ny|Grey|grey|Oran|oran).)+$'],
#                      ["Oran", "CY/NY/Grey", "Dsan"],
#                      ["Dsan", "CY/NY/Grey"])


#sv.snpAnscState('/users/bturne48/RogersLab/RemoteEdit/TempShellOutput/mergeSnpBlastRes.txt', True)

#popgen mapping 
#sv.popgenMap(yakubaInput=False, greyInput=True, objInput=False)


# # stat sig vars
# statSigLists = sv.statSigVars("/users/bturne48/RogersLab/RemoteEdit/FinalOutputs/PopgenMapping/GreyVsDsan/popgen_yak_greyVSdsan.txt",
#                                statIndexTarget=7, mathstuff=True, detail=False)


# sv.assMethod('/users/bturne48/RogersLab/RemoteEdit/FinalOutputs/AlleleFreqs/ass_dsanAF.txt', 
#             '/users/bturne48/RogersLab/RemoteEdit/FinalOutputs/AlleleFreqs/ass_oranAF.txt')
