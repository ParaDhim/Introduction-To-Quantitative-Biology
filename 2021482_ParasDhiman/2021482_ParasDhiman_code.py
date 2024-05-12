# Chou and Fasman method of secondary structure prediction
# Chou and Fasman parameters

lst1 = ["E","A","L","H","M","Q","W","V","F","K","I","D","T","S","R","C","N","Y","P","G"]
lst2 = ["Glu","Ala","Leu","His","Met","Gln","Trp","Val","Phe","Lys","Ile","Asp","Thr","Ser","Arg","Cys","Asn","Tyr","Pro","Gly"]
pAlpha = [1.53,1.45,1.34,1.24,1.20,1.17,1.14,1.14,1.12,1.07,1.00,0.98,0.82,0.79,0.79,0.77,0.73,0.61,0.59,0.53]
pBeta = [0.26,0.97,1.22,0.71,1.67,1.23,1.19,1.65,1.28,0.74,1.60,0.80,1.20,0.72,0.90,1.30,0.65,1.29,0.62,0.81]
protSeq = "SGFRKMAFPSGKVEGCMVQVTCGTTTLNGLWLDDTVYCPRHVICTAEDMLNPNYEDLLIRKSNHSFLVQAGNVQLRVIGHSMQNCLLRLKVDTSNPKTPKYKFVRIQPGQTFSVLACYNGSPSGVYQCAMRPNHTIKGSFLNGSCGSVGF"

#window size for helix and Stand has set over here
winSizeHelix = 6
winSizeStrand = 5


#predicting the nucleation site for the helix
coil = 0
Hel = 1
Strand = 2

#generating the site of nucleation for helices
siteHel = []
siteHelIdxStart = []
siteHelIdxEnd = []

#
for i in range(0,len(protSeq)):
    if (i+winSizeHelix <= len(protSeq)):
        siteHel.append(protSeq[i:i+winSizeHelix])
        siteHelIdxStart.append(i)
        siteHelIdxEnd.append(i+winSizeHelix)
    


#generating the site of nucleation for strands
siteStr = []
siteStrIdxStart = []
siteStrIdxEnd = []

for i in range(0,len(protSeq)):
    if (i+winSizeStrand<= len(protSeq)):
        siteStr.append(protSeq[i:i+winSizeStrand])
        siteStrIdxStart.append(i)
        siteStrIdxEnd.append(i+winSizeStrand)


# Helices and Strands sequence have been checked over here that they must have there score greater than 1 has that too must also has to have 4 in those series
siteHelChecked = []
siteStrChecked = []
siteHelCheckedIdxs= []
siteStrCheckedIdxs = []


# here for helics have been checked
for i in siteHel:
    count = 0
    for j in i:
        if pAlpha[(lst1.index(j))] >= 1:
            count += 1
    if count >= 4:
        siteHelChecked.append(i)
        ll = []
        for r in range(siteHelIdxStart[siteHel.index(i)],siteHelIdxStart[siteHel.index(i)] + winSizeHelix):
            ll.append(r)
        siteHelCheckedIdxs.append(ll)
        
# here for strands have been checked
for i in siteStr:
    count = 0
    for j in i:
        if pBeta[(lst1.index(j))] >= 1:
            count += 1
    if count >= 3:
        siteStrChecked.append(i)
        ll = []
        for r in range(siteStrIdxStart[siteStr.index(i)],siteStrIdxStart[siteStr.index(i)] + winSizeStrand):
            ll.append(r)
        siteStrCheckedIdxs.append(ll)
        
# a dummy list has been created that stores the the indices has to have helices or strands and then these will be compared
HAdder = []
for i in siteHelCheckedIdxs:
    for j in i:
        if j not in HAdder:
            HAdder.append(j)

SAdder = []
for i in siteStrCheckedIdxs:
    for j in i:
        if j not in SAdder:
            SAdder.append(j)

# here helix and strand has been made which will store the H and S acoording to the inde position
hel = []
stra = []

for i in range (0,len(protSeq)):
    hel.append("-")
    stra.append("-")

#extending the sequence for helix
for i in siteHelChecked:
    if(siteHelCheckedIdxs[siteHelChecked.index(i)][0] == 0):
        # extending this to the right
        left = siteHelCheckedIdxs[siteHelChecked.index(i)][0]+3
        right = left + 4
        temp = protSeq[left:right]
        calc = 4
        # here we have calculated  that the score is greater than or if not then it will stop
        while(calc >= 4):
            calc = 0
            for j in temp:
                calc += pAlpha[lst1.index(j)]
            if (calc >= 4):# we can say that the average P(H) < 1  if this condition come out as false
                if (right - 1) not in HAdder:
                    HAdder.append(right-1)
                if(right - 1 < len(protSeq)-1):
                    left += 1
                    right += 1
                    temp = protSeq[left:right]
                else:
                    break
    elif (siteHelCheckedIdxs[siteHelChecked.index(i)][len(i)-1] == len(protSeq)-1):
        # extention to the left
        right = siteHelCheckedIdxs[siteHelChecked.index(i)][len(i)-1]-2

        left = right - 4
        temp = protSeq[left:right]
        calc = 4
        # here we have calculated  that the score is greater than or if not then it will stop
        while(calc >= 4):
            calc = 0
            for j in temp:
                calc += pAlpha[lst1.index(j)]
            if (calc >= 4):# we can say that the average P(H) < 1  if this condition come out as false
                if left not in HAdder:
                    HAdder.append(left)
                if(left - 1 >= 0):
                    left -= 1
                    right -= 1
                    temp = protSeq[left:right]
                else:
                    break
    else:
        # extention to the left for general cases
        right = siteHelCheckedIdxs[siteHelChecked.index(i)][len(i)-1]-2
        left = right - 4
        temp = protSeq[left:right]
        calc = 4
        # here we have calculated  that the score is greater than or if not then it will stop
        while(calc >= 4):
            calc = 0
            for j in temp:
                calc += pAlpha[lst1.index(j)]
            if (calc >= 4):# we can say that the average P(H) < 1  if this condition come out as false
                if left not in HAdder:
                    HAdder.append(left)
                if(left - 1 >= 0):
                    left -= 1
                    right -= 1
                    temp = protSeq[left:right]
                else:
                    break
        #extending this only to the right
        left = siteHelCheckedIdxs[siteHelChecked.index(i)][0]+3
        right = left + 4
        temp = protSeq[left:right]
        calc = 4
        # here we have calculated  that the score is greater than or if not then it will stop
        while(calc >= 4):
            calc = 0
            for j in temp:
                calc += pAlpha[lst1.index(j)]
            if (calc >= 4):# we can say that the average P(H) < 1  if this condition come out as false
                if (right - 1) not in HAdder:
                    HAdder.append(right-1)
                if(right - 1 < len(protSeq)-1):
                    left += 1
                    right += 1
                    temp = protSeq[left:right]
                else:
                    break
                    
                    
for k in HAdder:
    hel[k] = "H"


for i in siteStrChecked:
    if(siteStrCheckedIdxs[siteStrChecked.index(i)][0] == 0):
        #extending this only to the right
        left = siteStrCheckedIdxs[siteStrChecked.index(i)][0]+2
        right = left + 4
        temp = protSeq[left:right]
        calc = 4
        # here we have calculated  that the score is greater than or if not then it will stop
        while(calc >= 4):
            calc = 0
            for j in temp:
                calc += pBeta[lst1.index(j)]
            if (calc >= 4):# we can say that the average P(S) < 1 if this condition come out as false
                if (right - 1) not in SAdder:
                    SAdder.append(right-1)
                if(right - 1 < len(protSeq)-1):
                    left += 1
                    right += 1
                    temp = protSeq[left:right]
                else:
                    break
    
    elif (siteStrCheckedIdxs[siteStrChecked.index(i)][len(i)-1] == len(protSeq)-1):
        #extention to the left
        right = siteStrCheckedIdxs[siteStrChecked.index(i)][len(i)-1]-1
        left = right - 4
        temp = protSeq[left:right]
        calc = 4
        # here we have calculated  that the score is greater than or if not then it will stop
        while(calc >= 4):
            calc = 0
            for j in temp:
                calc += pBeta[lst1.index(j)]
            if (calc >= 4):# we can say that the average P(S) < 1 if this condition come out as false
                if left not in SAdder:
                    SAdder.append(left)
                if(left - 1 >= 0):
                    left -= 1
                    right -= 1
                    temp = protSeq[left:right]
                else:
                    break
    else:
        #extending to the right
        left = siteStrCheckedIdxs[siteStrChecked.index(i)][0]+2
        right = left + 4
        temp = protSeq[left:right]
        calc = 4
        # here we have calculated  that the score is greater than or if not then it will stop
        while(calc >= 4):
            calc = 0
            for j in temp:
                calc += pBeta[lst1.index(j)]
            if (calc >= 4):# we can say that the average P(S) < 1 if this condition come out as false
                if (right - 1) not in SAdder:
                    SAdder.append(right-1)
                if(right - 1 < len(protSeq)-1):
                    left += 1
                    right += 1
                    temp = protSeq[left:right]
                else:
                    break
                
        #extention to the left
        right = siteStrCheckedIdxs[siteStrChecked.index(i)][len(i)-1]-1
        left = right - 4
        temp = protSeq[left:right]
        calc = 4
        # here we have calculated  that the score is greater than or if not then it will stop
        while(calc >= 4):
            calc = 0
            for j in temp:
                calc += pBeta[lst1.index(j)]
            if (calc >= 4):# we can say that the average P(S) < 1 if this condition come out as false
                if left not in SAdder:
                    SAdder.append(left)
                if(left - 1 >= 0):
                    left -= 1
                    right -= 1
                    temp = protSeq[left:right]
                else:
                    break
               
               
# here the value of S has been submitted in the list for each index where there must be the strand 
for k in SAdder:
    stra[k] = "S"

simiIdx = []
finalSeq = []
# new final list has been made and now appended with - for denoting the full list as -
for i in range(0,len(protSeq)):
    finalSeq.append("-")

#here helix sequence has been put in the final sequence
for i in range(0,len(protSeq)):
    if (hel[i] != "-"):
        finalSeq[i] = hel[i]
        
#here strand sequence has been put in the final sequence
for i in range(0,len(protSeq)):
    if (stra[i] != "-"):
        finalSeq[i] = stra[i]


# Now finding out the intersection for the helix and strand
i = 0
while(i != len(protSeq)):
    if(hel[i] == "H" and stra[i] == "S"):
        f = i
        ll = []
        while(hel[f] == "H" and stra[f] == "S" and f < len(protSeq)):
            ll.append(f)
            if(f+1 == len(protSeq)):
                break
            f += 1
        simiIdx.append(ll)
        i = f
    else:
        if(hel[i] != "-"):
            finalSeq[i] = hel[i]
        else:
            finalSeq[i] = stra[i]
    
    i += 1


# similar index has been noted in the simiIdx list than check the score for each one and then handled those where they both share the common index
for i in simiIdx:
    sumHel = 0
    sumStr = 0
    for idx in i:
        sumHel += pAlpha[(lst1.index(protSeq[idx]))]
        sumStr += pBeta[(lst1.index(protSeq[idx]))]
    if(sumHel > sumStr):
        for idx1 in i:
            finalSeq[idx1] = "H"
    else:
        for idx2 in i:
            finalSeq[idx2] = "S"
            
            
finalSeqString = ""

#final string has been soted now in the form of string 
for i in finalSeq:
    finalSeqString += i


# here the stored string has been putted up in set of 50-50 just to display them properly
print("These are the sequences printed 50 in each line:")
print("1","     ",protSeq[0:50])
print(" ","     ",finalSeqString[0:50])
print("2","     ",protSeq[50:100])
print(" ","     ",finalSeqString[50:100])
print("3","     ",protSeq[100:150])
print(" ","     ",finalSeqString[100:150])