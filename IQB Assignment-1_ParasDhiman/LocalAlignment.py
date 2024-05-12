txt1 = "GATGCGCAG"
txt2 = "GGCAGTA"
rowNum = len(txt2) + 1
colNum = len(txt1) + 1
print(rowNum)
print(colNum)
totalRowCol = rowNum * colNum
print("total number of rows and columns would be",totalRowCol)

lst = []

for i in range(0,rowNum):
    lst1 = []
    for j in range(0,colNum):
        lst1.append(0)
    lst.append(lst1)



match = 2
missMatch = -1
gap = -3

for i in range(0,rowNum):
    for j in range(0,colNum):
        
        if (i == 0 or j == 0):
            # here we have checked if i and j both are zero then nothing to be done as already initialised with zero
            if (j == 0 and i == 0):
                lst[i][j] = 0
            #Assigning the columns with there values subtracting gap from them
            elif(i == 0):
                lst[i][j] = 0
            #Assigning the rows with there values subtracting gap from them
            else:
                lst[i][j] = 0
        else:
            #checking if both the nucleotide pair is same then its a match and the score of match is added among them
            if(txt1[j - 1] == txt2[i - 1]):
                lst[i][j] = lst[i - 1][j - 1] + match
            #checking which is more bigger the diagonal one after subtracting the penalty of missmatch or upper or left one after subtracting the gap so maximum value is assigned to the list
            else:
                diagonal = lst[i-1][j-1] + missMatch
                up = lst[i-1][j] + gap
                left = lst[i][j-1] + gap
                lst[i][j] = max(diagonal, up, left)
        if(lst[i][j] < 0):
            lst[i][j] = 0
for i in lst:
    print(i)

str1 = ""
str2 = ""
row = len(txt2)
col = len(txt1)

highest = -1
idx1 = []
idx2 = []

for i in range (0,row + 1):
    for j in range(0,col + 1):
        if (lst[i][j] > highest):
            highest = lst[i][j]

for i in range (0,row + 1):
    for j in range(0,col + 1):
        if (lst[i][j] == highest):
            idx1.append(i)
            idx2.append(j)

alignment = []
for p in range(0,len(idx1)):
    trackLst = [(idx1[p], idx2[p], str1, str2)]
    k1 = -1
    k2 = -1
    while trackLst:
        #Assigning them the values 
        i, j, str1, str2 = trackLst.pop()
        if ((lst[i][j] == 0 or (i == 0 and j == 0))):
            # trackLst.append((0,0, str1 + txt2[i - 1], str2 + txt1[j - 1]))
            alignment.append((str1[::-1], str2[::-1]))

        else:
            #here all the four possibiliies have been tracked if fount then create a new path else continue evaluating the certain trees
            if (i > 0 and lst[i][j] == lst[i - 1][j] + gap ):
                trackLst.append((i - 1, j, str1 + txt2[i - 1], str2 + "-"))
            if (j > 0 and lst[i][j] == lst[i][j - 1] + gap ):
                trackLst.append((i, j - 1, str1 + "-", str2 + txt1[j - 1]))
            if (i > 0 and j > 0 and lst[i][j] == lst[i - 1][j - 1] + match and txt1[j - 1] == txt2 [i - 1]):
                trackLst.append((i - 1, j - 1, str1 + txt2[i - 1], str2 + txt1[j - 1]))
            if (i > 0 and j > 0 and lst[i][j] == lst[i - 1][j - 1] + missMatch ):
                trackLst.append((i - 1, j - 1, str1 + txt2[i - 1], str2 + txt1[j - 1]))
    trackLst = []
    str1 = ""
    str2 = ""

count = 0

#If there is more than one possibility
if(len(alignment) > 1):
    print("Yes, there will be more than one optical alignment")
else:
    print("No,There is only one Optical alignment possible")

#Printing of the different optical Alignments
for i in alignment:
    print(i[0])
    print(i[1])
    #calculating the score
    score = 0
    for j in range(0,len(i[0])):
        if(i[0][j] == i[1][j]):
            score += match
        elif(i[0][j] != i[1][j] and (i[0][j] == "-" or i[1][j] == "-" )):
            score += gap
        else:
            score += missMatch
    print("The score is: ", score)
    print("------------------")
    count += 1

print("Total Alignments are", count)
