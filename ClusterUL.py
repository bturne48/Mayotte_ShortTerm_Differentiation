import sys, gzip
InFile=gzip.open(sys.argv[1], 'r')
thresh=325
Dict={}
prev=[0,0,0]
count=1




for line in InFile:
  dummy=0
  A=line.split()
  rearrtype=A[0]
  if A[1][0]!='G' and A[1][0]!='M'and A[4][0]!='G' and A[4][0]!='M':
    if A[1]!='Y':
        chrom1=A[1]
        start1=int(A[2])
        chrom2=A[4]
        start2=int(A[5])
    else:
        chrom2=A[1]
        start2=int(A[2])
        chrom1=A[4]
        start1=int(A[5])
    if abs(start1-start2)>0:
      if Dict.has_key(chrom1):
        dummy=0
        for i in range(0, len(Dict[chrom1])):
            X=Dict[chrom1]
            prev=X[i]
            if prev[2]==chrom2 and abs(prev[0]-start1)<thresh and abs(prev[3]-start2)<thresh and rearrtype==prev[6]:
                newstart1=min(prev[0], start1)
                newstop1=max(prev[1], start1)
                newstart2=min(start2, prev[3])
                newstop2=max(start2, prev[4])
                count=prev[5]+1
                X[i]=[newstart1, newstop1, chrom2, newstart2, newstop2, count, rearrtype]
                Dict[chrom1]=X
                dummy=1
                break
        if dummy==0:
                Dict[chrom1].append([start1, start1, chrom2, start2,start2 ,1, rearrtype])
      else:
        Dict[chrom1]=[[start1, start1, chrom2, start2, start2, 1, rearrtype]]
#for item in Dict['11']:
 # if item[2]=='11':
  #  print item
#print Dict['10']


for kk in range(0,5):
  for key in Dict.keys():
    for i in range(0,len(Dict[key])):
        X=Dict[key]
        chrom1=key
        for j in range(i+1, len(Dict[key])):
            prev=X[i]
            check=X[j]
            chrom2=check[2]
            rearrtype=check[6]
            if abs(prev[0]-check[0])<thresh and abs(prev[3]-check[3])<thresh and prev[2]==chrom2 and prev[6]==rearrtype:
                newstart1=min(prev[0], check[0])
                newstop1=max(prev[1], check[1])
                newstart2=min(check[3], prev[3])
                newstop2=max(check[4], prev[4])
                count=prev[5]+check[5]
                X[i]=[newstart1, newstop1, chrom2, newstart2,newstop2, count, rearrtype]
                del X[j]
                Dict[chrom1]=X
                break
            elif abs(prev[1]-check[1])<thresh and abs(prev[4]-check[4])<thresh and prev[2]==chrom2 and prev[6]==rearrtype:
                newstart1=min(prev[0], check[0])
                newstop1=max(prev[1], check[1])
                newstart2=min(check[3], prev[3])
                newstop2=max(check[4], prev[4])
                count=prev[5]+check[5]
                X[i]=[newstart1, newstop1, chrom2, newstart2,newstop2, count, rearrtype]
                del X[j]
                Dict[chrom1]=X
                break
            elif abs(prev[1]-check[1])<thresh and abs(prev[3]-check[4])<thresh and prev[2]==chrom2 and prev[6]==rearrtype:
                newstart1=min(prev[0], check[0])
                newstop1=max(prev[1], check[1])
                newstart2=min(check[3], prev[3])
                newstop2=max(check[4], prev[4])
                count=prev[5]+check[5]
                X[i]=[newstart1, newstop1, chrom2, newstart2,newstop2, count, rearrtype]
                del X[j]
                Dict[chrom1]=X
                break
            elif abs(prev[1]-check[1])<thresh and abs(prev[4]-check[3])<thresh and prev[2]==chrom2 and prev[6]==rearrtype:
                newstart1=min(prev[0], check[0])
                newstop1=max(prev[1], check[1])
                newstart2=min(check[3], prev[3])
                newstop2=max(check[4], prev[4])
                count=prev[5]+check[5]
                X[i]=[newstart1, newstop1, chrom2, newstart2,newstop2, count, rearrtype]
                del X[j]
                Dict[chrom1]=X
                break
                

OutFile=gzip.open(sys.argv[2], 'w')
for key in Dict.keys():
    for item in Dict[key]:
        if item[5]>2 and item[5]<100:
            print >>OutFile, item[6], key, item[0], item[1], item[2], item[3], item[4], item[5]
