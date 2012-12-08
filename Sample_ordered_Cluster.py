##ordered sample cluster

from optparse import OptionParser
import subprocess
import time
import os
import numpy

parser=OptionParser()
parser.add_option("-k","--kcluster",dest="k",type="int")#higher cluster number 
parser.add_option("-o","--output",dest="output",type="str")#output_cluster
parser.add_option("-i","--input",dest="input",type="str")#input genelist
parser.add_option("-d","--directory",dest="directory",type="str")#directory
parser.add_option("-c","--chrom",dest="chrom",type="str")#directory
parser.add_option("-r","--resolution",dest="resolution",type="int")#resolution
(options,args)=parser.parse_args()



def seedinitiate(seed,list): ##variance mean
   length=len(list[0])
   var=[0 for i in range(length)]
   mean=[float(list[0][i]) for i in range(length)]
   seed.append([var,mean])  ####[[var1,var2....vark] [mean1,mean2....meank]]


def seedfresh(seed,list):
   length=len(seed)
   for loop in range(length):
      for component in range(len(seed[0][0])):
         seed[loop][0][component]+=((length-loop)*(list[length][component]-seed[loop][1][component])*(list[length][component]-seed[loop][1][component]))/float((length-loop)+1)
         seed[loop][1][component]=((length-loop)*seed[loop][1][component]+list[length][component])/float((length-loop)+1)
   
   var=[0 for i in range(len(list[0]))]
   mean=list[length]

   seed.append([var,mean])
   #print seed
      
def seedcalculate(seed,error,index,k):#k=higher number of clusters
   length=len(seed)
   error[length][1]=sum(seed[0][0])###minvalue minindex
   index[length][1]=1
   for loop in range(1,min(k,length)):
      inside=[]
      for upper in range(loop+1,length+1):
         inside.append(sum(seed[upper-1][0])+error[upper-1][loop])
      minvalue=min(inside)
      minindex=inside.index(minvalue)

      error[length][loop+1]=minvalue
      index[length][loop+1]=minindex+loop+1



def read(filename):
   input=open(filename)
   list=[]
   number=0
   while True:
      line=input.readline()
      if not line:break
      line=line.split()
      list.append([float(line[i]) for i in range(len(line))])
      number+=1
   return list,number


def error_initation(n,k):
   a=numpy.array([0.0]*(n+2)*(k+2))
   a.shape=n+2,k+2
   return a

def index_initation(n,k):
   a=numpy.array([0]*(n+2)*(k+2))
   a.shape=n+2,k+2
   return a




def main():
   list=read(options.input)[0]
   total_number=read(options.input)[1]
   error=error_initation(total_number,options.k)
   index=index_initation(total_number,options.k)
   seed=[]
   seedinitiate(seed,list)
   for fresh in range(total_number):
      
      starttime=time.time()
      seedcalculate(seed,error,index,options.k)
      if fresh!=total_number-1:
         seedfresh(seed,list)
      else:
         pass
      endtime=time.time()
      print 'Fresh '+str(fresh)+'\'time: '+str(endtime-starttime)
   del seed[:]

   #print error
   #print index
   #for j in range(2,(options.k)+1):
      #for i in range(j,(options.n)+1):
        # print dic[(i,j)]
   #print error
   #print index
   Pfile=open(options.directory+'/picture.R','w')

   index2=[]
   value=[]
   for i in range(1,options.k+1):
      index2.append(str(i))
      value.append(str(error[total_number][i]))               
   Rindex=','.join(index2)
   Rvalue=','.join(value)
   Pfile.writelines("""pdf("%s")\n"""%(options.directory+'/picture.pdf'))
   Pfile.writelines('plot(c('+Rindex+'),c('+Rvalue+'),'+"type='l',main='Variance picture',ylab='The variance in cluster',xlab='The number of clusters')")
   #Pfile.writelines('plot(c('+','.join(index)+')~c('+','.join(value)+'))')

   Pfile.writelines("""\ndev.off()""")
   Pfile.close()
   
   ncluster=input('Picture has been generated, Please input the number of clusters'+'(<'+str(options.k)+'): ')
   nordered=total_number
   compart=[]
   while True:
      if ncluster<2:break
      #file=open(options.directory+'/clusterresult/'+str(ncluster)+'.txt','w')
      start=index[nordered][ncluster]
      compart.insert(0,int(start))
      ncluster-=1
      nordered=int(start)-1

   compart.insert(0,1)
   compart.append(total_number+1)

   print compart
   output=open(options.output,'w')
   for i in range(len(compart)-1):
      content=[options.chrom,str((compart[i]-1)*int(options.resolution)),str((compart[i+1]-1)*int(options.resolution)),str(i%2)]
      output.writelines("\t".join(content)+'\n')

   output.close()


if __name__=='__main__':
    main() 
