##ordered sample cluster

from optparse import OptionParser
import subprocess
import time
import os
import numpy

parser=OptionParser()
parser.add_option("-k","--kcluster",dest="k",type="int")#highest cluster number 
#parser.add_option("-o","--output",dest="output",type="str")#output_cluster
parser.add_option("-i","--input",dest="input",type="str")#input genelist
parser.add_option("-d","--directory",dest="directory",type="str")#directory
parser.add_option("-c","--chrom",dest="chrom",type="str")#directory
parser.add_option("-r","--resolution",dest="resolution",type="int")#resolution
(options,args)=parser.parse_args()



def seedinitiate(seed,list): ##D
   seed.append(list[0][0])  


def seedfresh(seed,list):
   length=len(seed)
   seed.append(list[length][length]) #add new D
   for loop in range(length):
      seed[loop]+=sum(list[length][loop:])

      
def seedcalculate(seed,error,index,k):#k=higher number of clusters
   length=len(seed)
   error[length][1]= seed[0]###maxvalue maxindex
   index[length][1]=1

   for loop in range(1,min(k,length)):
      inside=[]
      for upper in range(loop+1,length+1):
         inside.append(seed[upper-1]+error[upper-1][loop])
      maxvalue=max(inside)
      maxindex=inside.index(maxvalue)

      error[length][loop+1]=maxvalue
      index[length][loop+1]=maxindex+loop+1


####read_human
def read_human(filename):
   input=open(filename)
   list=[]
   number=0
   while True:
      line=input.readline()
      if not line:break
      number+=1
      line=line.split()
      list.append([float(line[i]) for i in range(3,3+number)])
      
   return list,number

####read_fly
def read_fly(filename):
   input=open(filename)
   input.readline()
   list=[]
   number=0
   while True:
      line=input.readline()
      if not line:break
      number+=1
      line=line.split()
      list.append([float(line[i]) for i in range(1,1+number)])
      
   return list,number

def error_initation(n,k):
   a=numpy.array([0.0]*(n+2)*(k+2)) ###pay attention to initatied value assignment
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
   seedinitiate(seed)
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
   
   #ncluster=input('Picture has been generated, Please input the number of clusters'+'(<'+str(options.k)+'): ')
   ncluster=options.k
   nordered=total_number
   
   while True:
      if ncluster<=0:break
      compart=[]
      current_cluster=ncluster
      current_ordered=nordered
      while True:
         if current_cluster<2:break
         #file=open(options.directory+'/clusterresult/'+str(ncluster)+'.txt','w')
         start=index[current_ordered][current_cluster]
         compart.insert(0,int(start))
         current_cluster-=1
         current_ordered=int(start)-1

      compart.insert(0,1)
      compart.append(total_number+1)

      #print compart
      output=open(options.directory+'/'+options.chrom+'_'+str(ncluster)+'.bed','w')
      for i in range(len(compart)-1):
         content=[options.chrom,str((compart[i]-1)*int(options.resolution)),str((compart[i+1]-1)*int(options.resolution)),str(i%2)]
         output.writelines("\t".join(content)+'\n')
      output.close()
      
      ncluster-=10


if __name__=='__main__':
    main() 
