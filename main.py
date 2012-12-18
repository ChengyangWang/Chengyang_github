###peak_association

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-b","--boundary",dest="boundary",type="str")#boundary_file
parser.add_option("-p","--peak",dest="peak",type="str")###peaks
parser.add_option("-e","--expression",dest="expression",type="str")###expression
parser.add_option("-o","--output",dest="output",type="str")##output.csv
(options,args)=parser.parse_args()


###output {chrom:[[start,end,boundary_score]]}
def read_boundary_score(boundary_score_file):
   file=open(boundary_score_file) 
   file.readline()
   file.readline()
   boundary={"chr19":[]}
   while True:
      line=file.readline()
      if not line:break
      line=line.strip().split()
      boundary["chr19"].append([int(line[0]),int(line[0])+5000,float(line[1])])
   boundary["chr19"].sort(key=lambda x:x[2],reverse=True)
   for i in range(len(boundary["chr19"])):
      boundary["chr19"][i][2]=1-float(i)/len(boundary["chr19"])
   boundary["chr19"].sort(key=lambda x:x[0])

   print "finish boundary file reading\n"
   return boundary
      
##output {chrom:[[start,end,signal]]}
def read_peak(peak_file):
   file=open(peak_file)
   a={}
   while True:
      line=file.readline()
      if not line:break
      line=line.split()
      try:
         a[line[0]].append([int(line[1]),int(line[2]),float(line[4])])
      except KeyError:
         a[line[0]]=[[int(line[1]),int(line[2]),float(line[4])]]

   print "finish peak file reading\n"
   return a

### output {chrom:[[symbol,TSS,TTS,FC]]}
def read_expression(expression_file):
   file=open(expression_file)
   file.readline()
   output_list={}
   while True:
      line=file.readline()
      if not line:break
      line=line.split()
      if line[2]=="+":
         try:
            output_list[line[1]].append([line[0],int(line[3]),int(line[4]),float(line[7])])
         except KeyError:
            output_list[line[1]]=[[line[0],int(line[3]),int(line[4]),float(line[7])]]
      else:
         try:
            output_list[line[1]].append([line[0],int(line[4]),int(line[3]),float(line[7])])
         except KeyError:
            output_list[line[1]]=[[line[0],int(line[4]),int(line[3]),float(line[7])]]

   print "finish expression file reading\n"
   return output_list


##rapid searching
def searchend(sequence,number,lower,upper):
    if lower==upper-1:
        if sequence[lower][0]<=number<sequence[upper][0]: return lower
        elif number<sequence[lower][0]: return lower-1   #if data<all
        else: return upper
    else:
        middle=(lower+upper)/2
        if number>sequence[middle][0]:
            return searchend(sequence,number,middle,upper)
        else:
            return searchend(sequence,number,lower,middle)

###core_function scope=absolute distance
def core_function(peak_position_index,peak_signal,chrom_boundary_list,scope=100000):
   basic_domain_width=chrom_boundary_list[0][1]-chrom_boundary_list[0][0]
   extent_index=scope/basic_domain_width
   
   right=[[chrom_boundary_list[peak_position_index][0],chrom_boundary_list[peak_position_index][1],peak_signal]]
   for i in range(1,min(extent_index,len(chrom_boundary_list)-peak_position_index-1)+1):
      impact=-right[i-1][2]*(0.1*chrom_boundary_list[peak_position_index+i-1][2]-1)
      right.append([chrom_boundary_list[peak_position_index+i][0],chrom_boundary_list[peak_position_index+i][1],impact])

   left=[[chrom_boundary_list[peak_position_index][0],chrom_boundary_list[peak_position_index][1],peak_signal]]
   for i in range(1,min(extent_index,peak_position_index)+1):
      impact=-left[i-1][2]*(0.1*chrom_boundary_list[peak_position_index-i+1][2]-1)
      left.append([chrom_boundary_list[peak_position_index-i][0],chrom_boundary_list[peak_position_index-i][1],impact])

   left.reverse()
   
   return left[:-1]+right  ###[[start,end,influence]]
      
def main():
   boundary_list=read_boundary_score(options.boundary)
   peak_list=read_peak(options.peak)
   expression_list=read_expression(options.expression)
   output_file=open(options.output,'w')
   influence_dic={} ##{chrom:[]}
   
   chrom_list=['chr'+v for v in map(str,range(1,30))+['X','Y']]


   ###influence umbrella
   for chrom in chrom_list:
      try:
         chrom_peak_list=peak_list[chrom]
         chrom_boundary_list=boundary_list[chrom]
         length_chrom_boundary_list=len(chrom_boundary_list)
         influence_dic[chrom]=[]
      except KeyError:
         continue
      for peak in chrom_peak_list:
         
         peak_middle=(int(peak[0])+int(peak[1]))/2
         peak_signal=float(peak[2])
         peak_position_index=searchend(chrom_boundary_list,peak_middle,0,length_chrom_boundary_list-1)
         influence=core_function(peak_position_index,peak_signal,chrom_boundary_list)
         influence_dic[chrom]+=influence
         "handle one peak"


   ###integrate influence umbrella
   for chrom in influence_dic:
      influence_dic[chrom].sort(key=lambda x:int(x[0]))
      i=-1
      inf=0
      result=[]
      for v in influence_dic[chrom]:
         if int(v[0])!=i and i!=-1:
            result[-1].append(inf)
            inf=float(v[2])
            result.append([v[0],v[1]])
            i=int(v[0])
         elif int(v[0])!=i and i==-1:
            inf=float(v[2])
            result.append([v[0],v[1]])
            i=int(v[0])
         else:
            inf+=float(v[2])
      
      influence_dic[chrom]=result
      
   ###influences are mapped to gene {chrom:[[symbol,TSS,FC]]}
   for chrom in expression_list:
      try:
         chrom_influence_list=influence_dic[chrom]
      except KeyError:
         continue

      
      for sample in expression_list[chrom]:
         TSS_position=searchend(chrom_influence_list,int(sample[1]),0,len(chrom_influence_list)-1)
         if TSS_position==-1 or TSS_position>len(chrom_influence_list)-1 or sample[3]<=1:
            continue
         output_file.writelines(sample[0]+' '+str(chrom_influence_list[TSS_position][2])+' '+str(sample[1])+' '+
                                str(sample[2])+' '+str(sample[3])+'\n')
   output_file.close()

main()

         
            
         
         
         





         
   
