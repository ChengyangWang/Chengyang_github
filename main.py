###peak_association

from optparse import OptionParser
from bx.bbi.bigwig_file import BigWigFile
from scipy import spatial

parser = OptionParser()
parser.add_option("-r","--refseq",dest="refseq",type="str")#refseq_file
parser.add_option("--header",dest="header",type="int")#header
parser.add_option("-b","--boundary",dest="boundary",type="str")#boundary_file
parser.add_option("-p","--peak",dest="peak",type="str")###peaks
parser.add_option("-e","--expression",dest="expression",type="str")###expression
parser.add_option("-w","--bw",dest="wig",type="string",action="append")
#parser.add_option("-o","--output",dest="output",type="str")##output.csv
(options,args)=parser.parse_args()


###output {chrom:[0,1000,5000...]}
def read_boundary_score(boundary_score_file):
   file=open(boundary_score_file) 
   boundary={}
   while True:
      line=file.readline()
      if not line:break
      line=line.strip().split()
      try:
         boundary[line[0]].append(int(line[2]))
      except KeyError:
         boundary[line[0]]=[0,int(line[2])]

   print "finish boundary file reading\n"
   return boundary
      
##output {chrom:[[peak_ID, summit, signal]]}
def read_peak(peak_file):
   file=open(peak_file)
   peak_chrom_dic={}
   intensity=[]
   while True:
      line=file.readline()
      if not line:break
      line=line.split()
      try:
         peak_chrom_dic[line[0]].append([line[3],int(line[1]),float(line[4])])
         intensity.append(float(line[4]))
      except KeyError:
         peak_chrom_dic[line[0]]=[[line[3],int(line[1]),float(line[4])]]
         intensity.append(float(line[4]))

   intensity.sort()
   intensity_dic={}

   for i in range(len(intensity)):
      intensity_dic[intensity[i]]=float(i+1)/len(intensity)   
   

   for key in peak_chrom_dic.keys():
      peak_chrom=peak_chrom_dic[key]
      for value in peak_chrom:
         value[2]=intensity_dic[value[2]]

   print "finish peak file reading\n"
   return peak_chrom_dic

####output {refseq_ID:[chrom,TSS,gene_ID]}
def read_anotation(refseq_file):
   file=open(refseq_file)
   refseq_dic={}
   while True:
      line=file.readline()
      if not line:break
      line=line.strip().split("|")
      if line[2]=="+":
         refseq_dic[line[1]]=[line[0],int(line[3]),line[-1]]
      else:
         refseq_dic[line[1]]=[line[0],int(line[4]),line[-1]]
   return refseq_dic

### output {chrom:[[TSS,refseq_ID,gene_ID,log(FC)]]}
def read_expression(expression_file,refseq_file,header,upregulate=True):
   refseq_dic=read_anotation(refseq_file)
   file=open(expression_file)
   file.readline()
   output_list={}
   i=0
   while True:
      line=file.readline()
      i+=1
      if (not line) or (i==header+1):break
      line=line.strip().split("\t")

      refseq_ID=line[1][1:-4]
      logFC=float(line[2])

      if upregulate:
         if logFC<0:
            continue
      else:
         if logFC>0:
            continue

      try:
         refseq_information=refseq_dic[refseq_ID] ##[chrom,TSS,gene_ID]
      except KeyError:
         print refseq_ID
         continue

      try:
         output_list[refseq_information[0]].append([refseq_information[1],refseq_ID,refseq_information[2],logFC])
         ### output {chrom:[[TSS,refseq_ID,gene_ID,log(FC)]]}
      except KeyError:
         output_list[refseq_information[0]]=[[refseq_information[1],refseq_ID,refseq_information[2],logFC]]

   print "finish expression file reading\n"
   return output_list


##rapid searching
def searchend_for_boundary(sequence,number,lower,upper):
    if lower==upper-1:
        if sequence[lower]<=number<sequence[upper]: return lower
        elif number<sequence[lower]: return lower-1   #if data<all
        else: return upper
    else:
        middle=(lower+upper)/2
        if number>sequence[middle]:
            return searchend_for_boundary(sequence,number,middle,upper)
        else:
            return searchend_for_boundary(sequence,number,lower,middle)

def searchend_for_peak(sequence,number,lower,upper):
    if lower==upper-1:
        if sequence[lower][1]<=number<sequence[upper][1]: return lower
        elif number<sequence[lower][1]: return lower-1   #if data<all
        else: return upper
    else:
        middle=(lower+upper)/2
        if number>sequence[middle][1]:
            return searchend_for_peak(sequence,number,middle,upper)
        else:
            return searchend_for_peak(sequence,number,lower,middle)

##open wig
def open_wigs(wig):
  open_wig=[]
  for i in range(len(wig)):
    open_wig.append(BigWigFile(open(wig[i], 'rb')))
  return open_wig

def profile(open_wig,chrom,interval):
  profile=[]
  for i in range(len(open_wig)):
    try:
      summary=open_wig[i].summarize(chrom,interval[0],interval[1],1)
    except:
      print chrom,interval[0],interval[1]

    if not summary:
      profile.append(0.0)
    else:
      profile.append(float((summary.sum_data / summary.valid_count)[0]))
  return profile

def close_wigs(open_wig):
  for i in range(len(open_wig)):
    open_wig[i].close()

###gene_regulation_dic (refseq_ID,gene_ID):[[logFC,peak,distance,boundary_number,similarity]]
def normalization(gene_regulation_dic):
   ###for Eucluid distance
   similarity=[]
   for value in dictionary.iteritems():
      for infor_value in value[1]:
        similarity.append(float(infor_value[4]))

   similarity.sort(reverse=True)
   similarity_dic={}

   for i in range(len(similarity)):
      similarity_dic[similarity[i]]=float(i+1)/len(similarity)   

   for key in gene_regulation_dic.keys():
      for value in gene_regulation_dic[key]:
        value[4]=similarity_dic[value[4]]

   del similarity_dic
   del similarity

   print "normalization has been finished"
      

def main():
   boundary_dic=read_boundary_score(options.boundary)
   ### {chrom:[0,1000,5000...]}

   peak_dic=read_peak(options.peak)
   ## {chrom:[[peak_ID, summit, signal]]}

   expression_dic=read_expression(options.expression,options.refseq,options.header)
   ### {chrom:[[TSS,refseq_ID,gene_ID,log(FC)]]}

   open_wig=open_wigs(options.wig)
   
   chrom_list=['chr'+v for v in map(str,range(1,30))+['X','Y']]

   

   gene_regulation_dic={} ###(refseq_ID,gene_ID):[[logFC,peak,distance,boundary_number,similarity]]
   for chrom in chrom_list:
      if chrom!="chr1":
         break
      boundary=boundary_dic[chrom]
      peak=peak_dic[chrom]
      expression=expression_dic[chrom]

      for expression_value in expression:
         #print expression_value

         refseq_ID=expression_value[1]
         gene_ID=expression_value[2]
         TSS=expression_value[0]
         logFC=expression_value[3]

         gene_regulation_dic[(refseq_ID,gene_ID)]=[]

         peak_start=searchend_for_peak(peak,TSS-100000,0,len(peak)-1)+1
         peak_end=searchend_for_peak(peak,TSS+100000,0,len(peak)-1)

         TSS_boundary=searchend_for_boundary(boundary,TSS,0,len(boundary)-1)
         TSS_profile=profile(open_wig,chrom,[TSS-1000,TSS+1000])

         for peak_index in range(peak_start,peak_end+1):
            peak_ID=peak[peak_index][0]
            summit=peak[peak_index][1]
            summit_boundary=searchend_for_boundary(boundary,summit,0,len(boundary)-1)
            if summit <= TSS:
               boundary_number=TSS_boundary-summit_boundary
            else:
               boundary_number=summit_boundary-TSS_boundary

            distance=abs(TSS-summit)/float(100000)

            peak_profile=profile(open_wig,chrom,[summit-1000,summit+1000])
            similarity_score=spatial.distance.euclidean(TSS_profile,peak_profile)

            gene_regulation_dic[(refseq_ID,gene_ID)].append([logFC,peak_ID,distance,boundary_number,similarity_score])

          gene_regulation_dic[(refseq_ID,gene_ID)].sort(key=lambda x:x[2],reverse=True)

  close_wigs(open_wig)
  normalization(gene_regulation_dic)
  print gene_regulation_dic

main()

         
            
         
         
         





         
   
