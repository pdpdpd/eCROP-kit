'''

import ssl
from functools import wraps
def sslwrap(func):
    @wraps(func)
    def bar(*args, **kw):
        kw['ssl_version'] = ssl.PROTOCOL_TLSv1
        return func(*args, **kw)
    return bar

ssl.wrap_socket = sslwrap(ssl.wrap_socket)


'''
#---------------------------------import---------------------------------------

import urllib2
import math
import time
import random
import socket
import re
import sys, csv, operator
import ssl


#socket.setdefaulttimeout(20)


from threading import Timer

#------------------------------------------------------------------------------
#02.10 fixed some bugs of starting point and ending point. Improve the robustness.
#06.07 use this script to find the two-mutations
#08 for SNP
#08.22 correct the SNP localizing function

#12062018 Try to make this script suitable for anyone to use

'''
    SNP gRNA generating script
    Step:  1  remove the title, sort the raw data based on focal SNP. Only keep the first 5 (?) columns
    Store to a new txt
    Create a txt for final data. Add title.
    For Sorted table:
    1 Read the first line
    2 Store the focal SNP to a variable
    3 Check the gRNA of focal SNP
    Get Chr location based on rsID (there's no position in the table)
    Get the window (+-60)sequence
    Create Focal SNPID.txt, add title
    Sense: find NGG in the range 20-60, store in matrix[0][0] (or array.)
                Get the position, store the 30-(pos-2) in matrix[0][2]
                [0][1]="+"
                Get to COSMID.
                Use the script I wrote for counting: on target[0][3], 1mm[0][4],2 mm[0][5]
                If |[0][2]|>10 or [0][2]>1, [0][6]=0 else =1
                append to the txt
                Repeat till no more NGG
                Antisense: find NGG
                Complementary reverse the seq
                Then the same as sense.(except for [n][1]="-")
                Store all of [0][6]=1 to the final data.
                Format: Focal, control, array [0]->[6]
                Count: open the txt, check [n][6]
                If all =0, DO not NEED TO CHECK THE CONTROLS of IT AT ALL..
                MEWMEWMEWMEWMWEMMMMMEEEEWWWW
                SKIP all the lines with first column = focal SNPID
                MEWEMEWEMEWMEWEWMEWEMWEMWE
                Else,
                For line  with first column = focal SNPID, do the same thing in 
                3 Check the gRNA of focal SNP
'''

def handler(fh):
    fh.close()


################################################
def read_by_line(url): # read by line
    lines = url.splitlines()
    return lines
###########################################

################################################
def my_sq_range(start, end, step):
    while start < end:
        yield start
        start += step
###########################################

################################################
def SNP_position(seq):#get the pos of given SNP in hg19
    userMainUrl = "http://genome.ucsc.edu/cgi-bin/hgTracks?hgsid=506882255_Jr5qJGeYD0g3gRD2xc2vj6D99Zwv&org=Human&db=hg19&position="+seq+"&pix=1256"
    req = urllib2.Request(userMainUrl)
    resp = urllib2.urlopen(req)
    respHtml = resp.read()
    
    _gotcha=0
    pos=respHtml.find(' at chr')
    if pos != -1:
    #print "There is no "+seq+" in hg19."
        end_pos=respHtml.find('</A>')
        #print end_pos
        start_pos=pos+4
        chr_pos_SNP=respHtml[start_pos:end_pos]
        _gotcha=1
        return chr_pos_SNP
    else:
        print "There is no "+seq+" in hg19 database. Please check your RSID."
        return ''

#rs2590375 at chr10:59967228-59967728</A>
#search:" at chr"

#    print "respHtml=",respHtml # you should see the ouput html

###############################################################################


#######
def flanking_seq(chr_,snp_):#find the sequence for gRNA design
    
    # true location now

    
    try:
        start_= int(snp_)-30+250
        end_=int(snp_)+30+250#now the position is not at the begining of where the SNP starts
    except:
        print "Wrong pos"
    else:
        url_="http://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment="+chr_+"%3A"+str(start_)+","+str(end_)

        #print url_
        req_ = urllib2.Request(url_)
        resp_ = urllib2.urlopen(req_)
        respHtml_ = resp_.read()
        DNA_seq=respHtml_[respHtml_.find('</DNA>')-63:respHtml_.find('</DNA>')]
        DNA_=DNA_seq.replace('\n','')
        #print "DNA_"+DNA_
    return DNA_

#######


#########
def search_gRNA(seq,strand):#find gRNAs that can be sent to COSMID
    #print seq
    pos=[]
    guides=[]
    output=[]
    #print len(seq)
    if strand == '-':
        seq=''.join(["atcg"["tagc".index(n)] for n in seq[::-1]])#reverse complementary
        #print seq
    pos = [m.start() for m in re.finditer('(?=gg)', seq[20:len(seq)])]

    if len(pos) == 0:
        print "No matches in strand "+strand+" of "+seq
        return []
    else:
        #print pos
        for i in range (0,len(pos)):
            #print "i"
            #print i
            
            
# cutting site: xxx(here) xxx N[the position given]GG
            pos[i]+=20
            sub_seq=seq[pos[i]-20:pos[i]-1]


#add GC content,
#can add quanta"X"
            gc_rate=100*(sub_seq.count("g")+ sub_seq.count("c"))/len(sub_seq)
            print gc_rate
            if gc_rate >=40 and gc_rate <=80 and abs(31-(pos[i]-4))<=10: # and sub_seq.find("aaaa")==sub_seq.find("tttt")==sub_seq.find("cccc")==sub_seq.find("gggg")==-1: 

#depends
           
              if pos[i] >= 35:
                output = [sub_seq,-31+(pos[i]-4),strand]
              else:
                output = [sub_seq,31-(pos[i]-4),strand]
              guides.append(output)

#https://www.abmgood.com/marketing/knowledge_base/CRISPR_Cas9_gRNA_Design.php
        
#2019HappynewYear

        return guides


'''
referred from http://stackoverflow.com/questions/4664850/find-all-occurrences-of-a-substring-in-python
'''
#########

################################################
def COSMID(seq): #for COSMID checking, I used the sequence gRNA+NGG (NGG is appended by COSMID) because gRNA itself does not include NGG
    #print "sleeping"
    slp=random.random()
    time.sleep(slp*5)
    #print "sleeping done"
    
    print "searching "+seq+" now......."
    ssl._create_default_https_context = ssl._create_unverified_context
    userMainUrl = "https://crispr.bme.gatech.edu/cgi-bin/crispr/CRISPER_form.cgi?target_database=hg19&tag_type=seq&tag=N"+seq+"&tag_suffix=NGG&CheckBox_no_indel=on&mismatch_no_indel=2&mismatch_1_del=1&mismatch_1_ins=1&minSeparationUncleavedToCleaved=110&minCleavageProductSizeDifference=0&minAmpliconLength=220&maxAmpliconLength=330&optimalAmpliconLength=275"
    #print userMainUrl
    req = urllib2.Request(userMainUrl)
    resp=urllib2.urlopen(req)
    print req

    #print "here"
    
    
    t = Timer(5, handler,[resp])#we only need 1 perfect match site -- so short time.
    t.start()
        
    try:
        respHtml = resp.read()
    except Exception,e:
        respHtml = "no_data"
        pass
    t.cancel()
    
    #print "hhhhere"

    return respHtml
################################################


################################################
def get_info(html):#and write in a file:3
    #output_file="flanking_new.txt"#name of the output
    output=[]
    if html == "no_data":
        return['','','','']
    else:

        #outfile.write(search_seq+"  ")
        lines = read_by_line(html)
        length = len(lines)

        on_target=0
        one_mm=0
        two_mm=0
        total_=0

        for i in my_sq_range(0, length-4, 1):
            try:
                lines[0]
            except:
                print lines[i]
            else:
                if lines[i]== "<td nowrap>No indel</td>" and lines[i+1]== "<td nowrap><center>0</center></td>" and lines[i+2] == "<td nowrap><center>Yes</center></td>":
                    on_target+=1
                elif lines[i]== "<td nowrap>No indel</td>" and lines[i+1]== "<td nowrap><center>1</center></td>" and lines[i+2] == "<td nowrap><center>Yes</center></td>":
                    one_mm+=1
                elif lines[i]== "<td nowrap>No indel</td>" and lines[i+1]== "<td nowrap><center>2</center></td>" and lines[i+2] == "<td nowrap><center>Yes</center></td>":
                    two_mm+=1

        total_=on_target+one_mm+two_mm
        #outfile.write(str(on_target)+"   "+str(one_mm)+"  "+str(two_mm)+"  "+str(total_)+"\n")
        output=[on_target,one_mm,two_mm,total_]
        print output
        return output
         
###############################################################################




'''
#sort:
unsorted_data = csv.reader(open('Target_test1.csv'),delimiter=',')
sorted_data = sorted(unsorted_data, key = lambda x: (x[0], int(x[1])))

with open("Sorted.csv", "w", newline = '') as f:
    fileWriter = csv.writer(f, delimiter=',')
    for row in sorted_data:
        fileWriter.writerow(row)
'''

#0 get file info
#SNP_file = sys.argv[0]
SNP_file="SNP1.csv"

#1 Get the list of SNPs
SNP_list=[]
Cant_find_SNP=[]
output_SNP=[]
output_SNP_mm=[]#The ones do not fit the mm requirements or distance requirements (but maybe acceptable)


with open(SNP_file) as g:
    g_csv = csv.DictReader(g)
    for row in g_csv:
        print row['RSID']
        SNP_list.append(row['RSID'])

####file that contains no_gRNA SNP
file_object = open('no_gRNA_SNP.txt','w')
file_object.write("no_gRNA_SNP\n")
file_object.close( )
####

#cant find the RSID
file_object = open('no_SNP.txt','w')
file_object.close( )

######outputfile######Focal output######
outfile = open('SNP_test.txt', 'w');
outfile.write("RSID    gRNA    SNP-gRNA Distance    Strand  On target    1mm    2mm    total_mm"+"\n");
outfile.close();
######outputfile######Focal output######

for i in range(0,len(SNP_list)):
    print i
    result=[]
    result_filter=[]
               
    chr_pos_SNP= SNP_position(SNP_list[i])
    if len(chr_pos_SNP)==0:
    # didn't find anything based on the RSID
    #maybe store these to a new array? could be useful.
        file_object = open('no_SNP.txt','a')
        file_object.write(SNP_list[i]+"\n")
        file_object.close( )
               
    else:
        print chr_pos_SNP
        end_pos = chr_pos_SNP.find(':')
        chr_pos = chr_pos_SNP[0:end_pos]
        seq_pos = chr_pos_SNP[end_pos+1:chr_pos_SNP.find('-')]
        target_seq = flanking_seq(chr_pos,seq_pos)
        gRNA_seq_sense = search_gRNA(target_seq,'+')
        gRNA_seq_anti = search_gRNA(target_seq,'-')
        gRNA_seq=gRNA_seq_sense+gRNA_seq_anti
        print "All the gRNAs for "+SNP_list[i]+":"
        print gRNA_seq
               
    try:
        len(gRNA_seq[0])
    except:
        Cant_find_SNP.append(SNP_list[i])
        print "No matched gRNA found for " + SNP_list[i]
        file_object = open('no_gRNA_SNP.txt','a')
        file_object.write(SNP_list[i]+"\n")
        file_object.close( )
        pass
    else:
        for j in range(0,len(gRNA_seq)):
            print gRNA_seq[j][0]
            html = COSMID(gRNA_seq[j][0])
            print gRNA_seq[j]
            print get_info(html)
            result.append([SNP_list[i]]+gRNA_seq[j]+get_info(html))
        for j in range (0, len(result)):
            if result[j][2]<=10 and result[j][4]==1 and result[j][5]<=5: #just choose the # I like..
               #result_filter.append(result[j])
               result_filter.append(result[j])
               #print "get output!"
            elif result[j][4] != '':
               output_SNP_mm.append(result[j])
            print "result_filter for "+SNP_list[i]
            print result_filter
               
        if len(result_filter)==0:#no suitable gRNA for this SNP.
            Cant_find_SNP.append(SNP_list[i])
            print "No matched gRNA found for " + SNP_list[i]
            file_object = open('no_gRNA_SNP.txt','a')
            file_object.write(SNP_list[i]+"\n")
            file_object.close( )
        else:
            SNP_list.append(SNP_list[i])
            output_SNP=output_SNP+result_filter
            outfile = open('SNP_test.txt','a')
            for j in range(0,len(result_filter)):
                for k in range(0,len(result_filter[j])):
                    outfile.write(str(result_filter[j][k])+"	")
                    outfile.write("\n")
            outfile.close( )

print "SNP_list"
print SNP_list
print "Cant_find_SNP"
print Cant_find_SNP

print "output_SNP is "
print output_SNP

print "output_SNP_mm"
print output_SNP_mm

######outputfile######Focal output######
outfile = open("SNP_list_test_1003.txt", "w");
outfile.write("Focal_RSID	gRNA	SNP-gRNA Distance	Strand  On target	1mm	2mm	total_mm"+"\n");
for i in range(0,len(output_SNP)):
    for j in range(0,len(output_SNP[i])):
        outfile.write(str(output_SNP[i][j])+"	")
    outfile.write("\n")
outfile.close();
######outputfile######Focal output######


####
#life is beautiful
