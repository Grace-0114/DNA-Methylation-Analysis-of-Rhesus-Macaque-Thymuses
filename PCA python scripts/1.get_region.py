import os,re,random
gtf_file_path=os.path.join(os.path.dirname(__file__),'Mmul_8.0.1.92.chr.gtf')
with open(gtf_file_path,'r') as gtf_file:
    result=[]
    promotor_map={}
    last_exon={}
    for line in gtf_file:
        if line.startswith('#'):
            continue
        arr=line.split('\t')
        chr=arr[0]
        rtype=arr[2]
        start=int(arr[3])
        end=int(arr[4])
        strand=arr[6]
        temp=arr[8]
        gene_id = re.match(r'gene_id "(.*?)";', temp).group(1)

        if rtype=='gene':
            pass
        elif rtype=='exon':
            transcript_id = re.match(r'.*?transcript_id "(.*?)";', temp).group(1)
            if last_exon and last_exon['strand']==strand and last_exon['transcript_id']==transcript_id:
                intron_start=0
                intron_end=0
                if strand=='-':
                    intron_start=end
                    intron_end=last_exon['start']
                elif strand=='+':
                    intron_start=last_exon['end']
                    intron_end=start
                    
                result.append('\t'.join((chr,gene_id,str(intron_start+1),str(intron_end-1),strand,'intron'))+'\n')
            last_exon={'strand':strand,'transcript_id':transcript_id,'start':start,'end':end}
        elif rtype=='transcript':
            if gene_id in promotor_map.keys():
                continue
            rtype='promotor'
            if strand=='+':
                end=start
                start=start-1000
            elif strand=='-':
                start=end
                end=end+1000
            promotor_map[gene_id]=True
        else:
            continue
        result.append('\t'.join((chr,gene_id,str(start),str(end),strand,rtype))+'\n')



region_file_path=os.path.join(os.path.dirname(__file__),'output/1.region/region.txt')
with open(region_file_path,'w') as region_file:
    region_file.writelines(result)          


for sample in ('YA09','YB09','AA09','AB09'):
    arr=[]
    for item in result:
        arr.append(item.replace('\n','')+ '\t%f\t%f\t%f\t%f\t0\t0\t0\t0\n' % (random.random(),random.random(),random.random(),random.random()) )
    meth_file_path=os.path.join(os.path.dirname(__file__),'output/2.meth/'+sample+'.txt')
    with open(meth_file_path,'w') as meth_file:
        meth_file.writelines(arr)       
        
