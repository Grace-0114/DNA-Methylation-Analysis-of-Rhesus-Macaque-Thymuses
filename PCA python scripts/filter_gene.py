import os,re,random

gene_name_id_map={}
with open(os.path.join(os.path.dirname(__file__),'2.txt'),'r') as file2:
    for line in file2:
        arr=line.strip().split('\t')
        gene_id=arr[0]
        gene_name=gene_id
        if(len(arr)==2):
            gene_name=arr[1]
        gene_name_id_map[gene_name]=gene_id


result_lines=[]
with open(os.path.join(os.path.dirname(__file__),'1.txt'),'r') as file1:
    for line in file1:
        arr=line.strip().split('\t')
        chr=arr[0]
        strand=arr[1]
        gene_name=arr[2]
        gene_id=gene_name
        if gene_name in gene_name_id_map.keys():
            gene_id=gene_name_id_map[gene_name]
        else:
            print(gene_name)
        result_lines.append( '%s\t%s\t%s\t%s\n' % (chr,strand,gene_name,gene_id)  )


with open(os.path.join(os.path.dirname(__file__),'3.txt'),'w') as output_file:
    output_file.writelines(result_lines)          


   
        
