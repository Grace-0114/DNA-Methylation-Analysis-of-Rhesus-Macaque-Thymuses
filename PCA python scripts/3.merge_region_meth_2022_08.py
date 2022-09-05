import os
from random import sample
from unittest import result

samples=('AA09','AB09','YA09','YB09')
rtypes=('gene','promotor','intron','exon')

def deal_sample(sample):
    meth_file_path=os.path.join(os.path.dirname(__file__),'samples/'+sample+'.txt')
    gene_map={}
    with open(meth_file_path,'r') as meth_file:
        index=0
        for line in meth_file:
            if index==0:
                index+=1
                continue
            arr=line.strip().split('\t')
            chr=arr[0]
            gene_id=arr[1]
            strand=arr[4]
            rtype=arr[5]

            meth=int(arr[8])
            not_meth=int(arr[9])

            key='%s_%s_%s_%s' % (chr,gene_id,strand,rtype)
            if key not in gene_map.keys():
                gene_map[key]={'chr':chr,'gene_id':gene_id,'strand':strand,'rtype':rtype,'count':0,'meth':0,'not_meth':0}
            gene_map[key]['count']+=1
            gene_map[key]['meth']+=meth
            gene_map[key]['not_meth']+=not_meth
            # print(line)

    result_lines=['chr\tgene_id\tstrand\trtype\tcount\tmeth\tnot_meth\ttotal\tmeth_ratio\n']
    for key,val in gene_map.items():
        val['total']=val['meth']+val['not_meth']
        val['mr']='NA'
        if val['total']!=0:
            val['mr']=val['meth']/val['total']
        result_lines.append( '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (val['chr'],val['gene_id'],val['strand'],val['rtype'],val['count'],val['meth'],val['not_meth'],val['total'],val['mr'],) )
    
    #写中间结果
    dir_path=os.path.join(os.path.dirname(__file__),'result1')
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
    with open(os.path.join(dir_path,sample+'.txt'),'w') as output_file:
        output_file.writelines(result_lines)

    return gene_map

def init_merge_item(item):
    for sample in samples:
        for rtype in rtypes:
            attr_name='%s-%s' % (sample,rtype)
            item[attr_name]='NA'
def gen_merge_header():
    header=''
    for sample in samples:
        for rtype in rtypes:
            attr_name='%s-%s' % (sample,rtype)
            header+='\t'+attr_name
    return header
def gen_merge_row(item):
    s=''
    for sample in samples:
        for rtype in rtypes:
            attr_name='%s-%s' % (sample,rtype)
            s+='\t%s' % item[attr_name]
    row='%s\t%s\t%s%s\n' % (item['chr'],item['gene_id'],item['strand'],s)

    return row

protein_coding_gene_map={}
with open(os.path.join(os.path.dirname(__file__),'protein_coding_gene.txt'),'r') as protein_coding_gene_file:
    for line in protein_coding_gene_file:
        protein_coding_gene_map[line.strip()]=True

filter_gene_map={}
with open(os.path.join(os.path.dirname(__file__),'filter_gene/3.txt'),'r') as filter_gene_file:
    for line in filter_gene_file:
        arr=line.strip().split('\t')
        key='%s_%s_%s' % (arr[0],arr[3],arr[1])
        filter_gene_map[key]=True


def valid(item):
    # for sample in samples:
    #     for rtype in rtypes:
    #         attr_name='%s-%s' % (sample,rtype)
    #         if item[attr_name]=='NA':
    #             return False
    #         # if item[attr_name]==0:
    #         #     return False
    # if item['gene_id'] not in protein_coding_gene_map.keys():
    #     return False

    key='%s_%s_%s' % (item['chr'],item['gene_id'],item['strand'])
    if key not in filter_gene_map.keys():
        return False


    return True
    
merge_map={}
for sample in samples:
    gene_map=deal_sample(sample)
    for key,val in gene_map.items():
        merge_key='%s_%s_%s' % (val['chr'],val['gene_id'],val['strand'])
        if merge_key not in merge_map.keys():
            merge_map[merge_key]={'chr':val['chr'],'gene_id':val['gene_id'],'strand':val['strand']}
            init_merge_item(merge_map[merge_key])
        attr_name='%s-%s' % (sample,val['rtype'])
        merge_map[merge_key][attr_name]=val['mr']
    
result_lines=['chr\tgene_id\tstrand%s\n' % gen_merge_header() ]
for key,val in merge_map.items():
    if valid(val):
        result_lines.append( gen_merge_row(val) )

dir_path=os.path.join(os.path.dirname(__file__),'result2')
if not os.path.exists(dir_path):
    os.makedirs(dir_path)
with open(os.path.join(dir_path,'result.txt'),'w') as output_file:
    output_file.writelines(result_lines)