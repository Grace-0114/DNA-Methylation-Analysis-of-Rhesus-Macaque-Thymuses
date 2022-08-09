
# -*- coding: utf-8 -*-
import argparse
import os

parser = argparse.ArgumentParser(add_help = False, usage = '\npython3 deal_vcf.py -i [input] -s [snv] \n\n')
required = parser.add_argument_group('必选项')
optional = parser.add_argument_group('可选项')
required.add_argument('-i', '--input', metavar = '[input]', help = '输入文件', required = True)
required.add_argument('-s', '--snv', metavar = '[snv]', help = 'snv文件', required = True)
required.add_argument('-a', '--sample', metavar = '[sample]', help = 'sample', required = True)
optional.add_argument('-h', '--help', action = 'help', help = '帮助信息')
args = parser.parse_args()

#读取depth文件
input_file_name=args.input
snv_file_name=args.snv
sample_name=args.sample

input_file_name=os.path.join(os.path.dirname(__file__), input_file_name)
snv_file_name=os.path.join(os.path.dirname(__file__), snv_file_name)

output_dir=os.path.join(os.path.dirname(__file__), os.path.join('vcf_output',sample_name))
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

output_file_name_result=os.path.join(output_dir, 'all.txt') 
output_file_name_statistic=os.path.join(output_dir, 'statistic.txt')

class SnvFilter():
    def __init__(self,file_name):
        super().__init__()

        self.filter_map={}
        with open(file_name,'r') as input_file:
            result_lines=[]
            self.all_total=0
            self.total=0
            for line in input_file:
                if line.startswith('#'):
                    continue
                self.all_total+=1
                line=line.rstrip()
                line_arr=line.split()
                key='chr:%s,nuc:%s,pos:%s' % ( str(line_arr[0]),str(line_arr[1]),str(line_arr[2]) )
                p_value=float(line_arr[-1])
                if p_value<0.000000001 and line_arr[-2].find('/') == -1:
                    self.filter_map[key]=True
                    self.total+=1

    def valid(self,key):
        if key in self.filter_map.keys():
            return True
        return False

snv_filter=SnvFilter(snv_file_name)


def valid(el_str):
    el_arr=el_str.split(',')
    for el in el_arr:
        if el not in ('A','T','C','G'):
            return False
    return True

all_total=0
total=0
dw=0
ljyx=0
bdw=0
chz=0
zhz=0
zhs=0
dhs=0
zdb=0
snv3=0
not_snv=0
with open(input_file_name,'r') as input_file:
    chr_result_lines_map={}
    result_lines=[]
    comment_lines=[]
    
    for line in input_file:
        if line.startswith('#'):
            result_lines.append(line)
            comment_lines.append(line)
            continue

        all_total+=1


        line_arr=line.split('\t')
        NA00001=line_arr[9]
        REF=line_arr[3]
        ALT=line_arr[4]
        alt_arr=ALT.split(',')

        chrom=str(line_arr[0])
        pos=line_arr[1]
        nuc=line_arr[3]
        snv_key='chr:%s,nuc:%s,pos:%s' % ( str(chrom),str(nuc),str(pos) )
        if not snv_filter.valid(snv_key):
            continue

        if not valid(REF)  or not valid(ALT):
            continue
        if not int(NA00001.split(':')[2])>=15:
            continue

        #等位基因
        if len(ALT)==1 and REF != ALT:
            dw+=1
            if (not ( NA00001.startswith('0/0') or NA00001.startswith('1/1') ) ) :
                result_lines.append(line)
                if chrom not in chr_result_lines_map:
                    chr_result_lines_map[chrom]=[]
                chr_result_lines_map[chrom].append(line)
        #两基因型位点
        val_arr=NA00001.split(':')[0].split('/')
        if val_arr[0]!=val_arr[1]:
            ljyx+=1
        else:
            bdw+=1
        #纯合子 杂合子
        # if len(REF)==1 and len(ALT)==1 and REF==ALT:
        #     chz+=1
        # else:
        #     zhz+=1
        if NA00001.startswith('0/0') or NA00001.startswith('1/1'):
            chz+=1
        else:
            zhz+=1
            
        #转换数
        if (REF=='G' and ALT=='C') or (REF=='C' and ALT=='G') or (REF=='A' and ALT=='T') or (REF=='T' and ALT=='A'):
            zhs+=1
        #颠换数
        if (REF=='G' and (ALT=='A' or ALT=='T') ) or \
            (REF=='C' and (ALT=='A' or ALT=='T') ) or \
            (REF=='A' and (ALT=='G' or ALT=='C') ) or \
            (REF=='T' and (ALT=='G' or ALT=='C') ) :
            dhs+=1
        #SNV三种
        snv_arr=[REF]
        if len(alt_arr)==2:
            if 'R' not in alt_arr and 'Y' not in alt_arr:
                snv_arr+=alt_arr
        if len(set(snv_arr))==3:
            snv3+=1
        #未准确预测SNV数
        if 'R' in alt_arr or 'Y' in alt_arr:
            not_snv+=1

        #总数
        total+=1

    for key,val in chr_result_lines_map.items():
        val=comment_lines+val

if dhs!=0:
    zdb=round(zhs/float(dhs)*100,2)

with open(output_file_name_result,'w') as output_file:
    output_file.writelines(result_lines)

for key,val in chr_result_lines_map.items():
    chr_output_file_name_result=os.path.join(output_dir, key+'.txt') 
    with open(chr_output_file_name_result,'w') as output_file:
        output_file.writelines(val)


with open(output_file_name_statistic,'w') as output_file:
    # lines=['snv总数\tsnv总数(过滤后)\t总数\t总数(过滤后)\t纯合子\t杂合子\t转换数\t颠换数\t转换/颠换比\t两基因型位点\t二等位+三等位\t二等位位点\tSNV三种\t未准确预测SNV数\n']
    # lines.append('%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%d' % (snv_filter.all_total,snv_filter.total,all_total,total,chz,zhz,zhs,dhs,str(zdb)+'%',ljyx,(dw+snv3),dw,snv3,not_snv) )
    lines=['总数\t总数(过滤后)\t纯合子\t杂合子\t转换数\t颠换数\t转换/颠换比\t二等位+三等位\t二等位位点\tSNV三种\t未准确预测SNV数\n']
    lines.append('%d\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d' % (all_total,total,chz,zhz,zhs,dhs,str(zdb)+'%',(dw+snv3),dw,snv3,not_snv) )
    output_file.writelines(lines)
