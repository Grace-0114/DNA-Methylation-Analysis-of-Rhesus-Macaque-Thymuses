# -*- coding: utf-8 -*-
import argparse
import os
from turtle import right

parser = argparse.ArgumentParser(add_help = False, usage = '\npython3 asm_cluster.py -i [input] -o [output] \n\n')
required = parser.add_argument_group('必选项')
optional = parser.add_argument_group('可选项')
required.add_argument('-i', '--input', metavar = '[input]', help = 'asm文件', required = True)
required.add_argument('-o', '--output', metavar = '[output]', help = '输出文件', required = True)

optional.add_argument('-h', '--help', action = 'help', help = '帮助信息')
args = parser.parse_args()

#参数定义
input=args.input
output=args.output

in_file_path=os.path.join(os.path.dirname(__file__),input)
out_file_path=os.path.join(os.path.dirname(__file__),output)


#提取Clustered S-ASM，每个ASM点+50bp内出现5个asm点称作为Clustered S-ASM
chr_asm_points={}
with open(in_file_path) as input_file:
    for line in input_file:
        arr=line.strip().split('\t')
        if len(arr)<12:
            continue
        chr=arr[0]
        pos=int(arr[1])
        is_asm=arr[11]
        if is_asm=='TRUE':
            if chr not in chr_asm_points.keys():
                chr_asm_points[chr]=[]
            chr_asm_points[chr].append(pos)

chr_asm_clusters={}
for chr,asm_points in chr_asm_points.items():
    print(chr)
    print(asm_points)

    asm_clusters=[]
    left=0
    length=len(asm_points)
    while left<length:
        start_pos=asm_points[left]

        num=1
        right=left+1
        last_vaild_index=-1
        end_pos=-1
        while right<length:
            pos=asm_points[right]
            if pos-start_pos<50:
                last_vaild_index=right
                end_pos=pos
                num+=1
                right+=1
                
            else:
                break
        #满足超过5个点时，设置left为下一个起始位置;否则left设置为下一个点
        if num>=5:
            asm_clusters.append({'start':start_pos,'end':end_pos,'num':num})
            #最后的点满足条件则该点加一，否则则设置为该点
            if last_vaild_index==right:
                left=right+1
            else:
                left=right
        else:
            left+=1
    
    chr_asm_clusters[chr]=asm_clusters

    
print(chr_asm_clusters)

#提取aDMR，当Clustered S-ASM距离小于10bp，合并成aDMR，并计算长度
lines=['chr\tstart\tend\tlen\tnum\tcombined_line\n']
for chr,asm_clusters in chr_asm_clusters.items():
    cur_start=0
    cur_end=0
    cur_num=0
    combine_line=0
    for asm_cluster in asm_clusters:
        start=asm_cluster['start']
        end=asm_cluster['end']
        num=asm_cluster['num']
        #如果当前start为0表明读第一行，直接赋值后返回
        if cur_start==0:
            cur_start=start
            cur_end=end
            cur_num=num
            combine_line=1
            continue

        #当前end与该行start相距小于10，则设置当前end为该行end，并累加num，合并行数加1；否则将当前数据加入结果，并用该行赋值当前
        if start-cur_end<10:
            cur_end=end
            cur_num+=num
            combine_line+=1
        else:
            lines.append('%s\t%d\t%d\t%d\t%d\t%d\n' % (chr,cur_start,cur_end,cur_end-cur_start,cur_num,combine_line) )
            cur_start=start
            cur_end=end
            cur_num=num
            combine_line=1

    lines.append('%s\t%d\t%d\t%d\t%d\t%d\n' % (chr,cur_start,cur_end,cur_end-cur_start,cur_num,combine_line) )

with open(out_file_path,'w') as output_file:
    output_file.writelines(lines)
