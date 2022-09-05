# -*- coding: utf-8 -*-
import argparse
import os
from multiprocessing import Process

parser = argparse.ArgumentParser(add_help = False, usage = '\npython3 region_meth_ratio.py -g [gtf] -m [meth] \n\n')
required = parser.add_argument_group('必选项')
optional = parser.add_argument_group('可选项')
required.add_argument('-g', '--gtf', metavar = '[gtf]', help = 'gtf promotor文件', required = True)
required.add_argument('-m', '--meth', metavar = '[meth]', help = '输入文件,meth文件前缀', required = True)
optional.add_argument('-c', '--chr', metavar = '[chr]', help = '染色体名称')
optional.add_argument('-h', '--help', action = 'help', help = '帮助信息')
args = parser.parse_args()

#读取depth文件
gtf_file_name=args.gtf
meth_prefix=args.meth
chr_names=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','X','Y','MT']#'MT'
if args.chr!=None:
    chr_names=args.chr.split(',')

output_dir=os.path.join(os.path.dirname(__file__), 'region4_meth_ratio_output')
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

###################################################################################################################
def cal_ratio(arr):
    total=arr[0]+arr[1]
    if total!=0:
        return arr[0]/total
    return 0

def cal(gene_obj,areas):
    chr=gene_obj['chr']
    gene=gene_obj['gene']
    start=gene_obj['start']
    end=gene_obj['end']
    strand=gene_obj['strand']
    rtype=gene_obj['type']

    all=[0,0]
    chh=[0,0]
    chg=[0,0]
    cpg=[0,0]
    point_num=0
    chh_point_num=0
    chg_point_num=0
    cg_point_num=0

    #使用start，end计算起始结束区域
    start_area=int(start/1000)
    end_area=int(end/1000)+1
    for area in range(start_area,end_area):
        #如果区域没在区域列表（即该区域无位点）则直接下个区域
        if area not in areas.keys():
            continue
        #获取到区域，遍历区域内每个点
        area_arr=areas[area]
        for obj in area_arr:
            obj_pos=obj[0]
            obj_strand=obj[1]
            obj_meth=obj[2]
            obj_not_meth=obj[3]
            obj_context=obj[4]

            #位点在start，end的区间内，并且链类型相同
            if obj_strand==strand and (obj_pos>=start and obj_pos<=end):
                #start，end的区间内的点需要满足，测序深度大于等于3，有一个点不满足则该基因不满足条件，直接返回空
                # if obj_meth+obj_not_meth<3:
                #     return None

                point_num+=1
                all[0]+=obj_meth
                all[1]+=obj_not_meth

                if obj_context=='CHH':
                    chh_point_num+=1
                    chh[0]+=obj_meth
                    chh[1]+=obj_not_meth

                if obj_context=='CHG':
                    chg_point_num+=1
                    chg[0]+=obj_meth
                    chg[1]+=obj_not_meth

                if obj_context=='CG':
                    cg_point_num+=1
                    cpg[0]+=obj_meth
                    cpg[1]+=obj_not_meth
    
    # 没有一个点落在start，end区间则直接返回空
    if point_num==0:
        return None
     
    mr=str(cal_ratio(all))
    chh_mr=str(cal_ratio(chh))
    chg_mr=str(cal_ratio(chg))
    cpg_mr=str(cal_ratio(cpg))

        
    return '%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\n' % (chr,gene,start,end,strand,rtype,all[0],all[1],cpg[0],cpg[1],mr,chh_mr,chg_mr,cpg_mr,point_num,chh_point_num,chg_point_num,cg_point_num)

def deal_one_chr(chr_name,region_file_name,meth_file_name,output_file_name):
    genes=[]
    chr_areas={}
    # 初始化areas，遍历所有基因，获取start，end，每1000个点为一个区域，这样初始化数据时只需要要看区域key是否在这个范围内决定是否要这个点，减少不需要的信息计算存储
    with open(region_file_name,'r') as input_file:
        for line in input_file:
            if line.startswith('#'):
                continue
            line=line.rstrip()
            line_arr=line.split()

            chr=line_arr[0]
            start=int(line_arr[2])
            end=int(line_arr[3])
            strand=line_arr[4]

            #计算启动子区域
            # start,end=cal_promotor_region(start,end,strand)

            if chr!=chr_name:
                continue

            genes.append({'chr':chr,'gene':line_arr[1],'start':start,'end':end,'strand':strand,'type':line_arr[5]})

            start_area=int(start/1000)
            end_area=int(end/1000)+1
            for area in range(start_area,end_area):
                if area not in chr_areas.keys():
                    chr_areas[area]=[]
    # 初始化areas内数据,遍历甲基化文件的所有点，若在区域列表内则留下，否则舍去
    with open(meth_file_name,'r') as input_file:
        for line in input_file:
            line=line.rstrip()
            line_arr=line.split()

            pos=int(line_arr[1])
            strand=line_arr[2]
            meth=int(line_arr[3])
            not_meth=int(line_arr[4])
            context=line_arr[5]
            area=int(pos/1000)
            if area not in chr_areas.keys():
                continue
            line_obj=[pos,strand,meth,not_meth,context]
            chr_areas[area].append(line_obj)

    # 计算比例
    result_lines=['chr\tgene\tstart\tend\tstrand\ttype\tmeth\tnot_meth\tcpg_meth\tcpg_not_meth\tmr\tchh_mr\tchg_mr\tcpg_mr\tpoint_num\tchh_point_num\tchg_point_num\tcg_point_num\n']
    for gene in genes:
        cal_result=cal(gene,chr_areas)
        if cal_result:
            result_lines.append( cal_result )
    with open(output_file_name,'w') as output_file:
        output_file.writelines(result_lines) 

# 每个染色体一个进程运行
process_list=[]
for chr_name in chr_names:
    region_file_name=os.path.join(os.path.dirname(__file__), gtf_file_name)
    meth_file_name=os.path.join(os.path.dirname(__file__), meth_prefix+chr_name+'.CX_report.txt')
    output_file_name=os.path.join(output_dir, chr_name+'.txt')

    process=Process(target=deal_one_chr,args=(chr_name,region_file_name,meth_file_name,output_file_name))
    process_list.append(process)
    process.start()
for process in process_list:
    process.join()

# 组合结果
all_result_lines=['chr\tgene\tstart\tend\tstrand\ttype\tmeth\tnot_meth\tcpg_meth\tcpg_not_meth\tmr\tchh_mr\tchg_mr\tcpg_mr\tpoint_num\tchh_point_num\tchg_point_num\tcg_point_num\n']
for chr_name in chr_names:
    chr_result_file_name=os.path.join(output_dir, chr_name+'.txt')
    with open(chr_result_file_name,'r') as input_file:
        num=0
        for line in input_file:
            if num==0:
                num+=1
                continue
            all_result_lines.append(line)

with open(os.path.join(output_dir, 'all.txt'),'w') as output_file:
    output_file.writelines(all_result_lines)    
