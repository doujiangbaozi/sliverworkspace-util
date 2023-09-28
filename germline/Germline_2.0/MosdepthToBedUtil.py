#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division

__author__ = '豆浆包子'
__all__ = ['MosdepthToBedUtil', 'process', 'doProcess', 'usage']

import os
import sys
import re
import pandas
import numpy
import getopt
import logging


class MosdepthToBedUtil(object):
    """ 
        输入mosdepth计算bam文件获取的 SRR9993255.per-base.bed.gz 文件反向推算并输出bed文件
        mosdepth --threads 32 \
            /opt/result/SRR9993255/mosdepth/SRR9993255 \
            /opt/result/SRR9993255/align/SRR9993255_marked.bam
    """

    def __init__(self):
        '''
            初始化:验证工具,获取脚本参数,Usage使用说明等
        '''
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(filename)s[line:%(lineno)d] - %(levelname)s: %(message)s'
        )
        self.log       = logging.getLogger(__name__)

        self.input     = None
        self.out       = None
        self.chrs      = [
            'chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12',
            'chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY',
            'chrM'
        ]

        try:
            self.opts = getopt.getopt(
                sys.argv[1:],
                "i:o:dh",
                ["input=", "out=", "document", "help"]
            )
        except getopt.GetoptError:
            print('错误：请检查您的输入参数是否正确\n\n')
            self.usage()
            sys.exit()

        if len(self.opts[0]) == 0:
            self.usage()
            sys.exit()

        # 获取命令行参数值a
        for arg, value in self.opts[0]:
            if arg == "-i" or arg == "--input":
                if (value.endswith('.bed') or value.endswith('.bed.gz')) and os.path.isfile(value):
                    self.input = value
                else:
                    print('[-i]或[--input]参数值 %s 不是一个有效的文件(以%s结尾)' % (value, '.bed 或.bed.gz'))
                    sys.exit()
            elif arg == "-o" or arg == "--out":
                if value.endswith('.bed'):
                    self.out = value
                    if os.path.isfile(value):
                        self.log.warn('[-o]或[--out]参数值 %s 文件已经存在，输出结果将覆盖该文件' % (value))
                else:
                    print('[-o]或[--out]参数值 %s 不是一个有效的文件(以%s结尾)' % (value, '.bed'))
            elif arg == "-h" or arg == "--help":
                self.usage()
                sys.exit()
            elif arg == "-d" or arg == "--document":
                import MosdepthToBedUtil
                help(MosdepthToBedUtil)
                sys.exit()
        ps_stats = False
        if self.input is None:
            print('[-i]或[--input]参数值不能为空')
            ps_stats = True
        if self.out is None:
            print('[-o]或[--out]参数值不能为空')
            ps_stats = True
        if ps_stats:
            sys.exit()
        self.cpath = os.getcwd()
        self.log.info('WorkingDir is : %s' % self.cpath)

    def process(self):
        """输入mosdepth计算bam文件获取的 SRR9993255.per-base.bed.gz 文件反向推算并输出bed文件"""
        self.doProcess(self.input, self.out)

    def doProcess(self, ipt, out):
        """
        输入mosdepth计算bam文件获取的 SRR9993255.per-base.bed.gz 文件反向推算并输出bed文件
        args:
            ipt 输入mosdepth计算bam文件获取的 SRR9993255.per-base.bed.gz文件
            out 输出bed文件
        return:
            无返回值
        """
        if (os.path.exists(ipt)) and (os.path.getsize(ipt) > 0):
            if ipt.endswith('.bed.gz'):
                self.input_data = pandas.read_csv(
                    ipt,
                    names=['CHROM', 'START', 'END', 'DEPTH'],
                    sep='\t',
                    header=None,
                    compression='gzip',
                    skip_blank_lines=True,
                    comment='#')
            else:
                self.input_data = pandas.read_csv(
                    ipt,
                    names=['CHROM', 'START', 'END', 'DEPTH'],
                    sep='\t',
                    header=None,
                    skip_blank_lines=True,
                    comment='#')

            
            #获取DEPTH >0结果
            self.input_data = self.input_data[self.input_data['DEPTH'] > 0]
            self.input_data = self.input_data[self.input_data['CHROM'].apply(lambda x: self.rmOverChrom(x))]
            #self.input_data.sort_values(by=['CHROM', 'START'], inplace=True, ascending=[True, True])
            '''
            for index, row in self.input_data.iterrows():
                if row['DEPTH']=='NAN':
                    continue
                self.merge(index,row,self.input_data)
            
            self.input_data = self.input_data[self.input_data['DEPTH'] != 'NAN']
            '''
            #输出匹配结果
            if not self.input_data.empty:
                print(self.input_data)
                self.log.info('writting result to file %s', out)
                self.input_data.to_csv(out, index=False, encoding='utf-8', sep='\t')
            else:
                self.log.info('writting empty result to file %s', out)
                emptyFile = open(out, 'w')
                emptyFile.write('')
                emptyFile.close()
        else:
            self.log.info('writting empty result to file %s', out)
            emptyFile = open(out, 'w')
            emptyFile.write('')
            emptyFile.close()
    

    def merge(self,index,row,data):
        '''
            重新合并计算坐标
        '''
        try:
            for i, temp in data.iterrows():
                if index >= i:
                    continue
                else:
                    if  temp['DEPTH']=='NAN':
                        continue
                    elif  row['CHROM']==temp['CHROM'] and row['END']>=(temp['START']-1):
                        row['END']   =temp['END']
                        temp['DEPTH']='NAN'
                        continue
                    elif row['CHROM']!=temp['CHROM'] or (row['CHROM']==temp['CHROM'] and row['END']<(temp['START']-1)):
                        break


        except Exception as e:
            print(e)

    def rmOverChrom(self,chr):
        for c in self.chrs:
            if chr==c:
                return True
        return False

    def usage(self):
        '''
        打印输出Usage
        '''
        print(
            'Usage : ./MosdepthToBedUtil [OPTION]... 或者 python MosdepthToBedUtil [OPTION]')
        print('''
            根据vcf文件，annotation注释输出结果，输出最终分析结果文件，格式为.xls
            Example:\t
                ./MosdepthToBedUtil.py \\
                    -i /opt/result/SRR9993255/mosdepth/SRR9993255.per-base.bed.gz \\
                    -o /opt/result/SRR9993256.calculated.bed
            或者\t
                python MosdepthToBedUtil.py \\
                    --input=/opt/result/SRR9993255/mosdepth/SRR9993255.per-base.bed.gz \\
                    --out=/opt/result/SRR9993256.calculated.bed
        ''')
        print('''部分使用方法及参数如下：\n''')
        print('-i, --input=\t\t输入bed或bed.gz格式文件')
        print('-o, --out=\t\t输出处理过的结果文件')
        print('-h, --help\t\t显示帮助')
        print('-d, --document\t\t显示开发文档')
        print('\n')
        print('提交bug,to <6041738@qq.com>. 网站:https://sliverworkspace.com \n')


if __name__ == '__main__':
    f = MosdepthToBedUtil()
    f.process()