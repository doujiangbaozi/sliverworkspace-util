#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division

__author__ = '豆浆包子'
__all__ = ['GermlineVepAnnotationUtil', 'process', 'doProcess', 'usage']

import os
import sys
import re
import pandas
import numpy
import getopt
import logging


class GermlineVepAnnotationUtil(object):
    """ 
        输入vcf文件，和vep annotation注释文件,输出注释结果 
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

        self.vcf       = None
        self.out       = None
        self.anno      = None
        self.whitelist = None
        self.datalist  = None
        self.depth     = None

        try:
            self.opts = getopt.getopt(
                sys.argv[1:],
                "a:o:v:dh",
                ["annotation=", "out=", "vcf=","min-depth=", "document", "help"]
            )
        except getopt.GetoptError:
            print('错误：请检查您的输入参数是否正确\n\n')
            self.usage()
            sys.exit()

        if len(self.opts[0]) == 0:
            self.usage()
            sys.exit()

        # 获取命令行参数值
        for arg, value in self.opts[0]:
            if arg == "-v" or arg == "--vcf":
                if (value.endswith('.vcf') or value.endswith('.vcf.gz')) and os.path.isfile(value):
                    self.vcf = value
                else:
                    print('[-v]或[--vcf]参数值 %s 不是一个有效的文件(以%s结尾)' % (value, '.vcf 或 .vcf.gz'))
                    sys.exit()
            elif arg == "-a" or arg == "--annotation":
                if os.path.isfile(value):
                    self.anno = value
                else:
                    print('[-a]或[--annotation]参数值 %s 不是一个有效的Annovar注释输出文件' % (value))
                    sys.exit()
            elif arg == "-o" or arg == "--out":
                if value.endswith('.tsv'):
                    self.out = value
                    if os.path.isfile(value):
                        self.log.warn(
                            '[-o]或[--out]参数值 %s 文件已经存在，输出结果将覆盖该文件' % (value))
                else:
                    print('[-o]或[--out]参数值 %s 不是一个有效的文件(以%s结尾)' % (value, '.tsv'))
            elif arg == "--min-depth":
                if value.isdigit():
                    self.depth = int(value)
                    self.log.info('--min-depth=%s', value)
                else:
                    print('[--min-depth]参数值 %s 不是一个有效的正整数 ' % (value))
            elif arg == "-h" or arg == "--help":
                self.usage()
                sys.exit()
            elif arg == "-d" or arg == "--document":
                import germline.GermlineVepAnnotationUtil as GermlineVepAnnotationUtil
                help(GermlineVepAnnotationUtil)
                sys.exit()
        ps_stats = False
        if self.vcf is None:
            print('[-v]或[--vcf]参数值不能为空')
            ps_stats = True
        if self.out is None:
            print('[-o]或[--out]参数值不能为空')
            ps_stats = True
        if self.anno is None:
            print('[-a]或[--annotation]参数值不能为空')
            ps_stats = True
        if ps_stats:
            sys.exit()
        self.cpath = os.getcwd()
        self.log.info('WorkingDir is : %s' % self.cpath)

    def process(self):
        """处理传入的vcf文件和annotation注释文件，输出结果"""
        self.doProcess(self.vcf, self.anno, self.out)

    def doProcess(self, vcf, anno, out):
        """
        处理传入的vcf文件和annotation注释文件，输出结果
        args:
            vcf    vcf文件
            anno   vep annotation注释结果文件
            out    输出分析结果文件
        return:
            无返回值
        """
        if (os.path.exists(vcf)) and (os.path.getsize(vcf) > 0) and (os.path.exists(anno)) and (os.path.getsize(anno) > 0):
            if vcf.endswith('.vcf.gz'):
                self.vcf_data = pandas.read_csv(
                    vcf,
                    names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE'],
                    sep='\t',
                    header=None,
                    compression='gzip',
                    skip_blank_lines=True,
                    comment='#')
            else:
                self.vcf_data = pandas.read_csv(
                    vcf,
                    names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL','FILTER', 'INFO', 'FORMAT', 'SAMPLE'],
                    sep='\t',
                    header=None,
                    skip_blank_lines=True,
                    comment='#')

            self.vcf_data['POS']   = self.vcf_data['POS'].astype('int64')
            self.vcf_data['DP']    = self.vcf_data.apply(lambda row: self.calculateDP(row), axis=1)
            self.vcf_data['DP']    = self.vcf_data['DP'].astype('int64')

            self.vcf_data['VAF']   = self.vcf_data.apply(lambda row: self.calculateVAF(row), axis=1)
            self.vcf_data['VAF']   = self.vcf_data['VAF'].astype('float')
            self.vcf_data['VAF']   = self.vcf_data['VAF'].round(decimals=5)

            self.vcf_data          = self.vcf_data.apply(lambda row: self.adjustRefAndAlt(row),axis=1)

            #获取FILTER为PASS的结果
            self.vcf_data          = self.vcf_data[self.vcf_data['FILTER'] == 'PASS']
            print(self.vcf_data)

            cols      = self.processAnnoCols(anno)
            anno_data = pandas.read_csv(
                anno,
                header=0,
                names=cols,
                sep='\t',
                skip_blank_lines=True,
                comment='#'
            )
            #获取CLIN_SIG注释结果含有pathogenic字段，减少后续运算量
            anno_data             = anno_data.apply(lambda row:self.calculateAnnotationPosition(row),axis=1)
            anno_data['Start']    = anno_data['Start'].astype('int64')
            
            self.vcf_data         = self.vcf_data.apply(lambda row: self.processVcfWithAnnotation(row,anno_data),axis=1)
            print(self.vcf_data)
            index = [
                'Chr', 'Start', 'End', 'Ref', 'Alt', 'DP', 'VAF', 'Gene_ID', 'Gene', 'Type', 'Transcript', 'Exon', 'cHGVS', 'pHGVS', 'CLNSIG',
                'Existing_variation', 'BioType', 'PUBMED'
            ]
            self.vcf_data          = self.vcf_data[index]
            self.vcf_data          = self.vcf_data.reset_index(drop=True)
            self.vcf_data['VAF']   = self.vcf_data['VAF'].apply(lambda x: format(x, '.2%'))

            if self.depth is not None:
                self.vcf_data = self.vcf_data[self.vcf_data['DP'] >= self.depth]
            
            self.vcf_data         = self.vcf_data[self.vcf_data['CLNSIG'].apply(lambda x: str(x).find('pathogenic') != -1)]
            self.vcf_data         = self.vcf_data[self.vcf_data['Chr'].notnull()]
            self.vcf_data         = self.vcf_data[self.vcf_data['Start'].notnull()]
            self.vcf_data['Start']= self.vcf_data['Start'].astype('int64')

            #输出匹配结果
            if not self.vcf_data.empty:
                print(self.vcf_data)
                self.log.info('writting result to file %s', out)
                self.vcf_data.to_csv(out, index=False, encoding='utf-8', sep='\t')
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

    def processAnnoCols(self,filename):
        with open(filename,'r') as f:
            for line in f:
                if not line.startswith('##') and line.startswith('#'):
                    header = line[1:len(line)-1].split('\t')
                    print(header)
                    return header
                else:
                    continue
    
    def calculateVAF(self, row):
        '''根据row/行数据计算该行的突变频率'''
        try:
            info_fields = row['INFO'].split(';')
            for field in info_fields:
                if field.find('AF=') != -1:
                    t_vaf = field[field.find('AF=')+len('AF='):]
                    if t_vaf.find(',')!=-1:
                        return t_vaf.split(',')[0]
                    else: 
                        return field[field.find('AF=')+len('AF='):]
        except Exception as e:
            print(e)
            return 'NAN'

    def calculateDP(self, row):
        '''根据row/行数据，获取INFO字段中的DP值，即测序深度值'''
        try:
            info_fields = row['INFO'].split(';')
            for field in info_fields:
                if field.find('DP=') != -1:
                    return field[field.find('DP=')+len('DP='):]
        except Exception as e:
            print(e)
            return 'NAN'

    def calculateChrNum(self, chr):
        '''根据染色体编号chrx计算后面纯数字字母x'''
        try:
            field = chr.lower()
            if field.find('chr') != -1:
                return field[field.find('chr')+len('chr'):].upper()
        except Exception as e:
            print(e)
            return 'NAN'

    def adjustRefAndAlt(self,row):
        '''调整vcf数据中REF和ALT表示方法和vep注释保持一致'''
        try:
            pos = row['POS']
            ref = row['REF']
            alt = row['ALT']            
            if len(ref) > len(alt) and ref.startswith(alt):
                row['REF']=ref[len(alt):]
                row['ALT']='-'
                row['POS']=row['POS']+len(alt)
                #print(row['CHROM'],row['POS'],row['REF'],row['ALT'])
        except Exception as e:
            print(e)
        
        return row

    def matchWhiteList(self, gene):
        if self.datalist is not None:
            temp = self.datalist[self.datalist['GENE'] == gene]
            if len(temp) > 0:
                return True
            else:
                return False
    
    def calculateAnnotationPosition(self,row):
        '''
            vep注释文件中第2列Location拆分成新突变坐标
        '''
        try:
            fields = row['Location'].split(':')
            if len(fields)==2:
                row['Chr']=fields[0]
                pos       =fields[1].split('-')
                if len(pos)==1:
                    row['Start']=pos[0]
                    row['End']  ='-'
                elif len(pos)==2:
                    row['Start']=pos[0]
                    row['End']  =pos[1]
        except Exception as e:
            print(e)
        return row
    
    def processVcfWithAnnotation(self,row,anno):
        '''
            结合vcf文件和vep注释结果
        '''
        try:
            for i in range(len(anno)):
                r = anno.iloc[i]
                if row['CHROM']==r['Chr'] and row['POS']==r['Start']:
                    row['Chr']                = self.calculateChrNum(row['CHROM'])
                    row['Start']              = row['POS']
                    row['Ref']                = row['REF']
                    row['Alt']                = row['ALT']
                    row['Gene']               = r['SYMBOL']
                    row['Type']               = r['VARIANT_CLASS'].upper()
                    row['Transcript']         = r['Feature']
                    row['Exon']               = r['EXON']
                    row['cHGVS']              = r['HGVSc'] # if r['HGVSc'].find(':')==-1 else r['HGVSc'].split(':')[1]
                    row['pHGVS']              = r['HGVSp'] # if r['HGVSp'].find(':')==-1 else r['HGVSp'].split(':')[1]
                    row['Gene_ID']            = r['Gene']
                    row['Existing_variation'] = r['Existing_variation']
                    row['BioType']            = r['BIOTYPE']
                    row['CLNSIG']             = r['CLIN_SIG']
                    row['PUBMED']             = r['PUBMED']
                    row['End']                = self.calculateEndPosition(row)
                    #print(row)
                    return row
        except Exception as e:
            print(e)
            #return row

    def calculateEndPosition(self,row):
        '''
            根据row Start,Ref,Alt Type计算 End坐标
        '''
        try:
            start = row['Start']
            ref   = row['Ref']
            tp    = row['Type']
            if   tp=='DELETION':
                return start+len(ref)-1
            elif tp=='SNV':
                return start+len(ref)-1
            elif tp=='INSERTION':
                return '-'
        except Exception as e:
            print(e)
            return '-'
    
    

    def filterPathogenic(self,row):
        '''
            个过滤注释中pathogenic相关记录
        '''
        return str(row['CLNSIG']).find('pathogenic')!=-1


    def usage(self):
        '''
        打印输出Usage
        '''
        print(
            'Usage : ./GermlineVepAnnotationUtil [OPTION]... 或者 python GermlineVepAnnotationUtil [OPTION]')
        print('''
            根据vcf文件，annotation注释输出结果，输出最终分析结果文件，格式为.tsv
            Example:\t
                ./GermlineVepAnnotationUtil.py \\
                    -v result/SRR9993256_filtered.vcf \\
                    -a result/SRR9993256_filtered_vep.tsv \\
                    -o result/SRR9993256.result.tsv
            或者\t
                python GermlineVepAnnotationUtil.py \\
                    --vcf=/opt/result/SRR9993256_filtered.vcf \\
                    --annotation=/opt/result/SRR9993256_filtered_vep.tsv \\
                    --out=/opt/result/SRR9993256.result.tsv
        ''')
        print('''部分使用方法及参数如下：\n''')
        print('-v, --vcf=\t\t输入vcf或vcf.gz格式文件')
        print('-a, --annotation=\t\t输入annotation注释后的文件')
        print('-o, --out=\t\t输出处理过的结果文件')
        print('    --min-depth=\t最小测序深度该突变位点最小reads数（可选）')
        print('-h, --help\t\t显示帮助')
        print('-d, --document\t\t显示开发文档')
        print('\n')
        print('提交bug,to <6041738@qq.com>. 网站:https://sliverworkspace.com \n')


if __name__ == '__main__':
    f = GermlineVepAnnotationUtil()
    f.process()
