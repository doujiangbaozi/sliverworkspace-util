#!/usr/bin/python
# -*-coding:utf-8-*-
from __future__ import division
__author__ = '豆浆包子'
__all__ = ['RnaSeqSnvAnnotationFilter', 'process', 'doProcess', 'usage']

import os
import sys
import re
import pandas
import numpy
import getopt
import logging


class RnaSeqSnvAnnotationFilter(object):
    """ 输入vcf文件，和annovar注释文件,输出注释结果 """

    def __init__(self):
        '''
        初始化:验证工具,获取脚本参数,Usage使用说明等
        '''
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(filename)s[line:%(lineno)d] - %(levelname)s: %(message)s'
        )
        self.log      = logging.getLogger(__name__)

        self.vcf      = None
        self.out      = None
        self.anno     = None
        self.depth    = None

        self.vtype = [
            'synonymous SNV', 
            'NAN'
        ]

        try:
            self.opts = getopt.getopt(
                sys.argv[1:],
                "a:o:v:i:e:dh",
                ["annovar=", "out=", "vcf=","min-depth=", "document", "help"]
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
            elif arg == "-a" or arg == "--annovar":
                if os.path.isfile(value):
                    self.anno = value
                else:
                    print('[-a]或[--annovar]参数值 %s 不是一个有效的Annovar注释输出文件' % (value))
                    sys.exit()
            elif arg == "-o" or arg == "--out":
                if value.endswith('.xls'):
                    self.out = value
                    if os.path.isfile(value):
                        self.log.warn(
                            '[-o]或[--out]参数值 %s 文件已经存在，输出结果将覆盖该文件' % (value))
                else:
                    print('[-o]或[--out]参数值 %s 不是一个有效的文件(以%s结尾)' % (value, '.xls'))
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
                import RnaSeqSnvAnnotationFilter
                help(RnaSeqSnvAnnotationFilter)
                sys.exit()
        ps_stats = False
        if self.vcf is None:
            print('[-v]或[--vcf]参数值不能为空')
            ps_stats = True
        if self.out is None:
            print('[-o]或[--out]参数值不能为空')
            ps_stats = True
        if self.anno is None:
            print('[-a]或[--annovar]参数值不能为空')
            ps_stats = True
        if ps_stats:
            sys.exit()
        self.cpath = os.getcwd()
        self.log.info('WorkingDir is : %s' % self.cpath)

    def process(self):
        """处理传入的vcf文件和annovar注释文件，输出结果"""
        self.doProcess(self.vcf, self.anno, self.out)

    def doProcess(self, vcf, anno, out):
        """
        处理传入的vcf文件和annovar注释文件，输出结果
        args:
            vcf    vcf文件
            anno   annovar注释结果文件
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

            self.vcf_data['START'] = self.vcf_data.POS
            self.vcf_data['DP']    = self.vcf_data.apply(lambda row: self.calculateDP(row), axis=1)
            self.vcf_data['DP']    = self.vcf_data['DP'].astype('int64')

            self.vcf_data['VAF']   = self.vcf_data.apply(lambda row: self.calculateVAF(row), axis=1)
            self.vcf_data['VAF']   = self.vcf_data['VAF'].astype('float')
            self.vcf_data['VAF']   = self.vcf_data['VAF'].round(decimals=5)

            self.vcf_data = self.vcf_data[['CHROM', 'START', 'REF', 'ALT', 'DP', 'VAF', 'FILTER']]
            print(self.vcf_data)

            anno_data = pandas.read_csv(
                anno,
                header=0,
                names=[
                    'Chr','Start','End','Ref','Alt','Func.refGene','Gene.refGene','GeneDetail.refGene','ExonicFunc.refGene',
                    'AAChange.refGene','avsnp150','esp6500siv2_all','1000g2015aug_all','CLNALLELEID',
                    'CLNDN','CLNDISDB','CLNREVSTAT','CLNSIG','cosmic89_coding'
                ],
                sep='\t',
                skip_blank_lines=True,
                comment='#'
            )

            result_data = pandas.merge(self.vcf_data, anno_data, how='inner', left_index=True, right_index=True)
            print('merged size:',len(result_data))
            if not result_data.empty:
                result_data = result_data.apply(lambda row: self.calculateTrans(row), axis=1)
                result_data['Chr'] = result_data.apply(lambda row: self.calculateChrNum(row['Chr']), axis=1)
                result_data['COSMIC_ID']        = result_data.apply(lambda row:self.calculateCosmicID(row),axis=1)
                result_data['COSMIC_OCCURENCE'] = result_data.apply(lambda row:self.calculateCosmicOCCURENCE(row),axis=1)
                index = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'DP', 'VAF', 'FILTER', 'Transcript', 'Exon', 'cHGVS', 'pHGVS',
                        'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene',
                        'AAChange.refGene', 'avsnp150', 'CLNALLELEID',
                        'CLNDN', 'CLNDISDB', 'CLNREVSTAT', 'CLNSIG','COSMIC_ID','COSMIC_OCCURENCE']
                result_data = result_data.reindex(columns=index)
                #result_data = result_data.reset_index(drop=True)
                result_data['VAF'] = result_data.VAF.apply(lambda x: format(x, '.2%'))
                print('reindex size:',len(result_data))
                result_data = result_data[result_data['FILTER']=='PASS']
                print('filter pass size:',len(result_data))
                if self.depth is not None:
                    result_data = result_data[result_data['DP'] >= self.depth]
                    print('depth pass size:',len(result_data))

                # 过滤掉同义点不变和NAN突变类型
                for tp in self.vtype:
                    result_data = result_data[result_data['ExonicFunc.refGene'] != tp]
                print('varition size:',len(result_data))
                
                #根据CLNSIG字段过滤
                result_data['CLNSIG_STATUS'] = result_data.apply(lambda row:self.filterCLINIG(row['CLNSIG']),axis=1)
                result_data = result_data[result_data['CLNSIG_STATUS']==True]
                result_data.drop(['CLNSIG_STATUS','FILTER'],axis=1,inplace=True)
                print('CLNSIG size:',len(result_data))

                self.log.info('writting result to file %s', out)
                result_data.to_csv(out, index=False, encoding='utf-8', sep='\t')
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
        except Exception, e:
            print(e)
            return 'NAN'

    def calculateDP(self, row):
        '''根据row/行数据，获取INFO字段中的DP值，即测序深度值'''
        try:
            info_fields = row['INFO'].split(';')
            for field in info_fields:
                if field.find('DP=') != -1:
                    return field[field.find('DP=')+len('DP='):]
        except Exception, e:
            print(e)
        return '0'
    
    def calculateCosmicID(self,row):
        '''根据row/行数据，获取cosmic89_coding字段中的ID值'''
        try:
            cos_fields = row['cosmic89_coding'].split(';')
            for field in cos_fields:
                if field.find('ID=') != -1:
                    return field[field.find('ID=')+len('ID='):]
                else:
                    return 'NAN'
        except Exception,e:
            print(e)
            return 'NAN'

    def calculateCosmicOCCURENCE(self,row):
        '''根据row/行数据，获取cosmic89_coding字段中的OCCURENCE值'''
        try:
            cos_fields = row['cosmic89_coding'].split(';')
            for field in cos_fields:
                if field.find('OCCURENCE=') != -1:
                    return field[field.find('OCCURENCE=')+len('OCCURENCE='):]
            return 'NAN'
        except Exception,e:
            self.log.error(e)
            return 'NAN'
    
    def calculateChrNum(self, chr):
        '''根据染色体编号chrx计算后面纯数字字母x'''
        try:
            field = chr.lower()
            if field.find('chr') != -1:
                return field[field.find('chr')+len('chr'):].upper()
        except Exception, e:
            print(e)
            return 'NAN'

    def calculateTrans(self, row):
        '''根据Annovar注释信息中的AAChange.refGene字段拆分出信息'''
        try:
            array_fields = row['AAChange.refGene'].split(',')
            if len(array_fields) > 0:
                ps_fields = array_fields[0].split(':')
                if len(ps_fields) == 5:
                    row['Transcript'] = ps_fields[1]
                    row['Exon']       = ps_fields[2]
                    row['cHGVS']      = ps_fields[3]
                    row['pHGVS']      = ps_fields[4]
        except Exception, e:
            print(e)
        return row

    def filterCLINIG(self,CLNSIG):
      try:
          if CLNSIG=='NAN':
              return False
          elif CLNSIG=='Benign':
              return False
          elif CLNSIG=='Likely_benign':
              return False
          elif CLNSIG=='Benign/Likely_benign':
              return False
          elif CLNSIG=='not_provided':
              return False
          elif CLNSIG=='not_specified':
              return False
          elif CLNSIG=='Uncertain_significance':
              return False
          else:
              return True
      except Exception, e:
          print(e)
          return False
    
    def usage(self):
        '''
        打印输出Usage
        '''
        print(
            'Usage : ./RnaSeqSnvAnnotationFilter [OPTION]... 或者 python RnaSeqSnvAnnotationFilter [OPTION]')
        print('''
            根据vcf文件，annovar注释输出结果，输出最终分析结果文件，格式为.xls
            Example:\t
                ./RnaSeqSnvAnnotationFilter.py \\
                    -v result/B1701_filtered_snpeff.vcf \\
                    -a result/B1701.hg19_multianno.txt \\
                    -o result/B1701.result.xls
            或者\t
                python RnaSeqSnvAnnotationFilter.py \\
                    --vcf=/opt/result/B1701_filtered_snpeff.vcf \\
                    --annovar=/opt/result/B1701.hg19_multianno.txt \\
                    --out=/opt/result/B1701.result.xls
        ''')
        print('''部分使用方法及参数如下：\n''')
        print('-v, --vcf=\t\t输入vcf或vcf.gz格式文件')
        print('-a, --annovar=\t\t输入annovar注释后的文件')
        print('-o, --out=\t\t输出处理过的结果文件')
        print('    --min-depth=\t最小测序深度该突变位点最小reads数（可选）')
        print('    --database=\t\t突变过滤数据文件，如耳聋基因及突变位点')
        print('-h, --help\t\t显示帮助')
        print('-d, --document\t\t显示开发文档')
        print('\n')
        print('提交bug,to <6041738@qq.com>.\n')


if __name__ == '__main__':
    f = RnaSeqSnvAnnotationFilter()
    f.process()
