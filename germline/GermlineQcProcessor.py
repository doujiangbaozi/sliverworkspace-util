#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
import os,sys,collections,pandas,numpy,getopt,logging

__author__ = '豆浆包子'
__all__    = ['GermlineQcProcessor', 'process', 'doProcess','usage']


class GermlineQcProcessor(object):
    """ 
    根据以下软件输出结果，生成最终Qc文件
    samtools flagstat
    samtools depth
    gatk     CollectInsertSizeMetrics
    每个软件操作sample和matched的bam文件输出，共计3个文件
    """


    def __init__(self):
        '''
        初始化:验证工具,获取脚本参数,Usage使用说明等
        '''
        logging.basicConfig(
            level   = logging.INFO,
            format  = '%(asctime)s - %(filename)s[line:%(lineno)d] - %(levelname)s: %(message)s'
        )
        self.log = logging.getLogger(__name__)

    
        #测序深度分组统计
        self.cover_array = [
            1,
            4,
            10,
            20,
            30,
            100,
            200,
            300
        ]
        self.bed                = None
        self.out                = None
        self.sample_depth       = None
        self.sample_flagstat    = None
        self.sample_insertsize  = None

        try:
            self.opts = getopt.getopt(
                sys.argv[1:], 
                "b:o:dh", [
                    "bed=",
                    "out=",
                    "sample-depth=",
                    "sample-flagstat=",
                    "sample-insertsize=",
                    "document",
                    "help"
                ]
            )
        except getopt.GetoptError:
            print('错误：请检查您的输入参数是否正确\n\n')
            self.usage()
            sys.exit()
        
        if len(self.opts[0])==0:
            self.usage()
            sys.exit()

        #获取命令行参数值
        for arg, value in self.opts[0]:
            if arg=="-b" or arg == "--bed":
                if value.endswith('.bed') and os.path.isfile(value) and os.path.getsize(value)>0:
                    self.bed = value
                else:
                    print('[-b]或[--bed]参数值 %s 不是一个有效的文件(以%s结尾)' %(value,'.bed'))
            elif arg=="-o" or arg == "--out":
                if value.endswith('.xls'):
                    self.out = value
                    if os.path.isfile(value) :
                        self.log.warn('[-o]或[--out]参数值 %s 文件已经存在，输出结果将覆盖该文件' %(value))
                else:
                    print('[-o]或[--out]参数值 %s 不是一个有效的文件(以%s结尾)' %(value,'.xls'))
            elif arg == "--sample-depth":
                if os.path.isfile(value):
                    self.sample_depth=value
                else:
                    print('[--sample-depth]参数值 %s 不是一个有效的文件' % value)
                    sys.exit()
            elif arg == "--sample-flagstat":
                if os.path.isfile(value):
                    self.sample_flagstat=value
                else:
                    print('[--sample-flagstat]参数值 %s 不是一个有效的文件' % value)
                    sys.exit()
            elif arg == "--sample-insertsize":
                if os.path.isfile(value):
                    self.sample_insertsize=value
                else:
                    print('[--sample-insertsize]参数值 %s 不是一个有效的文件' % value)
                    sys.exit()
            elif arg=="-h" or arg == "--help":
                self.usage()
                sys.exit()
            elif arg=="-d" or arg=="--document":
                import GermlineQcProcessor
                help(GermlineQcProcessor)
                sys.exit()
        ps_stats = False
        if self.bed is None:
            print('[-b]或[--bed]参数值不能为空')
            ps_stats = True
        if self.out is None:
            print('[-o]或[--out]参数值不能为空')
            ps_stats = True
        if self.sample_depth is None:
            print('[--sample_depth]参数值不能为空')
            ps_stats = True
        if self.sample_flagstat is None:
            print('[--sample_flagstat]参数值不能为空')
            ps_stats = True
        if self.sample_insertsize is None:
            print('[--sample_insertsize]参数值不能为空')
            ps_stats = True
        if ps_stats:
            sys.exit()

        self.cpath = os.getcwd()
        self.log.info('WorkingDir is : %s' % self.cpath)


    def doProcess(self):
        """
        处理传入的vcf文件和annovar注释文件，输出结果
        args:
            sample_depth       samtools depth 输出文件
            sample_flagstat    samtools flagstat 输出文件
            sample_insertsize  gatk     CollectInsertSizeMetrics 输出文件
        return:
            无返回值
        """

        if os.path.exists(self.bed) and (os.path.getsize(self.bed)>0):

            data = pandas.read_csv(
                self.bed,
                names=['CHROM','START','END'],
                sep='\t',
                header=None,
                skip_blank_lines=True,
                comment='#')
            data['LEN'] = data.apply(lambda row:self.calculateLength(row['START'],row['END']),axis=1)
            tt          = data['LEN'].sum()


            sample = collections.OrderedDict()
            sample['name']='SAMPLE'
            sample = self.collectSamtoolsFlagstat(self.sample_flagstat,sample)
            sample = self.generateResultRow(self.sample_depth,  tt,sample)
            sample = self.collectInsertSizeMetrics(self.sample_insertsize,sample)

            result = pandas.DataFrame([sample],index=[0])
            result.drop(['covered','total'],axis=1,inplace=True)
            self.log.info('writting result to file %s',self.out)
            print(result)
            result.to_csv(self.out,index=False,encoding='utf-8',sep='\t')            
        else:
            self.log.info('writting result to file %s',self.out)
            emptyFile = open(self.out,'w')
            emptyFile.write('')
            emptyFile.close()

    def calculateLength(self,start,end):
        return end-start+1

    def generateResultRow(self,filename,total,dic):
        if os.path.isfile(filename):
            temp_data = pandas.read_csv(
                filename,
                names=['CHROM','POSI','DEPTH'],
                sep='\t',
                header=None,
                skip_blank_lines=True,
                comment='#'
            )
            dic['covered']        = temp_data.shape[0]
            dic['total']          = total
            dic['depth mean']     = int(temp_data['DEPTH'].mean())
            dic['depth median']   = int(temp_data['DEPTH'].median())
            dic['total coverage'] = format(dic['covered']/total, '.2%')

            for i in self.cover_array:
                property = str(i)+'x coverage'
                dic[property]=format((temp_data[temp_data['DEPTH']>i].shape[0])/dic['total'],'.2%')
        return dic

    def collectInsertSizeMetrics(self,filename,dic):
        if os.path.isfile(filename):
            temp_data = pandas.read_csv(
                filename,
                names=[
                    'MEDIAN_INSERT_SIZE',
                    'MODE_INSERT_SIZE',
                    'MEDIAN_ABSOLUTE_DEVIATION',
                    'MIN_INSERT_SIZE',
                    'MAX_INSERT_SIZE',
                    'MEAN_INSERT_SIZE',
                    'STANDARD_DEVIATION',
                    'READ_PAIRS',
                    'PAIR_ORIENTATION',
                    'WIDTH_OF_10_PERCENT',
                    'WIDTH_OF_20_PERCENT',
                    'WIDTH_OF_30_PERCENT',
                    'WIDTH_OF_40_PERCENT',
                    'WIDTH_OF_50_PERCENT',
                    'WIDTH_OF_60_PERCENT',
                    'WIDTH_OF_70_PERCENT',
                    'WIDTH_OF_80_PERCENT',
                    'WIDTH_OF_90_PERCENT',
                    'WIDTH_OF_95_PERCENT',
                    'WIDTH_OF_99_PERCENT',
                    'SAMPLE',
                    'LIBRARY',
                    'READ_GROUP'
                ],
                sep='\t',
                header=0,
                nrows=1,
                comment='#'
            )
            dic['mean insert size']=int(temp_data['MEAN_INSERT_SIZE'][0])
        return dic
    
    def collectDupPercent(self,filename,dic):
        if os.path.isfile(filename):
            temp_data = pandas.read_csv(
                filename,
                names=[
                    'LIBRARY',
                    'UNPAIRED_READS_EXAMINED',
                    'READ_PAIRS_EXAMINED',
                    'SECONDARY_OR_SUPPLEMENTARY_RDS',
                    'UNMAPPED_READS',
                    'UNPAIRED_READ_DUPLICATES',
                    'READ_PAIR_DUPLICATES',
                    'READ_PAIR_OPTICAL_DUPLICATES',
                    'PERCENT_DUPLICATION',
                    'ESTIMATED_LIBRARY_SIZE'
                ],
                sep='\t',
                header=0,
                nrows=1,
                comment='#'
            )
            temp_data['PERCENT_DUPLICATION']=temp_data['PERCENT_DUPLICATION'].apply(lambda x:format(x, '.2%'))
            dic['duplication percent'] = temp_data['PERCENT_DUPLICATION'][0]
        return dic

    def collectSamtoolsFlagstat(self,filename,dic):
        if os.path.isfile(filename):
            temp_file                  = open(filename)
            temp_lines                 = temp_file.readlines()
            total_reads                = int(temp_lines[5].split(' ')[0])
            mapped_reads               = int(temp_lines[4].split(' ')[0])
            duped_reads                = int(temp_lines[3].split(' ')[0])
            all_mapped_reads           = int(temp_lines[0].split(' ')[0])
            dic['total reads(M)']      = format(total_reads/1000000,'.2f')
            dic['mapped reads(M)']     = format(mapped_reads/1000000,'.2f')
            dic['mapped rate']         = format(mapped_reads/all_mapped_reads,'.2%')
            dic['duplication percent'] = format(duped_reads/mapped_reads,'.2%')
        return dic

    def usage(self):
        '''
        打印输出Usage
        '''
        print('Usage : ./GermlineQcProcessor [OPTION]... 或者 python GermlineQcProcessor [OPTION]')
        print('''
            根据vcf文件，annovar注释输出结果，输出最终分析结果文件，扩展名为.xls
            Example:
                ./GermlineQcProcessor.py  -b /opt/ref/Illumina_pt2.bed \\
                    -o /opt/result/2019.result.QC.xls \\
                    --sample-depth=/opt/result/2019_sorted.depth \\
                    --sample-flagstat=/opt/result/2019_bam.stat \\
                    --sample-insertsize=/opt/result/2019_insertsize_metrics.txt
            或者:
                ./GermlineQcProcessor.py  --bed /opt/ref/Illumina_pt2.bed \\
                    --out /opt/result/2019.result.QC.xls \\
                    --sample-depth=/opt/result/2019_sorted.depth \\
                    --sample-flagstat=/opt/result/2019_bam.stat \\
                    --sample-insertsize=/opt/result/2019_insertsize_metrics.txt
        ''')
        print('''部分使用方法及参数如下：\n''')
        print('-b, --bed=\t\t输入bed文件')
        print('-o, --out=\t\t输出处理结果文件')
        print('sample-depth=\t\tSamtools depth 分析bam文件的输出文件')
        print('sample-flagstat=\tSamtools flagstat 分析bam文件的输出文件')
        print('sample-insertsize=\tgatk CollectInsertSizeMetrics 统计bam文件的输出文件')
        print('-h, --help\t显示帮助')
        print('-d, --document\t显示开发文档')
        print('\n')
        print('提交bug,to <6041738@qq.com>. 网站:https://sliverworkspace.com \n')


if __name__ == '__main__':
    f = GermlineQcProcessor()
    f.doProcess()