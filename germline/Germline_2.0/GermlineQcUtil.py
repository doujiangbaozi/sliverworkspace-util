#!/root/mambaforge-pypy3/envs/gatk/bin/python
#-*-coding:utf-8-*-
import os,sys,collections,pandas,numpy,getopt,logging,json

__author__ = '豆浆包子'
__all__    = [
    'GermlineQcUtil', 
    'process', 
    'doProcess',
    'collectFastqInfo',
    'collectBamInfo',
    'collectInsertSize'
    'usage'
]


class GermlineQcUtil(object):
    """ 
        根据以下软件输出结果，生成最终Qc文件
        fastp    fastp.json
        bamdst   coverage.report
        gatk     CollectInsertSizeMetrics(Picard)
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
        
        self.sn                 ='Sample'
        self.out                = None
        self.sample_fastp       = None
        self.sample_bamdst      = None
        self.sample_insertsize  = None


        try:
            self.opts = getopt.getopt(
                sys.argv[1:], 
                "o:dh", 
                [
                    "out=",
                    "sample-fastp=",
                    "sample-bamdst=",
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
            if arg=="-o" or arg == "--out":
                if value.endswith('.xls'):
                    self.out = value
                    if os.path.isfile(value) :
                        self.log.warn('[-o]或[--out]参数值 %s 文件已经存在，输出结果将覆盖该文件' %(value))
                else:
                    print('[-o]或[--out]参数值 %s 不是一个有效的文件(以%s结尾)' %(value,'.xls'))
            elif arg == "--sample-fastp":
                if os.path.isfile(value) and value.endswith('.json'):
                    self.sample_fastp = value
                    fullpath     = os.path.realpath(value)
                    fileOnly     = os.path.basename(fullpath)
                    self.sn      = fileOnly[0:fileOnly.rindex('_fastp.json')]
                else:
                    print('[--sample-fastp]参数值 %s 不是一个有效的文件' % value)
                    sys.exit()
            elif arg == "--sample-bamdst":
                if os.path.isfile(value) :
                    self.sample_bamdst=value
                else:
                    print('[--sample-bamdst]参数值 %s 不是一个有效的文件' % value)
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
                import GermlineQcUtil
                help(GermlineQcUtil)
                sys.exit()
        ps_stats = False
        if self.out is None:
            print('[-o]或[--out]参数值不能为空')
            ps_stats = True
        if self.sample_fastp is None:
            print('[--sample-fastp]参数值不能为空')
        if self.sample_bamdst is None:
            print('[--sample-bamdst]参数值不能为空')
            ps_stats = True
        if self.sample_insertsize is None:
            print('[--sample-insertsize]参数值不能为空')
            ps_stats = True
        if ps_stats:
            sys.exit()

        self.cpath = os.getcwd()
        self.log.info('WorkingDir is : %s' % self.cpath)


    def doProcess(self):
        """
        处理传入的vcf文件和annovar注释文件，输出结果
        args:
            fastp              fastp    输出文件
            sample_bamdst      bamdst   输出文件
            sample_insertsize  gatk     CollectInsertSizeMetrics 输出文件
        return:
            无返回值
        """
        dic                                    = collections.OrderedDict()
        dic['Sample Number']                   = self.sn
        dic['Total Reads (M) Before Filtering']='NAN'
        dic['Q20 Rate Before Filtering']       ='NAN'
        dic['Q30 Rate Before Filtering']       ='NAN'
        dic['GC Rate  Before Filtering']       ='NAN'

        dic['Total Reads (M) After  Filtering']='NAN'
        dic['Q20 Rate After  Filtering']       ='NAN'
        dic['Q30 Rate After  Filtering']       ='NAN'
        dic['GC Rate  After  Filtering']       ='NAN'
        dic['Mapping Rate']                    ='NAN'
        dic['Duplicate Rate']                  ='NAN'
        dic['Target in Mapping Rate']          ='NAN'
        dic['Region Length']                   ='NAN'
        dic['Average Depth']                   ='NAN'
        dic['1× Coverage']                     ='NAN'
        dic['4× Coverage']                     ='NAN'
        dic['10× Coverage']                    ='NAN'
        dic['20× Coverage']                    ='NAN'
        dic['30× Coverage']                    ='NAN'
        dic['100× Coverage']                   ='NAN'
        #dic['500× Coverage']                   ='NAN'
        dic['mean insert size']                ='NAN'
        dic                                    = self.collectFastqInfo(self.sample_fastp,dic)
        dic                                    = self.collectBamInfo(self.sample_bamdst,dic)
        dic                                    = self.collectInsertSize(self.sample_insertsize,dic)
        
        result                                 = pandas.DataFrame([dic],index=[0])
        #result.drop(['covered','total'],axis=1,inplace=True)
        self.log.info('writting result to file %s',self.out)
        #print(result.stack())
        result.to_csv(self.out,index=False,encoding='utf-8',sep='\t')
        

    def collectFastqInfo(self,filename,dic):
        '''
            从fastp输出json文件中读取过滤前后Total Reads Q20比例，Q30比例，GC %百分比
        '''
        if os.path.isfile(filename):
            try:
                with open(filename,'r',encoding='utf-8') as load_f:
                    load_dict = json.load(load_f)
                    dic['Total Reads (M) Before Filtering'] = format(float(load_dict['summary']['before_filtering']['total_reads'])/1000000,'.2f')
                    dic['Q20 Rate Before Filtering'] = format(load_dict['summary']['before_filtering']['q20_rate'],'.2%')
                    dic['Q30 Rate Before Filtering'] = format(load_dict['summary']['before_filtering']['q30_rate'],'.2%')
                    dic['GC Rate  Before Filtering'] = format(load_dict['summary']['before_filtering']['gc_content'],'.2%')
                    dic['Total Reads (M) After  Filtering'] = format(float(load_dict['summary']['after_filtering']['total_reads'])/1000000,'.2f')
                    dic['Q20 Rate After  Filtering'] = format(load_dict['summary']['before_filtering']['q20_rate'],'.2%')
                    dic['Q30 Rate After  Filtering'] = format(load_dict['summary']['before_filtering']['q30_rate'],'.2%')
                    dic['GC Rate  After  Filtering'] = format(load_dict['summary']['before_filtering']['gc_content'],'.2%')
            except Exception as e:
                self.log.error(e)
        return dic

    def collectBamInfo(self,filename,dic):
        '''
            读取bamdst输出文件coverage.report，获取相关字段
            Mapping Rate,Duplicate Rate,Target Reads in Mapping Rate,
            Average Depth,1×,4×,10×,20,30×,100×,500× Coverage@
        '''
        if os.path.isfile(filename):
            f = open(filename,'r',encoding='utf-8')
            line = f.readline()
            while line:
                #print (line)
                line = f.readline()
                if line.startswith('##') or len(line)==0:
                    continue
                arrs = line.split('\t')
                key  = arrs[0].strip()
                val  = arrs[1].split('\n')[0]
                if key=='[Total] Fraction of Mapped Reads':
                    dic['Mapping Rate']=val
                elif key=='[Total] Fraction of PCR duplicate reads':
                    dic['Duplicate Rate']=val
                elif key=='[Target] Fraction of Target Reads in mapped reads':
                    dic['Target in Mapping Rate']=val
                elif key=='[Target] Len of region':
                    dic['Region Length']=val
                elif key=='[Target] Average depth(rmdup)':
                    dic['Average Depth']=val
                elif key=='[Target] Coverage (>0x)':
                    dic['1× Coverage']=val
                elif key=='[Target] Coverage (>=4x)':
                    dic['4× Coverage']=val
                elif key=='[Target] Coverage (>=10x)':
                    dic['10× Coverage']=val
                elif key=='[Target] Coverage (>=20x)':
                    dic['20× Coverage']=val
                elif key=='[Target] Coverage (>=30x)':
                    dic['30× Coverage']=val
                elif key=='[Target] Coverage (>=100x)':
                    dic['100× Coverage']=val
                #elif key=='[Target] Coverage (>=500x)':
                #  dic['500× Coverage']=val
            f.close()
        return dic

    def collectInsertSize(self,filename,dic):
        '''
            gatk CollectInsertSizeMetrics 输出文件获取mean insert size
        '''
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
                comment='#',
                encoding = 'utf-8'
            )
            dic['mean insert size']=int(temp_data['MEAN_INSERT_SIZE'][0])
        return dic
    
    def usage(self):
        '''
        打印输出Usage
        '''
        print('Usage : ./GermlineQcUtil [OPTION]... 或者 python GermlineQcUtil [OPTION]')
        print('''
            根据fastp,bamdst,gatk CollectInsertSizeMetrics(picard)
            输出质控分析结果文件，扩展名为.xls
            Example:
                ./GermlineQcUtil.py \\
                    -o /opt/result/2019.result.QC.xls \\
                    --sample-fastp=/opt/result/2019_fastp.json  \\
                    --sample-bamdst=/opt/result/coverage.report \\
                    --sample-insertsize=/opt/result/2019_insertsize_metrics.txt
            或者:
                python GermlineQcUtil.py \\
                    --out=/opt/result/2019.result.QC.xls \\
                    --sample-fastp=/opt/result/2019_fastp.json  \\
                    --sample-bamdst=/opt/result/coverage.report \\
                    --sample-insertsize=/opt/result/2019_insertsize_metrics.txt
        ''')
        print('''部分使用方法及参数如下：\n''')
        print('-o, --out=\t\t输出处理结果文件')
        print('--sample-fastp=\t\tfastp 处理后的输出文件')
        print('--sample-bamdst=\tbamdst分析bam文件的输出文件')
        print('--sample-insertsize=\tgatk CollectInsertSizeMetrics(Picard) 统计bam文件的输出文件')
        print('-h, --help\t\t显示帮助')
        print('-d, --document\t\t显示开发文档')
        print('\n')
        print('提交bug,to <6041738@qq.com>.\n')


if __name__ == '__main__':
    f = GermlineQcUtil()
    f.doProcess()