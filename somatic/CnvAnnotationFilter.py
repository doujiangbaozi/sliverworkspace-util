#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = '豆浆包子'
__all__    = ['CnvAnnotationFilter', 'filter', 'doFilter','usage']

import os,sys,re,pandas,numpy,getopt,logging

class CnvAnnotationFilter(object):
    """ 过滤vcf文件,输出Filter Annotation """

    def __init__(self):
        '''
        初始化:验证工具,获取脚本参数,Usage使用说明等
        '''
        logging.basicConfig(
            level   = logging.INFO,
            format  = '%(asctime)s - %(filename)s[line:%(lineno)d] - %(levelname)s: %(message)s'
        )
        self.log        = logging.getLogger(__name__)
        self.ref        = None
        self.extension  = ".call.cns"
        
        self.refGene    = None #'/opt/ref/hg19_regGene.txt'
        self.min        = -0.5 #-0.5 默认值
        self.max        =  0.5 # 0.5 默认值
        self.depth      =  500 # 500 默认值
        self.file       = None #cnvkit 输出文件B1701_sorted.call.cns
        self.out        = None #输出文件
        self.bed        = None #输出最终结果bed文件，用于图像生成（可选）

        try:
            self.opts = getopt.getopt(
                sys.argv[1:], 
                "f:o:r:i:x:D:b:hv", 
                ["file=","out=","ref=","min=","max=","depth=","bed=","version", "help","document"]
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
            if arg=="-f" or arg == "--file":
                if value.endswith(self.extension) and os.path.isfile(value) :
                    self.file = value
                else:
                    print('[-f]或[--file]参数值 %s 不是一个有效的文件(以%s结尾)' %(value,self.extension))
                    sys.exit()
            elif arg=="-o" or arg == "--out":
                if value.endswith('.tsv'):
                    self.out = value
                    if os.path.isfile(value) :
                        self.log.warn('[-o]或[--out]参数值 %s 文件已经存在，输出结果将覆盖该文件' %(value))
                else:
                    print('[-o]或[--out]参数值 %s 不是一个有效的文件(以%s结尾)' %(value,'.tsv'))
            elif arg=="-b" or arg == "--bed":
                if value.endswith('.bed'):
                    self.bed = value
                    if os.path.isfile(value) :
                        self.log.warn('[-b]或[--bed]参数值 %s 文件已经存在，输出结果将覆盖该文件' %(value))
                else:
                    print('[-o]或[--out]参数值 %s 不是一个有效的文件(以%s结尾)' %(value,'.tsv'))
            elif arg=="-r" or arg == "--ref":
                if os.path.isfile(value):
                    self.refGene = value
                else:
                    print('[-r]或[--ref]参数值 %s 不是一个有效的文件' %(value))
            elif arg=="-i" or arg == "--min":
                pattern = re.compile(r'^[-+]?[0-9]+\.[0-9]+$')
                result = pattern.match(value)
                if value.isdigit() or result:
                    self.min = float(value)
                else:
                    print('[-i]或[--min]参数值 %s 不是一个有效的数值 ' %(value))
            elif arg=="-x" or arg == "--max":
                pattern = re.compile(r'^[-+]?[0-9]+\.[0-9]+$')
                result = pattern.match(value)
                if value.isdigit() or result:
                    self.max = float(value)
                else:
                    print('[-i]或[--min]参数值 %s 不是一个有效的数值 ' %(value))
            elif arg=="-D" or arg == "--depth":
                pattern = re.compile(r'^[-+]?[0-9]+\.[0-9]+$')
                result = pattern.match(value)
                if value.isdigit() or result:
                    self.depth = float(value)
                else:
                    print('[-i]或[--min]参数值 %s 不是一个有效的数值 ' %(value))
            elif arg=="-h" or arg == "--help":
                self.usage()
                sys.exit()
            elif arg=="--document":
                import CnvAnnotationFilter
                help(CnvAnnotationFilter)
                sys.exit()
            elif arg == "--version" or arg=="-v":
                print("版本: 1.00")
                sys.exit()
        status = False
        if self.refGene is None:
            print('[-r]或[--ref]参数值不能为空')
            status = True
        if self.min is None:
            print('[-i]或[--min]参数值不能为空')
            status = True
        if self.max is None:
            print('[-x]或[--max]参数值不能为空')
            status = True
        if self.depth is None:
            print('[-D]或[--depth]参数值不能为空')
            status = True
        if self.file is None:
            print('[-f]或[--file]参数值不能为空')
            status = True
        if self.out is None:
            print('[-o]或[--out]参数值不能为空')
            status = True
        if status:
            sys.exit()
        self.cpath = os.getcwd()
        self.log.info('WorkingDir is : %s' % self.cpath)

    def calculateType(self,log2):
        if log2>0:
            return 'CNV gain'
        elif log2<0:
            return 'CNV loss'

    def calculateCopyNumber(self,log2):
        return 2**(log2+0.5)

    def calculateExon(self,row):
        startIndex = None
        endIndex   = None
        chrom      = row['chromosome']
        gene       = row['gene']
        start      = row['start']
        end        = row['end']
        exonStarts = []
        exonEnds   = []
        temp       = self.ref[(self.ref['chrom']==chrom)&(self.ref['name2']==gene)]
        
        for _index, _row in temp.iterrows():
            exonStarts = list(set(exonStarts)|set(_row['exonStarts'].split(',')))
            exonEnds   = list(set(exonEnds  )|set(_row['exonEnds'  ].split(',')))
            if '' in exonStarts:
                exonStarts.remove('')
            if '' in exonEnds:
                exonEnds.remove('')
            
            _starts = []
            _ends   = []
            
            for ts in exonStarts:
                try:
                    _starts.append(int(ts))
                except ValueError:
                    print('转换错误')
            for te in exonEnds:
                try:
                    _ends.append(int(te))
                except ValueError:
                    print('转换错误')
            
            #_starts.sort()
            #_ends.sort()

            if len(_starts) == len(_ends):        
                for i in range(len(_starts)):
                    if start <= _ends[i]:
                        startIndex = i+1
                        break
                for j in range(len(_starts)): 
                    if end  >_starts[j]:
                        endIndex = j+1
            if (startIndex is not None) and (endIndex is not None):
                return 'Exon '+str(startIndex)+'-'+str(endIndex)
            else:
                startIndex = None
                endIndex   = None
        return 'Exon '+str(startIndex)+'-'+str(endIndex)

    def calculateChrNum(self,chr):
        '''根据染色体编号chrx计算后面纯数字字母x'''
        try:
            field = chr.lower()
            if field.find('chr') != -1:
                return field[field.find('chr')+len('chr'):].upper()
        except Exception as e:
            print(e)
            return 'NAN'

    def filter(self):
        """执行过滤，过滤传入的-f参数文件"""
        self.doFilter(self.file,self.out,self.bed)

    def doFilter(self,filename,out,bed):
        """
        对一个文件名为filename的文件执行过滤,并输出结果
        args:
            filename 一个文件名
            out      输出文件名
        return:
            无返回值
        """
        
        if (os.path.exists(filename)) and (os.path.getsize(filename)>0):
            if self.ref is None:
                self.ref = pandas.read_csv(
                    self.refGene,
                    names=['bin','name','chrom','strand','txStart','txEnd','cdsStart','cdsEnd','exonCount','exonStarts','exonEnds',
                        'score','name2','cdsStartStat','cdsEndStat','exonFrame'],
                    sep='\t',
                    header=None,
                    skip_blank_lines=True,
                    comment='#')
                self.ref.drop(['bin','strand','txStart','txEnd','cdsStart','cdsEnd',
                        'score','cdsStartStat','cdsEndStat','exonFrame'],axis=1,inplace=True)

            
            data = pandas.read_csv(
                filename,
                sep='\t',
                header=0,
                skip_blank_lines=True,
                comment='#')
            
            data   = data[((data['log2']>self.max)|(data['log2']<self.min))&(data['depth']>self.depth)]
            if not data.empty:
                data['type']   = data.apply(lambda row:self.calculateType(row['log2']),axis=1)
                data['region'] = data.apply(lambda row:self.calculateExon(row),axis=1)
                #data['num']    = data.apply(lambda row:self.calculateCopyNumber(row['log2']),axis=1)
                data.drop(['log2','depth','probes','weight'],axis=1,inplace=True)
                #print(data)
                data = data[['chromosome','start','end','gene','region','type','cn']]
                if bed is not None:
                    self.log.info('writting bed to file %s',bed)
                    bed  = data[['chromosome','start','end']].copy()
                    bed.to_csv(self.bed,index=False,encoding='utf-8',sep='\t',header=False)
                data['chromosome']  = data.apply(lambda row:self.calculateChrNum(row['chromosome']),axis=1)
            print(data)
            self.log.info('writting result to file %s',out)
            data.to_csv(out,index=False,encoding='utf-8',float_format = '%.3f',sep='\t')
        else:
            self.log.info('writting empty result to file %s',out)
            emptyFile = open(out,'w')
            emptyFile.write('')
            emptyFile.close()

    def usage(self):
        '''
        打印输出Usage
        '''
        print('Usage : ./CnvAnnotationFilter [OPTION]... 或者 python CnvAnnotationFilter [OPTION]')
        print('''
            用于处理变异注释文件.vcf格式，过滤其中有效变异，与样本编号相符条件，输出格式为.tsv
            Example: 
            ./CnvAnnotationFilter.py  \\
                -i -0.50 \\
                -x 0.50 \\
                -D 1500 \\
                -f result/B1701_sorted.call.cns \\
                -r /opt/ref/hg19_regGene.txt \\
                -o result/B1701.result.Cnv.tsv
        ''')
        print('''部分使用方法及参数如下：\n''')
        print('-f, --file=\t处理一个扩展名为.vcf的文件')
        print('-o, --out=\t处理结果输出文件')
        print('-r, --ref=\t引用参考文件refGene.txt')
        print('-i, --min=\tcutoff值下限')
        print('-x, --max=\tcutoff值上限')
        print('-D, --depth=\tcutoff测序深度')
        #print('-b, --bed=\t输出cnv突变范围bed文件，用于绘图(可选)')
        print('-h, --help\t显示帮助')
        print('-v, --version\t显示版本号')
        print('--document\t显示开发文档')
        print('\n')
        print('提交bug,to <6041738@qq.com>. 网站:https://sliverworkspace.com \n')


if __name__ == '__main__':
    f = CnvAnnotationFilter()
    f.filter()