#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = '豆浆包子'
__all__    = [
    'VcfProcessUtil', 
    'process', 
    'doProcessHeader',
    'doProcess',
    'toChr',
    'usage'
]

import io,os,logging,time,sys,getopt,gzip,pandas

class VcfProcessUtil(object):
    """ 
        
    """

    def __init__(self):
        '''
            初始化:验证工具,获取脚本参数,Usage使用说明等
        '''
        logging.basicConfig(
            level   = logging.INFO,
            format  = '%(asctime)s - %(filename)s[line:%(lineno)d] - %(levelname)s: %(message)s'
        )
        self.log    = logging.getLogger(__name__)

        self.data   = None
        self.header = None
        self.skip   = 0

        try:
            self.opts = getopt.getopt(sys.argv[1:], "f:o:hv", ["file=", "output=", "help", "version", "document"])
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
                if (value.endswith(".vcf") or value.endswith(".vcf.gz")) and os.path.isfile(value) and os.path.exists(value) and os.path.getsize(value)>0 :
                    self.file = value
                else:
                    print('[-f]或[--file]参数值 %s 不是一个有效的文件(扩展名为%s)' %(value,'.vcf或vcf.gz'))
                    sys.exit()
            elif arg=="-o" or arg == "--output":
                if value.endswith('.vcf') or value.endswith('.vcf.gz') :
                    self.output = value
                else:
                    print('[-o]或[--output]参数值 %s 不是一个有效的文件(扩展名为%s)' %(value,'.vcf或vcf.gz'))
                    sys.exit()
            elif arg=="-h" or arg == "--help":
                self.usage()
                sys.exit()
            elif arg=="--document":
                import VcfProcessUtil
                help(VcfProcessUtil)
                sys.exit()
            elif arg == "--version" or arg=="-v":
                print("版本: 1.00")
                sys.exit()
        self.cpath = os.getcwd()

    def process(self):
        """执行处理过程，处理传入的输入文件"""
        self.doProcess(self.file,self.output)

    def doProcess(self,filename,output):
        """
        对一个文件名为filename的文件执行过滤
        args:
            filename 一个文件名(str),将处理结果输出到为output文件
        return:
            无返回值,直接输出文件 -o 或者 --output
        """
        if filename==output:
            print("输入文件名和输出文件名不能相同")
            exit()

        fullpath   = os.path.realpath(filename)
        workingDir = os.path.dirname(fullpath)
        if workingDir is None:
            workingDir = "."
        os.chdir(workingDir)

        start = time.time()

        self.header = self.doProcessHeader(filename)

        if os.path.exists(filename): 
            if filename.endswith('.vcf.gz'):
                self.data = pandas.read_csv(
                    filename,
                    skiprows=self.skip, 
                    #converters={'#CHROM':str},
                    header=0,
                    sep='\t',
                    compression ='gzip',
                    #skip_blank_lines=True,
                    engine='c'
                )
            else:
                self.data = pandas.read_csv(
                    filename,
                    skiprows=self.skip, 
                    #converters={'#CHROM':str},
                    header=0,
                    sep='\t',
                    #skip_blank_lines=True,
                    engine='c'
                )

            col1 = self.data.iloc[:, 0]
            col1 = col1.apply(lambda x:self.toChr(x))
            self.data.iloc[:, 0] = col1
            #self.data['#CHROM'] = self.data.apply(lambda row:self.toChr(row['#CHROM']),axis=1)

            print(self.data)

            if output.endswith('.vcf.gz'):    
                #写入文件头到output
                with gzip.open(output, 'wt') as file:
                    file.write(self.header)
                #追加模式
                self.data.to_csv(
                    output,
                    sep='\t',
                    index=False,
                    encoding='utf-8',
                    header=0,
                    mode='a'
                )
            elif output.endswith('.vcf'):
                #写入文件头到output
                with open(output, 'wt') as file:
                    file.write(self.header)
                #追加模式
                self.data.to_csv(
                    output,
                    sep='\t',
                    index=False,
                    encoding='utf-8',
                    header=0,
                    mode='a'
                )
        
        stop  = time.time()
        self.log.info('RUN time:'+str(round((stop-start),3))+' seconds')
    
    def doProcessHeader(self,filename):
        '''
            记录vcf文件头,如果文件头包含contig信息则直接pass,避免复杂坐标转换,后续如果需要则改进
        '''
        header = ''
        if filename.endswith('.vcf') :
            with open(filename,'r') as f:
                for line in f:
                    li = line.decode('utf-8')
                    if li.startswith('##'):
                        self.skip=self.skip+1
                        if li.startswith('##contig='):
                            continue
                        header+=li
                    else:
                        break
        elif filename.endswith('.vcf.gz'):
            with gzip.open(filename,'r') as f:
                for line in f:
                    li = line.decode('utf-8')
                    if li.startswith('##'):
                        self.skip=self.skip+1
                        if li.startswith('##contig='):
                            continue
                        header+=li
                    else:
                        break
        return header

    def toChr(self,chrom):
        '''
            修改CHROM字段，序号前加chr
        '''
        ch = str(chrom)
        if not (ch.startswith('MT') or ch.startswith('NW')) :
            return 'chr'+ch
        elif ch.startswith('chr'):
            return ch
        return ch


    def usage(self):
        '''
        打印输出Usage
        '''
        print('Usage : ./VcfProcessUtil.py [OPTION]... 或者 python VcfProcessUtil [OPTION]')
        print('''
            用于.vcf格式CHROM列坐标转换,将数字格式转换为chr+数字格式
            Example: VcfProcessUtil.py  -f clinvar_20190708.vcf -o clinvar_20190708.trans.vcf
        ''')
        print('''部分使用方法及参数如下：\n''')
        print('-f, --file=\t处理文件扩展名为.vcf或.vcf.gz')
        print('-o, --output=\t输出文件扩展名为.vcf或vcf.gz')
        print('-h, --help\t显示帮助')
        print('-v, --version\t显示版本号')
        print('--document\t显示开发文档')
        print('\n')
        print('提交bug,to <6041738@qq.com>. 网站:https://sliverworkspace.com \n')


if __name__ == '__main__':
    f=VcfProcessUtil()
    f.process()