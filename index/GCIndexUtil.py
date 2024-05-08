#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = '豆浆包子'
__all__    = [
    'GCIndexUtil', 
    'process', 
    'doProcess',
    'calcSkip',
    'usage'
]

import os,logging,time,sys,getopt,pandas

class GCIndexUtil(object):
    """ 
        从SampleSheet.csv文件中计算index ATGC占比
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
        self.skip   = -1

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
                if value.endswith(".csv") and os.path.isfile(value) and os.path.exists(value) and os.path.getsize(value)>0 :
                    self.file = value
                else:
                    print('[-f]或[--file]参数值 %s 不是一个有效的文件(扩展名为%s)' %(value,'.csv'))
                    sys.exit()
            elif arg=="-o" or arg == "--output":
                if value.endswith('.csv') :
                    self.output = value
                else:
                    print('[-o]或[--output]参数值 %s 不是一个有效的文件(扩展名为%s)' %(value,'.csv'))
                    sys.exit()
            elif arg=="-h" or arg == "--help":
                self.usage()
                sys.exit()
            elif arg=="--document":
                import index.GCIndexUtil as GCIndexUtil
                help(GCIndexUtil)
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

        self.skip = self.calcSkip(filename)
        if self.skip == -1:
            self.log.info('未找到有效数据行')
            return
        else:
            self.log.info('Find [Data] Line %d',self.skip)
            cols  = ['position','index1-A','index1-T','index1-G','index1-C','index1-GC','index2-A','index2-T','index2-G','index2-C','index2-GC']
            res   = pandas.DataFrame(columns=cols)
            print(res)
            if os.path.exists(filename): 
                self.data = pandas.read_csv(
                    filename,
                    names=[
                        'Sample_ID',
                        'Sample_Name',
                        'Sample_Plate',
                        'Sample_Well',
                        'I7_Index_ID',
                        'index1',
                        'I5_Index_ID',
                        'index2',
                        'Sample_Project',
                        'Description'
                    ],
                    skiprows=self.skip-1,
                    header=0,
                    sep=',',
                    engine='c'
                )
                print(self.data)
                
                index1_len_max = self.data['index1'].str.len().max()
                self.log.info('Index1 Max len : %d',index1_len_max)      
                index2_len_max = self.data['index2'].str.len().max()
                self.log.info('Index2 Max len : %d',index2_len_max)
                
                for i in range(max(index1_len_max,index2_len_max)):
                    tmp = [i+1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
                    A1  = 0
                    T1  = 0
                    G1  = 0
                    C1  = 0
                    GC1 = 0
                    MX1 = len(self.data.values)
                    A2  = 0
                    T2  = 0
                    G2  = 0
                    C2  = 0
                    GC2 = 0
                    MX2 = len(self.data.values)
                    for row in self.data.values:
                        index1_s = row[5].upper()
                        index2_s = row[7].upper()
                        if i<=len(index1_s):
                            base1 = index1_s[i]
                            if   base1=='A':
                                A1 +=1
                            elif base1=='T':
                                T1 +=1
                            elif base1=='G':
                                G1 +=1
                                GC1+=1
                            elif base1=='C':
                                C1 +=1
                                GC1+=1
                        if i<=len(index2_s):
                            base2 = index2_s[i]
                            if   base2=='A':
                                A2 +=1
                            elif base2=='T':
                                T2 +=1
                            elif base2=='G':
                                G2 +=1
                                GC2+=1
                            elif base2=='C':
                                C2 +=1
                                GC2+=1
                    
                    tmp[1]  = A1/MX1
                    tmp[2]  = T1/MX1
                    tmp[3]  = G1/MX1
                    tmp[4]  = C1/MX1
                    tmp[5]  = GC1/MX1
                    tmp[6]  = A2/MX2
                    tmp[7]  = T2/MX2
                    tmp[8]  = G2/MX2
                    tmp[9]  = C2/MX2
                    tmp[10] = GC2/MX2
                    res.loc[len(res.index)] = tmp
                res['position'] = res['position'].astype(int)
                for idx in range(1,len(cols)):
                    res[cols[idx]] = res[cols[idx]].apply('{:.2f}'.format)
                print(res)
                self.log.info('writting file to : %s',output)
                res.to_csv(output,index=False,encoding='utf-8',sep=',')
            stop  = time.time()
            self.log.info('RUN time : '+str(round((stop-start),3))+' seconds')

    def calcSkip(self,filename):
        if os.path.exists(filename):
            with open(filename,'r') as f:
                for index, line in enumerate(f, start=1):
                    if line.startswith('Sample_ID'):
                        return index
                return -1
        return -1

    def usage(self):
        '''
        打印输出Usage
        '''
        print('Usage : ./GCIndexUtil.py [OPTION]... 或者 python GCIndexUtil [OPTION]')
        print('''
            从SampleSheet.csv文件中计算index ATGC占比并输出计算结果为.csv文件
            Example: GCIndexUtil.py  -f samplesheet.csv -o samplesheetGc.csv
        ''')
        print('''部分使用方法及参数如下：\n''')
        print('-f, --file=\t输入文件扩展名为.csv')
        print('-o, --output=\t输出文件扩展名为.csv')
        print('-h, --help\t显示帮助')
        print('-v, --version\t显示版本号')
        print('--document\t显示开发文档')
        print('\n')
        print('提交bug,to <6041738@qq.com>. https://github.com/doujiangbaozi/sliverworkspace-util \n')


if __name__ == '__main__':
    f=GCIndexUtil()
    f.process()