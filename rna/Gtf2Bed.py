#!/usr/bin/python3
#-*-coding:utf-8-*-

__author__ = '豆浆包子'
__all__    = ['Gtf2Bed', 'process', 'doProcess','usage']

import os,sys,re,pandas,numpy,getopt,logging

class Gtf2Bed(object):
  '''过滤vcf文件,输出Filter Annotation '''

  def __init__(self):
    '''
    初始化:验证工具,获取脚本参数,Usage使用说明等
    '''
    logging.basicConfig(
      level   = logging.INFO,
      format  = '%(asctime)s - %(filename)s[line:%(lineno)d] - %(levelname)s: %(message)s'
    )
    self.log        = logging.getLogger(__name__)
    self.file       = None
    self.out        = None

    try:
      self.opts = getopt.getopt(sys.argv[1:], "f:o:hv", ["file=","out=","version", "help","document"])
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
        if os.path.isfile(value) and os.path.getsize(value)>0:
          self.file = value
        else:
          print('[-f]或[--file]参数值 %s 不是一个有效的文件(以%s结尾)' %(value,'.vcf 或 .vcf.gz'))
          sys.exit()
      elif arg=="-o" or arg == "--out":
        if value.endswith('.bed'):
          self.out = value
          if os.path.isfile(value) :
            self.log.warn('[-o]或[--out]参数值 %s 文件已经存在，输出结果将覆盖该文件' %(value))
        else:
          print('[-o]或[--out]参数值 %s 不是一个有效的文件(以%s结尾)' %(value,'.bed'))
      elif arg=="-h" or arg == "--help":
        self.usage()
        sys.exit()
      elif arg=="--document":
        import Gtf2Bed
        help(Gtf2Bed)
        sys.exit()
      elif arg == "--version" or arg=="-v":
        print("版本: 1.00")
        sys.exit()
    self.cpath = os.getcwd()
    ps_stats = False
    if self.file is None:
        print('[-f]或[--file]参数值不能为空')
        ps_stats = True
    if self.out is None:
        print('[-o]或[--out]参数值不能为空')
        ps_stats = True
    if ps_stats:
        sys.exit()


  
  def process(self):
    """执行过滤，过滤传入的-f参数文件"""
    self.doProcess(self.file,self.out)

  def doProcess(self,filename,out):
    """
    对一个文件名为filename的文件执行过滤
    args:
      filename 一个文件名(str)
      out      输出文件路径和文件名
    return:
      无返回值
    """

    if (os.path.exists(filename)) and (os.path.getsize(filename)>0):
      if filename.endswith('.gz'):
        self.data       = pandas.read_csv(
          filename,
          compression='gzip',
          names=['seqname','source','feature','start','end','score','strand','frame','attribute'],
          sep='\t',
          header=None,
          skip_blank_lines=True,
          comment='#'
        )
      else :
        self.data       = pandas.read_csv(
          filename,
          names=['seqname','source','feature','start','end','score','strand','frame','attribute'],
          sep='\t',
          header=None,
          skip_blank_lines=True,
          comment='#'
        )
      
      if not self.data.empty:
        self.data          = self.data[['seqname','start','end']]
        self.data['start'] = self.data.apply(lambda row:self.fix(row['start']),axis=1)
        print(self.data)
        self.log.info('writting file to out:%s',out)
        self.data.to_csv(out,index=False,encoding='utf-8',sep='\t',header=None)      
      else:
        self.log.info('save empty file to: %s',out)
        emptyFile = open(out,'w')
        emptyFile.write('')
        emptyFile.close()
    else:
      self.log.info('save empty file to: %s',out)
      emptyFile = open(out,'w')
      emptyFile.write('')
      emptyFile.close()

  def fix(self,posi):
    return posi-1

  def usage(self):
    '''
    打印输出Usage
    '''
    print('Usage : ./Gtf2Bed [OPTION]... 或者 python Gtf2Bed [OPTION]')
    print('''
      GTF文件，输出格式为.bed
      Example: python Gtf2Bed.py -f gencode.v19.annotation.gtf.gz -o RNA-Seq.bed
    ''')
    print('''部分使用方法及参数如下：\n''')
    print('-f, --file=\t处理一个扩展名为gtf或gtf.gz的文件')
    print('-o, --out=\t处理结果输出文件,扩展名为.bed')
    print('-h, --help\t显示帮助')
    print('-v, --version\t显示版本号')
    print('--document\t显示开发文档')
    print('\n')
    print('提交bug,to <6041738@qq.com>.\n')


if __name__ == '__main__':
  f = Gtf2Bed()
  f.process()