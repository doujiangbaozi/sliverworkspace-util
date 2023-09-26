#!/usr/bin/python
# -*- coding: utf-8 -*-

__author__ = '豆浆包子'
__all__    = ['SvAnnotationFilter', 'filter', 'doFilter','usage']

import os,sys,re,pandas,numpy,getopt,logging

class SvAnnotationFilter(object):
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
    self.refGene    = None
    self.file       = None
    self.ref        = None
    self.cutoff     = 200

    try:
      self.opts = getopt.getopt(sys.argv[1:], "f:o:r:s:hv", ["file=","out=","ref=","score=","version", "help","document"])
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
        if (value.endswith('.vcf.gz') or value.endswith('.vcf')):
          self.file = value
        else:
          print('[-f]或[--file]参数值 %s 不是一个有效的文件(以%s结尾)' %(value,'.vcf 或 .vcf.gz'))
          sys.exit()
      elif arg=="-o" or arg == "--out":
        if value.endswith('.xls'):
          self.out = value
          if os.path.isfile(value) :
            self.log.warn('[-o]或[--out]参数值 %s 文件已经存在，输出结果将覆盖该文件' %(value))
        else:
          print('[-o]或[--out]参数值 %s 不是一个有效的文件(以%s结尾)' %(value,'.xls'))
      elif arg=="-r" or arg == "--ref":
        if os.path.isfile(value):
          self.refGene = value
        else:
          print('[-r]或[--ref]参数值 %s 不是一个有效的文件' %(value))
      elif arg=="-s" or arg == "--score":
        pattern = re.compile(r'^[-+]?[0-9]+\.[0-9]+$')
        result = pattern.match(value)
        if value.isdigit() or result:
          self.cutoff = float(value)
        else:
          print('[-s]或[--score]参数值 %s 不是一个有效的数值 ' %(value))
      elif arg=="-h" or arg == "--help":
        self.usage()
        sys.exit()
      elif arg=="--document":
        import SvAnnotationFilter
        help(SvAnnotationFilter)
        sys.exit()
      elif arg == "--version" or arg=="-v":
        print("版本: 1.00")
        sys.exit()
    self.cpath = os.getcwd()
    ps_stats = False
    if self.refGene is None:
        print('[-r]或[--ref]参数值不能为空')
        ps_stats = True
    if ps_stats:
        sys.exit()

  
  def calculateSomaticScore(self,info):
    info_fields = info.split(';')
    for field in info_fields:
      if field.find('SOMATICSCORE=')!=-1:
        return int(field[field.find('=')+1:])

  def calculateEnd(self,info):
    info_fields = info.split(';')
    for field in info_fields:
      if field.find('END=')!=-1:
        return field[field.find('=')+1:]
      else:
        return -1

  def calculateType(self,info):
    info_fields = info.split(';')
    for field in info_fields:
      if field.find('SVTYPE=')!=-1:
        return field[field.find('=')+1:]

  def calculateEvent(self,info):
    info_fields = info.split(';')
    for field in info_fields:
      if field.find('EVENT=')!=-1:
        return field[field.find('=')+1:]


  def compareChrom(self,c1,c2):
    v1 = c1[c1.find('chr')+len('chr'):]
    v2 = c2[c2.find('chr')+len('chr'):]
    if v1=='X' or v1=='Y':
      v1=23
    if v2=='X' or v2=='Y':
      v2=23
    if int(v1)>int(v2):
      return 1
    elif int(v1)==int(v2):
      return 0
    return -1
    

  def calculateChrom1(self,row):
    ALT = row['ALT']
    if ALT.find('[')!=-1 and ALT.find('[',ALT.find('[')+1)!=-1:
      temp = ALT[ALT.find('[')+1:ALT.find('[',ALT.find('[')+1)]
      if len(temp)>0:
        alt_chrom = temp.split(':')[0]
        if self.compareChrom(row['CHROM'],alt_chrom)<=0:
          return row['CHROM']
        else:
          return alt_chrom
    if ALT.find(']')!=-1 and ALT.find(']',ALT.find(']')+1)!=-1:
      temp = ALT[ALT.find(']')+1:ALT.find(']',ALT.find(']')+1)]
      if len(temp)>0:
        alt_chrom = temp.split(':')[0]
        if self.compareChrom(row['CHROM'],alt_chrom)<=0:
          return row['CHROM']
        else:
          return alt_chrom
    if ALT.find('[')==-1 and ALT.find(']')==-1:
      return row['CHROM']

  def calculateChrom2(self,row):
    ALT = row['ALT']
    if ALT.find('[')!=-1 and ALT.find('[',ALT.find('[')+1)!=-1:
      temp = ALT[ALT.find('[')+1:ALT.find('[',ALT.find('[')+1)]
      if len(temp)>0:
        alt_chrom = temp.split(':')[0]
        if self.compareChrom(row['CHROM'],alt_chrom)>=0:
          return row['CHROM']
        else:
          return alt_chrom
    if ALT.find(']')!=-1 and ALT.find(']',ALT.find(']')+1)!=-1:
      temp = ALT[ALT.find(']')+1:ALT.find(']',ALT.find(']')+1)]
      if len(temp)>0:
        alt_chrom = temp.split(':')[0]
        if self.compareChrom(row['CHROM'],alt_chrom)>=0:
          return row['CHROM']
        else:
          return alt_chrom
    if ALT.find('[')==-1 and ALT.find(']')==-1:
      return row['CHROM']

  def calculateBreakpoint1(self,row):
    if row['CHROM1'] == row['CHROM']:
      return row['POS']
    ALT=row['ALT']
    if ALT.find('[')!=-1 and ALT.find('[',ALT.find('[')+1)!=-1:
      temp = ALT[ALT.find('[')+1:ALT.find('[',ALT.find('[')+1)]
      if len(temp)>0:
        return temp.split(':')[1]
    if ALT.find(']')!=-1 and ALT.find(']',ALT.find(']')+1)!=-1:
      temp = ALT[ALT.find(']')+1:ALT.find(']',ALT.find(']')+1)]
      if len(temp)>0:
        return temp.split(':')[1]
      
  def calculateBreakpoint2(self,row):
    if row['CHROM1'] != row['CHROM']:
      return row['POS']
    ALT=row['ALT']
    if ALT.find('[')!=-1 and ALT.find('[',ALT.find('[')+1)!=-1:
      temp = ALT[ALT.find('[')+1:ALT.find('[',ALT.find('[')+1)]
      if len(temp)>0:
        return temp.split(':')[1]
    if ALT.find(']')!=-1 and ALT.find(']',ALT.find(']')+1)!=-1:
      temp = ALT[ALT.find(']')+1:ALT.find(']',ALT.find(']')+1)]
      if len(temp)>0:
        return temp.split(':')[1]
    if ALT.find('[')==-1 and ALT.find(']')==-1:
      return row['END']

  def calculateVAF(self,row):
    keys = row['FORMAT']
    if keys.find('SR')!=-1 and keys.find('PR')!=-1:
      arrays = row['SAMPLE'].split(':')
      spaned = arrays[0].split(',')
      split  = arrays[1].split(',')
      result = (float(split[1])+float(spaned[1]))/(float(split[0])+float(spaned[0])+float(split[1])+float(spaned[1]))
      return result
    if keys.find('SR')!=-1 or keys.find('PR')!=-1:
      arrays = row['SAMPLE'].split(':')
      temp   = arrays[0].split(',')
      result = float(temp[1])/(float(temp[0])+float(temp[1]))
      return result
  
  def calculateGene(self,row,ch,position):
    if self.ref is not None:
      chrom = row[ch]
      point = row[position]
      result  = self.ref[(self.ref['chrom']==chrom)&(point>=self.ref['txStart'])&(point<=self.ref['txEnd'])]
      if len(result)>0:
        return result['name2'].values[0]

  def calculateAnnotation(self,row):
    if row['GENE1'] is not None and row['GENE2'] is not None:
      return row['GENE1']+'-'+row['GENE2']
  

  def calculateChrNum(self,chr):
      '''根据染色体编号chrx计算后面纯数字字母x'''
      try:
          field = chr.lower()
          if field.find('chr') != -1:
              return field[field.find('chr')+len('chr'):].upper()
      except Exception,e:
          print(e)
          return 'NAN'
  
  def filter(self):
    """执行过滤，过滤传入的-f参数文件"""
    self.doFilter(self.file,self.out)

  def doFilter(self,filename,out):
    """
    对一个文件名为filename的文件执行过滤
    args:
      filename 一个文件名(str)
      out      输出文件路径和文件名
    return:
      无返回值
    """

    if (os.path.exists(filename)) and (os.path.getsize(filename)>0):
      
      if self.ref is None:
        '''初始化数据refGene.txt'''
        self.ref = pandas.read_csv(
          self.refGene,
          names=['bin','name','chrom','strand','txStart','txEnd','cdsStart','cdsEnd','exonCount','exonStarts','exonEnds',
            'score','name2','cdsStartStat','cdsEndStat','exonFrame'],
          sep='\t',
          header=None,
          skip_blank_lines=True,
          comment='#')
        self.ref.drop(['bin','strand','cdsStart','cdsEnd','exonCount','exonStarts','exonEnds',
            'score','cdsStartStat','cdsEndStat','exonFrame'],axis=1,inplace=True)
        self.ref['range'] = self.ref['txEnd']-self.ref['txStart']
        '''
        self.ref          = self.ref.groupby(['name2'],as_index=False).max()
        '''
        self.ref.sort_values("name2", inplace=True)
      
      #'somaticSV.vcf.gz'
      if out.endswith('.vcf.gz'):
        self.data       = pandas.read_csv(
          filename,
          compression='gzip',
          names=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','NORMAL','SAMPLE'],
          sep='\t',
          header=None,
          skip_blank_lines=True,
          comment='#'
        )
      else :
        self.data       = pandas.read_csv(
          filename,
          names=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','NORMAL','SAMPLE'],
          sep='\t',
          header=None,
          skip_blank_lines=True,
          comment='#'
        )
      
      if not self.data.empty:
        self.data['SOMATICSCORE'] = self.data.apply(lambda row:self.calculateSomaticScore(row['INFO']),axis=1)
        self.data['END'         ] = self.data.apply(lambda row:self.calculateEnd(row['INFO']),axis=1)
        self.data['TYPE'        ] = self.data.apply(lambda row:self.calculateType(row['INFO']),axis=1)
        self.data['EVENT'       ] = self.data.apply(lambda row:self.calculateEvent(row['INFO']),axis=1)

        self.data['CHROM1'      ] = self.data.apply(lambda row:self.calculateChrom1(row),axis=1)
        self.data['BREAKPOINT1' ] = self.data.apply(lambda row:self.calculateBreakpoint1(row),axis=1)
        self.data['CHROM2'      ] = self.data.apply(lambda row:self.calculateChrom2(row),axis=1)
        self.data['BREAKPOINT2' ] = self.data.apply(lambda row:self.calculateBreakpoint2(row),axis=1)

        self.data['BREAKPOINT1']  = self.data['BREAKPOINT1'].astype('int64')
        self.data['BREAKPOINT2']  = self.data['BREAKPOINT2'].astype('int64')

        self.data = self.data[self.data['SOMATICSCORE']>self.cutoff]
        
        if not self.data.empty:
          self.data.sort_values("CHROM", inplace=True)
          self.data['VAF' ]  = self.data.apply(lambda row:self.calculateVAF(row),axis=1)
          self.data['VAF']   = self.data.VAF.round(decimals=3)
          self.data['VAF']   = self.data.VAF.apply(lambda x: format(x, '.1%'))
          self.data.drop(['EVENT','CHROM','QUAL','FILTER','END','ID','POS','REF','ALT','INFO','FORMAT','SAMPLE','NORMAL'],axis=1,inplace=True)
          self.data['GENE1'] = self.data.apply(lambda row:self.calculateGene(row,'CHROM1','BREAKPOINT1'),axis=1)
          self.data['GENE2'] = self.data.apply(lambda row:self.calculateGene(row,'CHROM2','BREAKPOINT2'),axis=1)
          self.data['ANNOTATION'] = self.data.apply(lambda row:self.calculateAnnotation(row),axis=1)
          
          self.data['CHROM1']= self.data.apply(lambda row:self.calculateChrNum(row['CHROM1']),axis=1)
          self.data['CHROM2']= self.data.apply(lambda row:self.calculateChrNum(row['CHROM2']),axis=1)
          self.data = self.data.reset_index(drop=True)
      if not self.data.empty:
        self.data = self.data.drop_duplicates(subset=['CHROM1','CHROM2','BREAKPOINT1','BREAKPOINT2'], keep='first')
        print(self.data)
        self.log.info('writting file to out:%s',out)
        self.data.to_csv(out,index=False,encoding='utf-8',sep='\t')
      
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

      
  def usage(self):
    '''
    打印输出Usage
    '''
    print('Usage : ./SvAnnotationFilter [OPTION]... 或者 python SvAnnotationFilter [OPTION]')
    print('''
      用于处理变异注释文件.vcf.gz格式，过滤其中有效变异，与样本编号相符条件，输出格式为.xls
      Example: ./SvAnnotationFilter.py -r /opt/ref/hg19_refGene.txt -s 200 -f B1701/results/variants/somaticSV.vcf.gz -o result/1701.result.SV.xls
    ''')
    print('''部分使用方法及参数如下：\n''')
    print('-f, --file=\t处理一个扩展名为vcf或vcf.gz的文件')
    print('-o, --out=\t处理结果输出文件,扩展名为.xls')
    print('-r, --ref=\t引用基因组文件 hg19_refGene.txt')
    print('-s, --score=\t过滤参数，score数值')
    print('-h, --help\t显示帮助')
    print('-v, --version\t显示版本号')
    print('--document\t显示开发文档')
    print('\n')
    print('提交bug,to <6041738@qq.com>.\n')


if __name__ == '__main__':
  f = SvAnnotationFilter()
  f.filter()