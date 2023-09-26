#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
__author__ = '豆浆包子'
__all__    = ['MatchedSnvVepAnnotationFilter', 'process', 'doProcess','usage']

import os,sys,re,pandas,numpy,getopt,logging


class MatchedSnvVepAnnotationFilter(object):
    """ 输入vcf文件，和annovar注释文件,输出注释结果 """

    def __init__(self):
        '''
        初始化:验证工具,获取脚本参数,Usage使用说明等
        '''
        logging.basicConfig(
            level   = logging.INFO,
            format  = '%(asctime)s - %(filename)s[line:%(lineno)d] - %(levelname)s: %(message)s'
        )
        self.log = logging.getLogger(__name__)

        self.filters = [
            'PASS',
            'map_qual',
            'strand_bias',
            'clustered_events',
            'duplicate_evidence',
            'multiallelic',
            'germline',
            'normal_artifact',
            'slippage',
            'haplotype',
            'base_qual',
            'str_contraction'
        ]
        self.exclude = set()
        self.include = set()

        self.include.add('PASS')

        self.vcf     = None
        self.out     = None
        self.anno    = None
        self.vaf     = None
        self.tlod    = None
        self.depth   = None

        try:
            self.opts = getopt.getopt(
                sys.argv[1:], 
                "a:o:v:i:e:dh", 
                ["annotation=","out=","vcf=","include=","exclude=","min-vaf=","min-tlod=","min-depth=","document","help"]
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
            if arg=="-v" or arg == "--vcf":
                if (value.endswith('.vcf') or value.endswith('.vcf.gz')) and os.path.isfile(value) :
                    self.vcf = value
                else:
                    print('[-v]或[--vcf]参数值 %s 不是一个有效的文件(以%s结尾)' %(value,'.vcf 或 .vcf.gz'))
                    sys.exit()
            elif arg=="-a" or arg == "--annotation":
                if os.path.isfile(value) :
                    self.anno = value
                else:
                    print('[-a]或[--annotation]参数值 %s 不是一个有效的Annovar注释输出文件' %(value))
                    sys.exit()
            elif arg=="-o" or arg == "--out":
                if value.endswith('.tsv'):
                    self.out = value
                    if os.path.isfile(value) :
                        self.log.warn('[-o]或[--out]参数值 %s 文件已经存在，输出结果将覆盖该文件' %(value))
                else:
                    print('[-o]或[--out]参数值 %s 不是一个有效的文件(以%s结尾)' %(value,'.tsv'))
            elif arg=="-e" or arg == "--exclude":
                for filter in self.filters:
                    if value==filter:
                        self.exclude.add(value)
            elif arg=="-i" or arg == "--include":
                for filter in self.filters:
                    if value==filter:
                        self.include.add(value)
            elif arg=="--min-vaf":
                pattern = re.compile(r'^[-+]?[0-9]+\.[0-9]+$')
                result = pattern.match(value)
                if value.isdigit() or result:
                    self.vaf = float(value)
                    self.log.info('--min-vaf=%s',value)
                else:
                    print('[--min-vaf]参数值 %s 不是一个有效的数值 ' %(value))
            elif arg=="--min-tlod":
                pattern = re.compile(r'^[-+]?[0-9]+\.[0-9]+$')
                result = pattern.match(value)
                if value.isdigit() or result:
                    self.tlod = float(value)
                    self.log.info('--min-tlod=%s',value)
                else:
                    print('[--min-tlod]参数值 %s 不是一个有效的数值 ' %(value))
            elif arg=="--min-depth":
                if value.isdigit():
                    self.depth = int(value)
                    self.log.info('--min-depth=%s',value)
                else:
                    print('[--min-depth]参数值 %s 不是一个有效的正整数 ' %(value))
            elif arg=="-h" or arg == "--help":
                self.usage()
                sys.exit()
            elif arg=="-d" or arg=="--document":
                import MatchedSnvVepAnnotationFilter
                help(MatchedSnvVepAnnotationFilter)
                sys.exit()
        ps_status = False
        if self.vcf is None:
            print('[-v]或[--vcf]参数值不能为空')
            ps_status = True
        if self.out is None:
            print('[-o]或[--out]参数值不能为空')
            ps_status = True
        if self.anno is None:
            print('[-a]或[--annotation]参数值不能为空')
            ps_status = True
        if ps_status:
            sys.exit()
        self.cpath = os.getcwd()
        self.log.info('WorkingDir is : %s' % self.cpath)



    def process(self):
        """处理传入的vcf文件和annovar注释文件，输出结果"""
        self.doProcess(self.vcf,self.anno,self.out)

    def doProcess(self,vcf,anno,out):
        """
        处理传入的vcf文件和annovar注释文件，输出结果
        args:
            vcf    vcf文件
            anno   annovar注释结果文件
            out    输出分析结果文件
        return:
            无返回值
        """
        #if len(self.include)==0:
            #self.include.add('clustered_events')
            #self.include.add('strand_bias')
        #if len(self.exclude)==0:
            #self.exclude.add('normal_artifact')

        self.log.info("Exclude Filters : %s",self.exclude)
        self.log.info("Include Filters : %s",self.include)

        if (os.path.exists(vcf)) and (os.path.getsize(vcf)>0) and (os.path.exists(anno)) and (os.path.getsize(anno)>0):
            
            if vcf.endswith('.vcf.gz'):
                self.vcf_data = pandas.read_csv(
                    vcf,
                    names=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','NORMAL','SAMPLE'],
                    sep='\t',
                    header=None,
                    compression ='gzip',
                    skip_blank_lines=True,
                    comment='#')
            else:
                self.vcf_data = pandas.read_csv(
                    vcf,
                    names=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','NORMAL','SAMPLE'],
                    sep='\t',
                    header=None,
                    skip_blank_lines=True,
                    comment='#')
            
            self.vcf_data            = self.vcf_data.apply(lambda row:self.adjustVariationPosition(row),axis=1)

            self.vcf_data['START']   = self.vcf_data.POS
            self.vcf_data['DP']      = self.vcf_data.apply(lambda row:self.calculateDP(row),axis=1)
            self.vcf_data['DP']      = self.vcf_data['DP'].astype('int64')

            self.vcf_data['VAF']     = self.vcf_data.apply(lambda row:self.calculateVAF(row),axis=1)
            self.vcf_data['VAF']     = self.vcf_data.VAF.round(decimals=5)
            self.vcf_data['TLOD']    = self.vcf_data.apply(lambda row:self.calculateTLOD(row),axis=1)

            self.vcf_data            = self.vcf_data[['CHROM','START','REF','ALT','DP','VAF','TLOD','FILTER']]
            
            anno_data = pandas.read_csv(
                anno,
                header=0,
                names=[
                    'Uploaded_variation','Location','Allele','Gene','Feature','Feature_type','Consequence','cDNA_position','CDS_position','Protein_position','Amino_acids','Codons','Existing_variation','IMPACT','DISTANCE','STRAND','FLAGS','VARIANT_CLASS','SYMBOL','SYMBOL_SOURCE','HGNC_ID','BIOTYPE','CANONICAL','MANE_SELECT','MANE_PLUS_CLINICAL','TSL','APPRIS','CCDS','ENSP','SWISSPROT','TREMBL','UNIPARC','UNIPROT_ISOFORM','REFSEQ_MATCH','REFSEQ_OFFSET','GIVEN_REF','USED_REF','BAM_EDIT','GENE_PHENO','SIFT','PolyPhen','EXON','INTRON','DOMAINS','miRNA','HGVSc','HGVSp','HGVS_OFFSET','AF','AFR_AF','AMR_AF','EAS_AF','EUR_AF','SAS_AF','AA_AF','EA_AF','gnomAD_AF','gnomAD_AFR_AF','gnomAD_AMR_AF','gnomAD_ASJ_AF','gnomAD_EAS_AF','gnomAD_FIN_AF','gnomAD_NFE_AF','gnomAD_OTH_AF','gnomAD_SAS_AF','MAX_AF','MAX_AF_POPS','CLIN_SIG','SOMATIC','PHENO','PUBMED','MOTIF_NAME','MOTIF_POS','HIGH_INF_POS','MOTIF_SCORE_CHANGE','TRANSCRIPTION_FACTORS'
                ],
                sep='\t',
                skip_blank_lines=True,
                comment='#'
            )
            anno_data = anno_data.apply(lambda row:self.processAnnotation(row),     axis=1)
            anno_data = anno_data.apply(lambda row:self.processAdjustVariation(row),axis=1)
            
            anno_data['POS1'] = anno_data['POS1'].astype('int64')
            
            self.log.info("Annotation prepared for calculate: %s",anno)
            
            self.vcf_data = self.vcf_data.apply(lambda row:self.processAnnotationDepAndVaf(row,anno_data),axis=1) 
            
            if not self.vcf_data.empty:

                #删除NaN值的行
                self.vcf_data = self.vcf_data.dropna(axis=0, how='any')

                if len(self.exclude)>0 and (not self.vcf_data.empty):
                    self.vcf_data['exclude']=self.vcf_data.apply(lambda row:False,axis=1)
                    self.vcf_data = self.vcf_data.apply(lambda row:self.processExclude(row),axis=1)
                    self.vcf_data = self.vcf_data[self.vcf_data['exclude']==False]
                    self.vcf_data.drop(['exclude'],axis=1,inplace=True)
                
                if len(self.include)>0 and (not self.vcf_data.empty):
                    self.vcf_data['include']=self.vcf_data.apply(lambda row:False,axis=1)
                    self.vcf_data = self.vcf_data.apply(lambda row:self.processInclude(row),axis=1)
                    self.vcf_data = self.vcf_data[self.vcf_data['include']==True]
                    self.vcf_data.drop(['include'],axis=1,inplace=True)

                self.vcf_data = self.vcf_data.reset_index(drop=True)
                
                if self.vaf is not None and (not self.vcf_data.empty):
                    self.vcf_data        = self.vcf_data[self.vcf_data['VAF']>=self.vaf]
                self.vcf_data['VAF']     = self.vcf_data.VAF.apply(lambda x: format(x, '.2%'))
                
                if self.tlod is not None and (not self.vcf_data.empty):
                    self.vcf_data['TLOD_FILTERED'] = self.vcf_data.apply(lambda row:self.filterTLOD(row['TLOD']),axis=1)
                    self.vcf_data                  = self.vcf_data[self.vcf_data['TLOD_FILTERED']==True]
                    self.vcf_data.drop(['TLOD_FILTERED'],axis=1,inplace=True)
                
                
                if self.depth is not None and (not self.vcf_data.empty):
                    self.vcf_data        = self.vcf_data[self.vcf_data['DP']>=self.depth]
                
                self.log.info('writting result to file %s',out)
                
                self.vcf_data = self.vcf_data[[
                    'Chr','Start','End','Ref','Alt','Gene','Type','Transcript', 'cHGVS', 'pHGVS','VAF','DP','FILTER','Existing_variation','BIOTYPE','CLIN_SIG','PUBMED'
                ]]
                self.vcf_data['Start'] = self.vcf_data['Start'].astype('int64')
                print(self.vcf_data)
                self.vcf_data.to_csv(out,index=False,encoding='utf-8',sep='\t')
                
            else:
                self.log.info('writting empty result to file %s',out)
                emptyFile = open(out,'w')
                emptyFile.write('')
                emptyFile.close()
            
        else:
            self.log.info('writting empty result to file %s',out)
            emptyFile = open(out,'w')
            emptyFile.write('')
            emptyFile.close()
    
    def filterTLOD(self,str_tlod):
        lods = str_tlod.split(',')
        if self.tlod is None:
            return False
        if len(lods)<=0:
            return False
        try:
            for tmp in lods:
                pattern = re.compile(r'^[-+]?[0-9]+\.[0-9]+$')
                result = pattern.match(tmp)
                if not (tmp.isdigit() or result):
                    return False
                elif float(tmp)<self.tlod:
                    return False
            return True
        except Exception as e:
            print(e)
            return False

    def calculateVAF(self,row):
        '''根据row/行数据计算该行的突变频率/丰度'''
        try:
            keys         = row['FORMAT'].split(':')
            vals         = row['SAMPLE'].split(':')
            for i in range(min(len(keys),len(vals))):
                if keys[i]=='AD':
                    vaf_str = vals[i]
                    if len(vaf_str)>0 and vaf_str.find(',')!=-1:
                        arry  = vaf_str.split(',')
                        ref_n = float(arry[0])
                        alt_n = float(arry[1])
                        if (ref_n+alt_n)>0:
                            return (alt_n/(ref_n+alt_n))
                        else:
                            return 0.00
        except Exception as e:
            print(e)
            return 0.00
    
    def calculateTLOD(self,row):
        '''根据row/行数据，获取INFO字段中的TLOD值，过滤cutoff值'''
        try:
            info_fields = row['INFO'].split(';')
            for field in info_fields:        
                if field.find('TLOD=') != -1:
                    tlod = field[field.find('TLOD=')+len('TLOD='):]
                    return tlod
        except Exception as e:
            print(e)
            return 'NAN'

    def calculateDP(self,row):
        '''根据row/行数据，获取INFO字段中的DP值，即测序深度值'''
        try:
            info_fields = row['INFO'].split(';')
            for temp in info_fields:
                if temp.find('DP=') != -1:
                    return int(temp[temp.find('DP=')+len('DP='):])
            return 0
        except Exception as e:
            print(e)
            return 0

    def calculateChrNum(self,chr):
        '''根据染色体编号chrx计算后面纯数字字母x'''
        try:
            field = chr.lower()
            if field.find('chr') != -1:
                return field[field.find('chr')+len('chr'):].upper()
        except Exception as e:
            print(e)
            return 'NAN'


    def findSamePrefix(self,ref,arr):
        '''
            查找ref和arr[str]数组中相同字符前缀
        '''
        tarr = arr+[ref]
        result = ''
        for i in range(len(ref)):
            test   = result+ref[i]
            status = True
            for cursor in range(len(tarr)):
                if tarr[cursor].startswith(test):
                    continue
                else:
                    status = False
                    return result
            result = test
        return result

    def findSameSuffix(self,ref,arr):
        '''
            查找ref和arr[str]数组中相同字符后缀
        '''
        tarr = arr+[ref]
        result = ''
        for i in range(len(ref)):
            test   = result+ref[i]
            status = True
            for cursor in range(len(tarr)):
                if tarr[cursor].endswith(test):
                    continue
                else:
                    status = False
                    return result
            result = test
        return result

    def processExclude(self,row):
        '''
            根据FILTER info字段过滤指定的过滤器，匹配返回True，不匹配返回False
        '''
        try:
            for filter in iter(self.exclude):
                if row['exclude']==True:
                    return row
                fields = row['FILTER'].split(';')
                for i in range(len(fields)):
                    if filter==fields[i]:
                        row['exclude']=True
                        return row
                continue
            return row
        except Exception as e:
            print(e)
            return row

    def processInclude(self,row):
        '''
            根据FILTER info字段过滤指定的过滤器，匹配返回True，不匹配返回False
        '''
        try:
            for filter in iter(self.include):
                if row['include']==True:
                    return row
                fields = row['FILTER'].split(';')
                for i in range(len(fields)):
                    if filter==fields[i]:
                        row['include']=True
                        return row
                continue
            return row
        except Exception as e:
            print(e)
            return row

    
    def adjustContinuityPosition(self,row):
        '''
            连续等效突变合并
        '''
        return

    def adjustVariationPosition(self,row):
        '''
            调整突变坐标使vep的突变坐标和gatk输出突变一致
        '''
        try:
            ref = row['REF']
            arr = row['ALT'].split(',')
            sta = self.findSamePrefix(ref,arr)
            if sta!='':
                row['POS']=row['POS']+len(sta)
                if ref==sta:
                    row['REF']='-'
                else:
                    row['REF']=row['REF'][row['REF'].find(sta)+len(sta):]
                for i in range(len(arr)):
                    if sta==arr[i]:
                        arr[i]='-'
                    else:
                        arr[i]=arr[i][arr[i].find(sta)+len(sta)]
                row['ALT']=','.join(arr)
            
            sta = self.findSameSuffix(ref,arr)
            if sta!='':
                row['POS']=row['POS']-len(sta)
                if ref==sta:
                    row['REF']='-'
                else:
                    row['REF']=row['REF'][:row['REF'].rfind(sta)-1]
                for i in range(len(arr)):
                    if sta==arr[i]:
                        arr[i]='-'
                    else:
                        arr[i]=arr[i][:arr[i].find(sta)-1]
                row['ALT']=','.join(arr)
            return row
        except Exception as e:
            print(e)
            return row

    def processAnnotation(self,row):
        '''
            vep注释文件中第一列Uploaded_variation拆分成和POSI与vcf文件匹配
        '''
        try:
            fields = row['Uploaded_variation'].split('_')
            if len(fields)==3:
                row['CHROM1']= fields[0]
                row['POS1']  = fields[1]
                vs           = fields[2].split('/')
                row['REF']   = row['USED_REF']
                row['ALT']   = row['Allele']
            source           = row['VARIANT_CLASS']
            if source=='SNV':
                row['type']  = source
            elif source=='deletion':
                row['type']  = 'DEL'
            elif source=='insertion':
                row['type']  = 'INS'
            else:
                row['type']  = 'Complex'

            row['Transcript']= row['Feature'][:row['Feature'].find('.')]
            
            cHGVS_fields     = row['HGVSc'].split(':')
            if len(cHGVS_fields)>1:
                row['cHGVS'] = cHGVS_fields[1]
            else:
                row['cHGVS'] = '-'
            
            HGVSp_fields     = row['HGVSp'].split(':')
            if len(HGVSp_fields)>1:
                row['pHGVS'] = HGVSp_fields[1]
            else:
                row['pHGVS'] = '-'
            return row
        except Exception as e:
            print(e)
            return None

    def processAdjustVariation(self,row):
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
            return row
        except Exception as e:
            print(e)
            print(row)
            return row

    def processCalculateEndPosition(self,row):
        '''
            根据row Start,Ref,Alt Type计算 End坐标
        '''
        try:
            start = row['Start']
            ref   = row['Ref']
            tp    = row['Type']
            if   tp=='DEL':
                return start+len(ref)-1
            elif tp=='Complex':
                return start+len(ref)-1
            elif tp=='SNV':
                return start+len(ref)-1
            elif tp=='INS':
                return '-'
        except Exception as e:
            print(e)

    def processAnnotationDepAndVaf(self,row,anno):
        '''
            结合vcf文件和vep注释结果，匹配为vep注释文件添加DP和VAF信息，并导出最终结果
        '''
        try:
            for i in range(len(anno)):
                r = anno.iloc[i]
                if (row['CHROM']==r['CHROM1']) and (row['START']==r['POS1']):
                    row['Chr']        = self.calculateChrNum(row['CHROM'])
                    row['Start']      = int(r['Start'])
                    row['Ref']        = r['REF']
                    row['Alt']        = r['ALT']
                    row['Gene']       = r['SYMBOL']
                    row['Type']       = r['type']
                    row['Transcript'] = r['Transcript']
                    row['cHGVS']      = r['cHGVS']
                    row['pHGVS']      = r['pHGVS']
                    '''
                        增加几列数据便于判断
                    '''
                    row['Existing_variation'] = r['Existing_variation']
                    row['BIOTYPE']            = r['BIOTYPE']
                    row['CLIN_SIG']           = r['CLIN_SIG']
                    row['PUBMED']             = r['PUBMED']
                    
                    #todo 此处待测试
                    row['End']        = self.processCalculateEndPosition(row)#r['End']
                    #print(row)
                    break
            #print(row)
            return row
        except Exception as e:
            print(e)
            return row


    def usage(self):
        '''
            打印输出Usage
        '''
        print('Usage : ./MatchedSnvVepAnnotationFilter [OPTION]... 或者 python MatchedSnvVepAnnotationFilter [OPTION]')
        print('''
            根据vcf文件，annovar注释输出结果，输出最终分析结果文件，格式为.tsv
            Example:\t
                ./MatchedSnvVepAnnotationFilter.py  \\
                    -e normal_artifact \\
                    -e germline   \\
                    -i strand_bias   \\
                    -i clustered_events   \\
                    --min-vaf=0.01   \\
                    --min-tlod=16.00 \\
                    --min-depth=500  \\
                    -v /opt/result/B1701_filtered_snpeff.vcf \\
                    -a /opt/result/202101/snv/202101_filtered_vep.tsv \\
                    -o /opt/result/B1701.result.tsv
            或者\t
                python MatchedSnvVepAnnotationFilter.py  \\
                    --exclude=normal_artifact \\
                    --exclude=germline   \\
                    --include=strand_bias   \\
                    --include=clustered_events   \\
                    --min-vaf=0.01   \\
                    --min-tlod=16.00 \\
                    --min-depth=500  \\
                    --vcf=/opt/result/B1701_filtered_snpeff.vcf \\
                    --annotation=/opt/result/202101/snv/202101_filtered_vep.tsv \\
                    --out=/opt/result/B1701.result.tsv
        ''')
        print('''部分使用方法及参数如下：\n''')
        print('-v, --vcf=\t\t输入vcf或vcf.gz格式文件')
        print('-a, --annotation=\t输入annotation注释后的文件')
        print('-o, --out=\t\t输出处理过的结果文件')
        print('-e, --exclude=\t\t排除Gatk过滤后某过滤器的结果，可以使用多个（可选）')
        print('-i, --include=\t\t包含Gatk过滤后某过滤器的结果，可以使用多个（可选）')
        print('    --min-vaf=\t\t最小突变频率，小数如：0.01（可选）')
        print('    --min-tlod=\t\t最小TLOD值，小数如：16.00（可选）')
        print('    --min-depth=\t最小测序深度，正整数如：1000（可选）')
        print('-h, --help\t\t显示帮助')
        print('-d, --document\t\t显示开发文档')
        print('\n')
        print('提交bug,to <6041738@qq.com>.\n')

        '''
        python MatchedSnvVepAnnotationFilter.py  \
            --exclude=normal_artifact \
            --exclude=germline   \
            --include=strand_bias   \
            --include=clustered_events   \
            --min-vaf=0.01   \
            --min-tlod=16.00 \
            --min-depth=500  \
            --vcf=/opt/result/202101/snv/202101_filtered.vcf.gz \
            --annotation=/opt/result/202101/snv/202101_filtered_vep.tsv \
            --out=/opt/result/202101/snv/202101_anno.result.tsv
        
        MatchedSnvVepAnnotationFilter.py  \
            --exclude=normal_artifact \
            --exclude=germline   \
            --include=strand_bias   \
            --include=clustered_events   \
            --min-vaf=0.01   \
            --min-tlod=16.00 \
            --min-depth=500  \
            --vcf=/opt/result/202101/snv/202101_filtered.vcf.gz \
            --annotation=/opt/result/202101/snv/202101_filtered_vep.tsv \
            --out=/opt/result/202101/snv/202101_anno.result.tsv
        '''


if __name__ == '__main__':
    f = MatchedSnvVepAnnotationFilter()
    f.process()