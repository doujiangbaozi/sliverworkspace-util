#!/usr/bin/env python
#-*-coding:utf-8-*-

__author__ = '豆浆包子'
__all__    = [
    'FastqGenerator', 
    'process', 
    'genTit',
    'genQve',
    'genPos',
    'genSeq',
    'revSeq',
    'usage'
]

import os,sys,pandas,numpy,getopt,logging,time,subprocess,random,gzip
reload(sys)
sys.setdefaultencoding("utf-8")

class FastqGenerator(object):
    """ 
        根据reference文件，测序深度depth,读长length，bed文件生成Tumor-Normal中的Normal文件
    """

    def __init__(self):
        '''
            初始化:验证工具,获取脚本参数,Usage使用说明等
        '''
        logging.basicConfig(
            level   = logging.INFO,
            format  = '%(asctime)s - %(filename)s[line:%(lineno)d] - %(levelname)s: %(message)s'
        )
        self.log      = logging.getLogger(__name__)
        self.ref      = None #参考序列文件reference
        self.bed      = None #bed文件和panel相关，染色体范围
        self.depth    = None #测序深度，测序区域覆次数
        self.length   = None #reads读长 eg: 100,150,151
        self.out      = None #输出文件路径,输出 {name}_R1.fastq.gz {name}_R2.fastq.gz
        self.name     = None #默认样本名称为Normal

        #一下为调用工具路径
        self.samtools = 'samtools' #path下软件名称；samtools faidx /opt/ref/hg38/hg38.fa chr1:1-1000获取序列
        #self.bgzip    = None #bgzip文件路径，用于压缩输出文件成gz格式
        
        
        try:
            self.opts = getopt.getopt(
                sys.argv[1:], 
                "r:b:d:l:o:n:hv", 
                ["ref=", "bed=", "depth=", "len=", "out=", "name=", "help", "version", "document"]
            )
        except getopt.GetoptError:
            print('错误：请检查您的输入参数是否正确\n\n')
            self.usage()
            sys.exit()
        
        #如果没有输入参数，输出usage()
        if len(self.opts[0])==0:
            self.usage()
            sys.exit()
        
        #获取命令行参数值
        for arg, value in self.opts[0]:
            if arg=="-r" or arg == "--ref":
                if (value.endswith(".fa") or value.endswith('.fasta')) and os.path.isfile(value) and os.path.exists(value) and os.path.getsize(value)>0 :
                    self.ref = value
                else:
                    print('[-r]或[--ref]参数值 %s 不是一个有效的文件(以%s结尾)' %(value,'.fa 或 .fasta'))
                    sys.exit()
            elif arg=="-b" or arg == "--bed":
                if value.endswith(".bed") and os.path.exists(value) and os.path.isfile(value) and os.path.getsize(value)>0 :
                    self.bed = value
                else:
                    print('[-b]或[--bed]参数值 %s 不是一个有效的文件(以%s结尾)' %(value,'.bed'))
                    sys.exit()
            elif arg=="-d" or arg=="--depth":
                if value.isdigit():
                    self.depth = int(value)
                    #self.log.info('--depth=%s',value)
                else:
                    print('[-d] 或者 [--depth]参数值 %s 不是一个有效的正整数 ' %(value))
                    sys.exit()
            elif arg=="-l" or arg=="--len":
                if value.isdigit():
                    self.length = int(value)
                    #self.log.info('--length=%s',value)
                else:
                    print('[-l] 或者 [--len]参数值 %s 不是一个有效的正整数 ' %(value))
                    sys.exit()
            elif arg=="-o" or arg == "--out":
                if os.path.isdir(value) and os.path.exists(value) :
                    self.out = value
                elif not os.path.exists(value):
                    print('[-o]或[--out]参数值 %s 目录不存在' %(value))
                    sys.exit()
                elif not os.path.isdir(value):
                    print('[-o]或[--out]参数值 %s 不是目录'   %(value))
                    sys.exit()
            elif arg=="-n" or arg == "--name":
                if len(value)>0 :
                    self.name = value
                else:
                    print('[-n]或[--name]参数值 %s 不是一个有效的值' %(value))
                    sys.exit()
            elif arg=="-h" or arg == "--help":
                self.usage()
                sys.exit()
            elif arg=="--document":
                import FastqGenerator
                help(FastqGenerator)
                sys.exit()
            elif arg == "--version" or arg=="-v":
                print("\n  版本: 1.00\n")
                exit()
            
        status = False
        if self.ref is None:
            print('[-r]或[--ref]\t参数值不能为空')
            status = True
        if self.bed is None:
            print('[-b]或[--bed]\t参数值不能为空')
            status = True
        if self.depth is None:
            print('[-d]或[--depth]\t参数值不能为空')
            status = True
        if self.length is None:
            print('[-l]或[--len]\t参数值不能为空')
            status = True
        if self.name is None:
            print('[-n]或[--name]\t参数值不能为空')
            status = True
        if self.out is None:
            print('[-o]或[--out]\t参数值不能为空')
            status = True
        if status:
            sys.exit()
        self.path = os.getcwd()
        #self.log.info(self.path)

    def validate(self):
        """
            检查self.samtools软件是否在系统变量路径中，方便后续调用
        """
        status = False
        for ps in os.environ['PATH'].split(':'):
            if os.path.isdir(ps) and self.samtools in os.listdir(ps):
                print('samtools founded  : %s' %ps)
                status = True
        if not status:
            print("samtools is not installed, or not in PATH , install it or set it in PATH")
            exit(2)
        return status

    def process(self):
        """
            执行处理过程，处理传入的输入文件
        """
        self.validate()
        
        if os.path.exists(self.out+self.name+'_R1.fastq.gz') and os.path.getsize(self.out+self.name+'_R1.fastq.gz')>0:
            self.log.warn(self.out+self.name+'_R1.fastq.gz exists,remove it before process.')
            os.remove(self.out+self.name+'_R1.fastq.gz')
        if os.path.exists(self.out+self.name+'_R2.fastq.gz') and os.path.getsize(self.out+self.name+'_R2.fastq.gz')>0:
            self.log.warn(self.out+self.name+'_R2.fastq.gz exists,remove it before process.')
            os.remove(self.out+self.name+'_R2.fastq.gz')
        
        bed = open(self.bed,'r')
        try:
            while True:
                line = bed.readline()
                if line:
                    line  = line.replace('\n', '')
                    line  = line.replace('\n', '')
                    args  = line.split('\t')
                    #print(args)
                    chrom = args[0]
                    start = int(args[1])
                    end   = int(args[2])
                    self.log.info('Processing :  %s     %s\t%s' % (format(chrom," <5"),format(start," >12"),format(end," >12")))
                    if not self.out.endswith('/'):
                        self.out=self.out+'/'
                    
                    with gzip.open(self.out+self.name+'_R1.fastq.gz', "ab")     as R1_out:
                        with gzip.open(self.out+self.name+'_R2.fastq.gz', "ab") as R2_out:
                            try:
                                #读长、测序深度
                                depth  = self.depth
                                lns    = self.length
                                if (end-start+1) > self.length:
                                    depth = int(((end-start+1)/self.length)*self.depth)
                                #print(depth)
                                for index in range(depth):
                                    posi   = self.genPos(start,end,lns)
                                    x      = posi[0]
                                    y      = posi[1]
                                    #print(x,y)
                                    seq = self.genSeq(chrom,x,y,lns)+'\n'
                                    qve = self.genQve(lns,33.00)+'\n'
                                    tit = self.genTit('1')+'\n'
                                    l3  = self.gen3d()+'\n'
                                    R1_out.write(tit)
                                    R1_out.write(seq)
                                    R1_out.write(l3)
                                    R1_out.write(qve)

                                    #经测试，生成reverse seq造成R2数据不正常，mpileup同样为0，不是此处问题
                                    temp = self.genSeq(chrom,x,y,lns)+'\n'
                                    seq2 = self.revSeq(temp)+'\n'
                                    qve2 = self.genQve(lns,32.00)+'\n'
                                    
                                    #此处替换tit的strand参数1为2，bwa mem会计算该坐标是否匹配。
                                    ars  = tit.split(":")
                                    strd = ars[6].split(' ')
                                    six  = strd[0]+' '+'2'
                                    ars[6]=six
                                    tit2 = ":".join(ars)

                                    R2_out.write(tit2)
                                    R2_out.write(seq2)
                                    R2_out.write(l3)
                                    R2_out.write(qve2)
                                    continue
                            finally:
                                R1_out.close()
                                R2_out.close()
                else:
                    break
        finally:
            bed.close()


    def genTit(self,strand='1'):
        """
            Illumina机型生成fastq文件格式
        """        
        instrumentID = 'C70108'                      #设备编号
        runNumber    = '18'                          #Run Number
        flowcellID   = '000000000-AP2P7'             #flow cell ID
        laneID       = '1'                           #lane ID
        tileNumber   = str(random.randint(0,99))     #tile Number
        x            = str(random.randint(0,99999))  #x 坐标
        y            = str(random.randint(0,9999))   #y 坐标
        #strand      = '1'                           #read方向 1read1，2 read2
        filter       = 'N'                           #过滤器N表示通过，Y表示未通过 
        controlNum   = '0'                           #control Number 一般为0
        index        = 'CAGGCAAG+CACGGCGG'           #拆分index 序列
        res= ":".join([
            '@'+instrumentID,
            runNumber,
            flowcellID,
            laneID,
            tileNumber,
            x,
            y+' '+strand,
            filter,
            controlNum,
            index
        ])
        #print(res)
        return res

    
    def genQve(self,len,baseQuality=35.00):
        """
            根据序列长度，生成对应的Q值
        """
        res = ''
        for i in range(len):
            ran = random.random()/50
            qua = baseQuality+33-i*ran
            res+= (chr(int(qua)).encode('ascii'))
        #print(res)
        return res


    def genPos(self,start,end,length):
        """
            根据传入start，end 随机生成start，end坐标
            返回原start，end左右4倍length的随机范围
        """
        overlt = random.randint(-4*length,(length*4+(end-start+1)))
        start  = start-overlt
        end    = start+length-1
        return [start,end]

    def genSeq(self,chrom,start,end,length):
        """
            根据传入坐标，从fa获取reads序列
        """
        cmd = "samtools faidx {ref} {chrom}:{start}-{end}".format(ref=self.ref,chrom=chrom,start=start,end=end)
        seq = self._execute(cmd)
        arr = seq.split("\n")
        res = ''
        for temp in arr:
            if not temp.startswith('>') and len(temp)>0:
                res+=temp
        #print(res)
        return res
    
    
    def revSeq(self,sequence):
        """
            根据已有序列生成反向互补序列
        """
        res = ''
        for tmp in sequence:
            if   tmp=='A':
                res+='T'
            elif tmp=='a':
                res+='t'
            elif tmp=='T':
                res+='A'
            elif tmp=='t':
                res+='a'
            elif tmp=='C':
                res+='G'
            elif tmp=='c':
                res+='g'
            elif tmp=='G':
                res+='C'
            elif tmp=='g':
                res+='c'
            elif tmp=='N':
                res+='N'
            elif tmp=='n':
                res+='n'
        #print(res)
        return res

    
    def gen3d(self):
        #print('+')
        return '+'

    
    def _execute(self,query,workingDir="."):
        """
            调用其他可执行程序，获取输出
        """
        try:
            sub = subprocess.Popen(
                query,
                shell=True,
                cwd=os.path.expanduser(workingDir),
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )
            stdout, stderr = sub.communicate()
            stdout = stdout.decode()
            stderr = stderr.decode()
            if  ((sub.returncode == 0 or sub.returncode == 141)  and
                (stderr == "" or (stderr.startswith("gof3r") and stderr.endswith("broken pipe")))):
                #print(stdout)
                return stdout
            else:
                self.log.error(stderr)
        except Exception as e:
            self.log.error(str(e))


    def usage(self):
        '''
            打印输出Usage
        '''
        print('Usage : ./FastqGenerator.py [OPTION]... 或者 python FastqGenerator.py [OPTION]')
        print('''
            根据输入参考序列Fasta格式文件、bed文件、depth测序深度、len序列长度、输出路径及文件前缀生成模拟的fastq文件 Example:
            
            FastqGenerator.py  -r hg38.fa -b langcancer.bed -d 500 -l 150 -o /opt/result/normal -n Normal

            FastqGenerator.py  --ref=hg38.fa \\
                --bed=/opt/ref/projects/langcancer.bed \\
                --depth=500 \\
                --len=150 \\
                --out=/opt/result/normal \\
                --name=normal
        ''')
        print('''部分使用方法及参数如下：\n''')
        print('-r, --ref=\t参考基因序列文件 .fa或.fasta')
        print('-b, --bed=\t获取序列范围文件 .bed')
        print('-d, --depth=\t生成文件的测序深度')
        print('-l, --len=\t生成文件序列读长')
        print('-n, --name=\t生成Normal样本名称')
        print('-o, --out=\t输出文件目录,输出文件名为{name}_R1.fastq.gz和{name}_R2.fastq.gz')
        print('-h, --help\t显示帮助')
        print('-v, --version\t显示版本号')
        print('--document\t显示开发文档')
        print('\n')
        print('提交bug,to <6041738@qq.com>.\n')


if __name__ == '__main__':
    f=FastqGenerator()
    f.process()