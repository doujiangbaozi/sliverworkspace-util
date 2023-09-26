#!/usr/bin/env python
#-*-coding:utf-8-*-

__author__ = '豆浆包子'
__all__    = ['bcl2fq', 'downSampleSheet', 'execute','usage']


import os,sys,csv,getopt,io,requests,json,subprocess,time,tempfile
from multiprocessing import cpu_count
from configparser    import ConfigParser
reload(sys)
sys.setdefaultencoding("utf-8")


class bcl2fq(object):
    """
        对接SliverWorkspace系统的自动化数据分拆工具
    """

    def __init__(self):
        """
        初始化,获取运行参数
        """
        try:
            #print('__file__:', __file__)
            #print('realpath of __file__:', os.path.realpath(__file__))
            #print('sys.executable:', sys.executable)
            #print('realpath of sys.executable:', os.path.realpath(sys.executable)[:os.path.realpath(sys.executable).rfind('/')])
            #print('sys.argv[0]:', sys.argv[0])
            #print('realpath of sys.argv[0]:', os.path.realpath(sys.argv[0]))
            #print('sys.path[0]:', sys.path[0])
            #print('realpath of sys.path[0]:', os.path.realpath(sys.path[0]))
            cfg = ConfigParser()
            
            if os.path.exists(os.path.dirname(os.path.realpath(sys.argv[0]))+os.sep+'configuration.cfg'):
                cfg.read(os.path.dirname(os.path.realpath(sys.argv[0]))+os.sep+'configuration.cfg')
            elif os.path.exists(os.path.realpath(sys.path[0])+os.sep+'configuration.cfg'):
                cfg.read(os.path.realpath(sys.path[0])+os.sep+'configuration.cfg')
            elif os.path.exists(os.path.realpath(sys.executable)[:os.path.realpath(sys.executable).rfind('/')]+os.sep+'configuration.cfg'):
                cfg.read(os.path.realpath(sys.executable)[:os.path.realpath(sys.executable).rfind('/')]+os.sep+'configuration.cfg')
            self.data       = {
                'userName':cfg.get('USER','username'),
                'userPass':cfg.get('USER','password')
            }
            self.headers    = {
                'User-agent':cfg.get('URL', 'User-agent')
            }
            self.url_base   = cfg.get('URL','host')
            self.url_login  = self.url_base+cfg.get('URL','login')
            self.url_logout = self.url_base+cfg.get('URL','logout')
            self.url_down   = self.url_base+cfg.get('URL','download')
            self.url_update = self.url_base+cfg.get('URL','update')
            self.verify     = cfg.getboolean('URL','verify')
            self.args       = {}
            self.opts       = getopt.getopt(
                sys.argv[1:],
                "R:o:i:r:p:w:hvl:", 
                [
                    "runfolder-dir=",
                    "output-dir=",
                    "input-dir=",
                    "loading-threads=",
                    "processing-threads=",
                    "writing-threads=",
                    "help",
                    "version",
                    "min-log-level=",
                    "sample-sheet=",
                    "intensities-dir=",
                    "interop-dir=",
                    "stats-dir=",
                    "reports-dir=",
                    "adapter-stringency=",
                    "barcode-mismatches=",
                    "create-fastq-for-index-reads",
                    "ignore-missing-bcls",
                    "ignore-missing-filter",
                    "ignore-missing-positions",
                    "minimum-trimmed-read-length=",
                    "mask-short-adapter-reads=",
                    "tiles=" ,
                    "use-bases-mask=",
                    "with-failed-reads",
                    "write-fastq-reverse-complement",
                    "no-bgzf-compression",
                    "fastq-compression-level=",
                    "no-lane-splitting",
                    "find-adapters-with-sliding-window",
                    "id=",
                    "document"
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
        threads = cpu_count()
        if threads>16:
            threads = 16
        for arg, value in self.opts[0]:
            if arg=="-R" or arg == "--runfolder-dir":
                if str.endswith(value.strip(),'/'):
                    value = value[:len(value)-1]
                if os.path.isdir(value.strip()) :
                    self.args['runfolder-dir']=value.strip()
                else:
                    print('[-R]或[--runfolder-dir]参数值 %s 目录不存在' %(value))
                    sys.exit()
            elif arg=="-o" or arg == "--output-dir":
                if str.endswith(value.strip(),'/'):
                    value = value[:len(value)-1]
                if len(value.strip())>0 :
                    self.args['output-dir']=value.strip()
                else:
                    self.args['output-dir']='.'
                    #print('[-o]或[--output-dir]参数值 %s 不能为空' %(value))
                    #sys.exit()
            elif arg=="-i" or arg == "--input-dir":
                if str.endswith(value.strip(),'/'):
                    value = value[:len(value)-1]
                if os.path.isdir(value.strip()) :
                    self.args['input-dir']=value.strip()
                else:
                    print('[-i]或[--input-dir]参数值 %s 目录不存在' %(value))
                    sys.exit()
            elif arg=="-r" or arg=="--loading-threads":
                if value.isdigit() and int(value)<=threads:
                    self.args['loading-threads']=value
                else:
                    self.args['loading-threads']=str(threads)
            elif arg=="-p" or arg=="--processing-threads":
                if value.isdigit() and int(value)<=threads:
                    self.args['processing-threads']=value
                else:
                    self.args['processing-threads']=str(threads)
            elif arg=="-w" or arg=="--writing-threads":
                if value.isdigit() and int(value)<=threads:
                    self.args['writing-threads'] = value
                else:
                    self.args['writing-threads'] = str(threads)
            elif arg=="--sample-sheet":
                if len(value.strip()):
                    self.args['sample-sheet'] = value
                else:
                    print('[--sample-sheet]参数值 %s 不能为空' %(value))
                    sys.exit()
            elif arg=="--no-lane-splitting":
                self.args['no-lane-splitting'] = '--no-lane-splitting'
            elif arg=="-h" or arg == "--help":
                self.usage()
                sys.exit()
            elif arg=="-v" or arg == "--version":
                print("版本: 1.00")
                exit()
            elif arg=="barcode-mismatches":
                if value.isdigit() and (int(value)==0 or int(value)==1 or int(value)==2):
                    self.args['barcode-mismatches']=value
                else:
                    print('[barcode-mismatches]参数值只能是: 0,1,2；默认值为1')
                    sys.exit()
            elif arg=="--id":
                if(len(value.strip())>0):
                    self.args['id']=value
                else:
                    print('[--id]参数值 %s 没有有效值' %(value))
                    sys.exit()
            elif arg=="--document":
                from bcl2fq import bcl2fq
                help(bcl2fq)
                sys.exit()
        if not self.args.has_key('id'):
            print('\n[--id]参数没有有效值\n')
            sys.exit()
        if not self.args.has_key('runfolder-dir'):
            self.args['runfolder-dir']='.'
        if not self.args.has_key('output-dir'):
            self.args['output-dir'] = self.args['runfolder-dir']+'/Data/Intensities/BaseCalls/'
        if not self.args.has_key('input-dir'):
            self.args['input-dir']  = self.args['runfolder-dir']+'/Data/Intensities/BaseCalls/'
        if not self.args.has_key('no-lane-splitting'):
            self.args['no-lane-splitting']=''
        if not self.args.has_key('loading-threads'):
            self.args['loading-threads']='8'
        if not self.args.has_key('processing-threads'):
            self.args['processing-threads']='8'
        if not self.args.has_key('writing-threads'):
            self.args['writing-threads']='8'
        if not self.args.has_key('barcode-mismatches'):
            self.args['barcode-mismatches']='1'
        self.query = '''
            bcl2fastq -r %(loading-threads)s -p %(processing-threads)s -w %(writing-threads)s %(no-lane-splitting)s \
                --barcode-mismatches %(barcode-mismatches)s \
                --runfolder-dir %(runfolder-dir)s \
                --sample-sheet  %(sample-sheet)s \
                --output-dir    %(output-dir)s
        '''


    def execute(self,workingDir="."):
        '''
            本地运行bcl2fastq
        '''
        start = time.time()
        try:
            query = self.query % self.args
            print(query)
            # 得到一个临时文件对象， 调用close后，此文件从磁盘删除
            #out_temp = tempfile.TemporaryFile(mode='w+')
            # 获取临时文件的文件号
            #fileno = out_temp.fileno()
            sub = subprocess.Popen(query,shell=True,cwd=os.path.expanduser(workingDir),stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
            while sub.poll() is None:
                line = sub.stdout.readline()
                line = line.strip()
                if line:
                    print(line)
                if sub.returncode == 0:
                    print('Demultiplex Run Succeed')
            #sub.wait()
            #out_temp.seek(0)
            #print(out_temp.read())
        except Exception as e:
            print(str(e))
        #finally:
            #if out_temp:
            #    out_temp.close()
        stop  = time.time()
        print('Run time:'+str(round((stop-start),3))+' seconds')
    


    def process(self):
        ''' 
            1、从web登录SliverWorkspace系统，获取Session ID之后下载SampleSheet 
            2、本地运行bcl2fastq，拆分数据
            3、Web发送请求，更新系统信息，标记已经拆分过数据的Sample
            4、远程退出账户，完成操作
        '''
        #构造Session
        session = requests.Session()
        #在session中发送登录请求，此后这个session里就存储了cookie
        #可以用print(session.cookies.get_dict())查看
        res = session.post(
            self.url_login,
            self.data,
            verify  = self.verify,
            headers = self.headers
        )
        state = json.loads(res.content)
        if state['status'] and state['status']==u'SUCCESS':
            #发送访问请求
            print('成功登录系统')
            res = session.get(
                self.url_down+'?id='+self.args['id'], 
                verify  = self.verify,
                headers = self.headers
            )
            print(res.headers)
            if (res.headers['Content-Disposition']) and res.headers['Content-Disposition']=='attachment;fileName=SampleSheet.csv':
                if not self.args.has_key('sample-sheet'):
                    self.args['sample-sheet']=self.args['runfolder-dir']+os.sep+'SampleSheet.csv'
                with open(self.args['sample-sheet'], "wb") as code:
                    code.write(res.content)
                if os.path.exists(self.args['sample-sheet']):
                    self.execute()
                    result = self.args['output-dir']+os.sep+'Stats'+os.sep+'Stats.json'
                    if(os.path.exists(result)):
                        with open(result) as file_object:
                            data = json.load(file_object)
                            ps   = {
                                'sampleRunId':self.args['id'],
                                'state':'true',
                                'stat':json.dumps(data)
                            }
                            print(ps)
                            res = session.post(
                                self.url_update,
                                ps,
                                verify  = self.verify,
                                headers = self.headers
                            )
        else:
            print('登录失败，请检查账户密码是否正确')
        res = session.get(
            self.url_logout,
            verify  = self.verify,
            headers = self.headers
        )
        state = json.loads(res.content)
        if state['status'] and state['status']==u'SUCCESS':
            print('成功退出系统')




    def usage(self):
        '''
        打印输出Usage
        '''
        print('\nUsage : ./bcl2fq.py [OPTION]... 或者 python bcl2fq.py [OPTION]  编译版本：bcl2fq [OPTION]')
        print('''
            用于配合系统SliverWorkspace自动拆分样本数据，并更新软件中样本信息
            使用参数与illumina官方工具bcl2fastq一致
        ''')
        print('''部分使用方法及参数如下：\n''')
        print('\n')
        print('   --id\t\t\t SliverWorkspace系统里对应于RUN ID字段,缺少脚本将无法正常运行')
        print('\n')
        print("-R --runfolder-dir\t Path to run folder directory Default:./")
        print("-o --output-dir\t\t Path to demultiplexed output Default:<runfolder-dir>/Data/Intensities/BaseCalls")
        print("-i --input-dir\t\t Path to input directory Default:<runfolder_dir>/Data/Intensities/BaseCalls")
        print("-r --loading-threads\t Number of threads used for loading BCL data.Default depends on architecture.")
        print("-p --processing-threads\t Number of threads used for processing demultiplexed data.Default depends on architecture.")
        print("-w --writing-threads\t Number of threads used for writing FASTQ data. This number should not be set higher than number of samples.Default depends on architecture")
        print("-h --help\t\t help document")
        print("-v --version\t\t display version")
        print('   --document\t\t display documents')
        print("   --sample-sheet\t Path to sample sheet,so you can specity the location and name of the sample sheet,if different from default.Default:<runfolder_dir>/SampleSheet.csv")
        print("   --no-lane-splitting\t Do not split FASTQ files by lane.")
        print("   --barcode-mismatches\t Number of allowed mismatches per index Multiple entries, comma delimited allowed. Each entry is applied to the corresponding index; last entry applies to all remaining indexes.Default: 1. Accepted values: 0, 1 or 2.")
        print('\n')
        print('提交bug,to <6041738@qq.com>.\n')

if __name__=="__main__":
    bcl2fq = bcl2fq()
    bcl2fq.process()