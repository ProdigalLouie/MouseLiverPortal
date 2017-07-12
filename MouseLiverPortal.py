# -*- coding: utf-8 -*-
import flask
from flask import Flask,render_template,url_for,request,jsonify
from sqlalchemy import create_engine,text 
import heapq,copy
from collections import OrderedDict
import math,scipy
from scipy import stats
app= Flask(__name__)
app.jinja_env.variable_start_string = '{{ '
app.jinja_env.variable_end_string = ' }}'
engine=create_engine('mysql://root:mysql122500@localhost/mouse_liver_portal')

@app.route('/compare')
def hello():
    return render_template('compare.html') 

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/browse')
def browse():
    return render_template('browse.html')
    
@app.route('/getData')
def getData():
    sql=text('select * from gene where identified_num>=1')
    conn=engine.connect()
    result=conn.execute(sql)
    genes=[]
    for record in result:
        gene=OrderedDict()
        gene['GeneID']=record[0]
        gene['GeneSymbol']=record[1]
        gene['Chromosome']=record[2]
        gene['Start']=record[3]
        gene['End']=record[4]
        gene['Identified num']=record[5]
        if record[6]==1:
            gene['TF']='Yes'
        else:
            gene['TF']='No'
        genes.append(gene)
    return jsonify({'data':genes})

@app.route('/typePairwise')
def typePairwise():
    #type_gene {'geneID':{'genename':,'median':,'data':}}
    type1,type2=lookup(request.args.get('type0')),lookup(request.args.get('type1'))
    num=200
    geneIDList=set(type1.keys())&set(type2.keys())
    diffHeap=[]
    for geneID in geneIDList:
        if len(diffHeap)<=num:
            heapq.heappush(diffHeap,(abs(type1[geneID]['median']-type2[geneID]['median']),geneID))
        else:
            minDiff=diffHeap[0][0]
            diff=abs(type1[geneID]['median']-type2[geneID]['median'])
            if diff>minDiff:
                heapq.heapreplace(diffHeap,(diff,geneID))
    diffHeap.sort(reverse=True)
    geneIDList=[x[1] for x in diffHeap]
#    geneName=[type1[x]['geneName'] for x in geneIDList]
    topGene1,topGene2=selectGene(geneIDList,type1),selectGene(geneIDList,type2)
    pvalue={}
    for geneID in geneIDList:
        pvalue[geneID]=significance([topGene1[geneID]['data'],topGene2[geneID]['data']])[0]
    return render_template('type.html',type1str=request.args.get('type0'),type2str=request.args.get('type1'),type1=topGene1,type2=topGene2,geneIDList=geneIDList,pvalue=pvalue)

def lookup(typeStr):
    condition=typeStr.split(',')
    strain,gender,genotype,celltype,organelle,treatment,time=condition
    if treatment=='none':
        treatment,time='',''
    sql=text('select genes.geneID,geneName,group_concat(median_normalized_iBAQ) from gene join (select geneID,median_normalized_iBAQ from gene_expression_level_2, (select exp_num from exp_metadata where strain=:strain and gender=:gender and genotype=:genotype and celltype=:celltype and organelle=:organelle and treatment=:treatment and time_after_operation=:time_after_operation and method="profiling") as exps where exps.exp_num=gene_expression_level_2.exp_num) as genes on genes.geneID=gene.geneID group by geneID ')
    conn=engine.connect()
    result=conn.execute(sql,strain=strain,gender=gender,genotype=genotype,celltype=celltype,organelle=organelle,treatment=treatment,time_after_operation=time).fetchall()
    typeData={}
    for record in result:
        #record=[geneID,geneName,group_concat(FOT)]
        typeData[record[0]]={}
        typeData[record[0]]['data']=sorted([math.log10(float(x)) for x in record[2].split(',')])
        typeData[record[0]]['geneName']=record[1]
        length=len(typeData[record[0]]['data'])
        if length%2==0:
            typeData[record[0]]['median']=sum(typeData[record[0]]['data'][length/2-1:length/2+1])/2
        else:
            typeData[record[0]]['median']=typeData[record[0]]['data'][length/2]
    conn.close()
    return typeData
def selectGene(geneIDList,typeData):
    topGene={}
    for geneID in geneIDList:
        topGene[geneID]={}
        topGene[geneID]['geneName']=typeData[geneID]['geneName']
        topGene[geneID]['median']=typeData[geneID]['median']
        topGene[geneID]['data']=typeData[geneID]['data']
    return topGene

@app.route('/searchGene/<keyword>')
def searchGene(keyword):
    sql=text('select geneID,geneName from gene where geneName like :keyword"%"')
    conn=engine.connect()
    result=conn.execute(sql,keyword=keyword).fetchall()
    genes={}
    for record in result:
        genes[record[0]]=record[1]
    conn.close()
    return jsonify(genes)

@app.route('/gene/<int:geneID>/<geneName>')
def gene(geneID,geneName):
    attribute={'strain':{'C57':[],'Balb-c':[]},'gender':{'M':[],'F':[]},'celltype':{'WTE':[],'HC':[],'KC':[],'HSC':[],'LSEC':[]},'organelle':{'WCE':[],'LDP':[],'NE':[],'Mitochondria':[],'Membrane Proteins':[]},'time_after_operation':{'0':[],'12':[],'24':[],'48':[],'72':[],'120':[],'168':[]}}
    pvalue={}
    for attributeKey in attribute:
        fot=getFot(geneID,attributeKey,attribute[attributeKey].keys())
        attribute[attributeKey]=fot['data']       
        num=len(attribute[attributeKey].keys())
        pvalue[attributeKey]=[[None for col in range(num)] for row in range(num)]
        idx=0
        for i in range(num):
            for j in range(i+1,num):
                pvalue[attributeKey][i][j]=fot['pvalue'][idx]
                idx+=1

        #pvalue[attributeKey]=fot['pvalue']
    return render_template('gene.html',geneID=geneID,geneName=geneName,attribute=attribute,pvalue=pvalue)
    #return jsonify(pvalue)

def getFot(geneID,attribute,value):
    defaultCondition={'geneID':geneID,'strain':'C57','gender':'M','genotype':'WT','celltype':'WTE','organelle':'WCE','treatment':'','time_after_operation':''}
    if attribute=='time_after_operation':
        defaultCondition['treatment']='PHx'
    conn=engine.connect()
    fot={'data':{}}
    for x in value:
        condition=copy.deepcopy(defaultCondition)
        condition[attribute]=x
        sql=text('select median_normalized_iBAQ from gene_expression_level_2, (select exp_num from exp_metadata where strain=:strain and gender=:gender and genotype=:genotype and celltype=:celltype and organelle=:organelle and treatment=:treatment and time_after_operation=:time_after_operation and method="profiling") as exps where geneID=:geneID and gene_expression_level_2.exp_num=exps.exp_num')
        result=conn.execute(sql,condition).fetchall()
        fot['data'][x]=[]
        for record in result:
            fot['data'][x].append(math.log10(record[0]))
    conn.close()
    fot['pvalue']=significance(fot['data'].values())
    return fot

def significance(data):
    #test for normality
    num=sum(range(1,len(data)))
    sig=[None]*num
    norm_flag=[True]*len(data)
    for idx,group in enumerate(data):
        if len(group)<3:
            norm_flag[idx]=False
        elif 3<=len(group)<=50:
            if scipy.stats.shapiro(group)[1]<0.05:
                norm_flag[idx]=False
        else:
            if scipy.stats.kstest(group,'norm')[1]<0.05:
                norm_flag[idx]=False
    #test significance of difference
    length=len(data)
    idx=0
    for i in range(length):
        for j in range(i+1,length):
            if len(data[i])>=3 and len(data[j])>=3:
                if norm_flag[i] and norm_flag[j]:
                    sig[idx]=scipy.stats.ttest_ind(data[i],data[j],equal_var=scipy.stats.bartlett(data[i],data[j])[1]>0.05)[1]
                else:
                    if len(data[i])<=20 or len(data[j])<=20:
                        sig[idx]=stats.ranksums(data[i],data[j])[1]
                    else:
                        sig[idx]=stats.mannwhitneyu(data[i],data[j])[1]
            idx+=1
    return sig

@app.template_filter('get_color_hex')
def get_color_hex(pvalue):
    if pvalue<0.01:
        return '#ff3e3e'
    elif pvalue<0.05:
        return '#ff9999'
    else:
        return '#fff'

@app.route('/guide')
def guide():
    return render_template('guide.html')
@app.route('/test')
@app.route('/test/<a1>/<a2>')
def test(a1='1',a2='2'):
    '''
    type1={}
    type1[2]={}
    type1[2]['geneName']='ert'
    type1[2]['data']=[1,2,3]
    type1[2]['median']='ert'
    return render_template('test.html',type1=type1,geneIDList=[2])
    '''
    return render_template('test.html')
@app.route('/test1')
def test1():
    return render_template('test1.html')
if __name__=='__main__':
    app.run(debug=True)

