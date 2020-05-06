

#__author__ = "ilker Tunc"
# NHLBI BCB Core
# Mar30, 2020
# edits : May,3,2020

import scipy.stats as stats
from collections import defaultdict
import pandas as pd 
import numpy as np 
import argparse
import sys

from rpy2 import robjects




# could use python for that, but R much better for stats
robjects.r(
"""
ttest <- function(x1,x2,cutoff=0.2)
  {

       ss=var.test(x1,x2)       
       if (is.na(ss$p.value)) {
         return(NULL)
       }
       if (ss$p.value > cutoff) 
         {
            zz=t.test(x1,x2,var.equal = TRUE)
            return(zz$p.value)
          }
       else {
          zz=t.test(x1,x2,var.equal = FALSE)
          return(zz$p.value)
          }
}
""")




r_ttest=robjects.r['ttest']



def get_score(df,gname):
    vv=df.loc[idx[gname,:],:]
    vv=vv.reset_index()
    vv=vv.set_index('segment')
    ss=vv.loc[vv.index[1:-1],'score'].sum()
    if ss>0:
        return 4
    if ss<0:
        return -4
    if ss==0:
        return 0 




def get_sum_verdict(df,gname):
    vv=df.loc[idx[gname,:],:]
    vv=vv.reset_index()
    vv=vv.set_index('segment')
    ss=vv.loc[vv.index[0],'score']+vv.loc[vv.index[-1],'score']
    return ss




# designate if first, last, only?
def mark_rows(df):
    for gname in df.index.get_level_values('GeneID'):
        zz=sorted(df.loc[(gname,),:].index)
        if len(zz)==1:
            df.loc[(gname,1),'is_first']=1
            df.loc[(gname,1),'is_only']=1
            df.loc[(gname,1),'is_last']=1
        else:
            df.loc[(gname,1),'is_first']=1  
            df.loc[(gname,zz[-1]),'is_last']=1      

    df['is_first']=df['is_first'].fillna(0)
    df['is_last']=df['is_last'].fillna(0)
    df['is_only']=df['is_only'].fillna(0)



# get the segment number before segment i
def get_pre_segment(df,gname,i):
    xx=df.loc[(gname,),:].index.tolist()
    return xx[xx.index(i)-1]


# get the segment number after segment i
def get_post_segment(df,gname,i):
    xx=df.loc[(gname,),:].index.tolist()
    return xx[xx.index(i)+1]


def do_ttest(x1,x2):
    try:

        if (x1.notna().sum())<2:
            return None

        if (x2.notna().sum())<2:
            return None

        pval=r_ttest(robjects.FloatVector(x1),robjects.FloatVector(x2))
        return pval[0]

    except:
        #print(sys.exc_info())
        #print(x1,x2,pval[0])
        return None





parser = argparse.ArgumentParser(description='APA analysis')
parser.add_argument('--table', '-t',dest="table",help='count table for all samples together, \
    GeneID and segment must be in the columns, this option is ignored if files_path specified')
parser.add_argument('--files_path', '-f',dest="files_path",help='make table from the given directory of sample files  ')
parser.add_argument('--group_control', '-g0',dest='control_group',help='comma seperated list of samples in control group')
parser.add_argument('--group_treatment', '-g1',dest='treatment_group',help='comma seperated list of samples in treatment group')
parser.add_argument('--label_control', '-l0',default='Control',dest='label_control',help='label for control group')
parser.add_argument('--label_treatment', '-l1',default='Treatment',dest='label_treatment', help='label for treatment group')
parser.add_argument('--out', '-o',dest='outfile', help='outfile for main results')
parser.add_argument('--out_main', '-x',dest='outfile2', help='outfile for main table before filtering and analysis')
parser.add_argument('--score_single',default=0, type=int,dest='score_single', help='score for single')
parser.add_argument('--score_first',default=3, type=int,dest='score_first', help='score for first')
parser.add_argument('--score_last',default=4, type=int,dest='score_last', help='score for last')
parser.add_argument('--score_seg',default=1, type=int,dest='score_seg', help='score for segment')
parser.add_argument('--score_seg_change',default=3, type=int,dest='score_seg_change', help='score for segment change')
parser.add_argument('--mean_filter',default=0.5, type=float,dest='mean_filter', help='mean filter')
parser.add_argument('--missing_filter',default=2, type=int,dest='missing_filter', help='missing filter')
parser.add_argument('--pval_filter',default=0.05, type=float,dest='pval_filter', help='pval filter')
parser.add_argument('--cutoff',default=0.2, type=float,dest='cutoff', help='cutoff used if number of replicates less than 3 and no t-test is applied')


args = parser.parse_args()


try:
    g0=args.control_group.split(',')
    g1=args.treatment_group.split(',')
except:
    print("error in group sample names")



if args.files_path:

    print("reading sample files...")
    
    df_list=[]
    for g in g0+g1:
        df=pd.read_csv('%s/%s_abundance.tsv' %(args.files_path,g),sep='\t',usecols=['target_id','tpm'],index_col=0)
        df.rename(columns={'tpm':g},inplace=True)
        df_list.append(df)
        print("read the sample file for %s" %g)


    print("merging sample files...")

    df=pd.concat(df_list,axis=1)
    df=df.reset_index()
    df['GeneID']=df['target_id'].str.split("|",expand=True)[0]
    df['GeneName']=df['target_id'].str.split("|",expand=True)[1]
    df['leftcoord']=df['target_id'].str.split("|",expand=True)[2]
    df['rightcoord']=df['target_id'].str.split("|",expand=True)[3]
    df['pAsite']=df['target_id'].str.split("|",expand=True)[4]
    df['TC']=df['target_id'].str.split("|",expand=True)[5]
    df['usage']=df['target_id'].str.split("|",expand=True)[6]
    df['strand']=df['target_id'].str.split("|",expand=True)[7]
    df['segment']=df['strand'].str.split("(",expand=True)[0]
    df['segment']=df['segment'].astype(int)
    df['strand']=df['strand'].str.extract(r'\(([-+])\)',expand=True)
       
    # segment sorting seems working, otherwise use this to be sure sorted correct for negativs strand 
    #df.loc[df['strand']=="-","leftsort"]=-df.loc[df['strand']=="-","leftcoord"]
    #df.loc[df['strand']=="+","leftsort"]=df.loc[df['strand']=="+","leftcoord"]
    #df=df.sort_values(['GeneID','leftsort'])

    df=df.set_index(['GeneID','segment'])
    df=df.sort_index()


    print("merging finished, total reads=%s" %len(df))


elif args.table:
    print("reading main data table...")
    try:
        df=pd.read_csv(args.table,sep='\t',index_col=['GeneID','segment'])
        df=df.sort_index()
        print("main table read, total reads=%s" %len(df))
    except:
        print("Problem with opening the datafile, data table is required to have GeneID and segment as columns")

else:
    print("no table file or files path specified")
    sys.exit(0)



idx=pd.IndexSlice


lb0=args.label_control
lb1=args.label_treatment


mean_filter=args.mean_filter
missing_filter=args.missing_filter
pval_filter=args.pval_filter

score_single=args.score_single
score_first=args.score_first
score_last=args.score_last
score_seg=args.score_seg
score_seg_change=args.score_seg_change

cutoff=args.cutoff 

outfile=args.outfile
outfile2=args.outfile2



dct_catg={
    'single_1': "no APA, single seg more in %s" %lb1,
    'single_0': "no APA, single seg more in %s" %lb0,
    'single_eq': "no APA, single seg equal amts",
    'first_1': "first seg, more in %s" %lb1,
    'first_0': "first seg, more in %s" %lb0,
    'first_eq': "first seg, equal amts",
    'last_1': "last segment, more in %s" %lb1,
    'last_0': "last segment, more in %s" %lb0,
    'last_eq': "last seg, same change",
    'seg_eq': "this seg, same change",
    'seg_1': "more in %s" %lb1,
    'seg_0': "more in %s" %lb0,
    'seg_change_1': "changepoint to more in %s" %lb1,
    'seg_change_0': "changepoint to more in %s" %lb0,
    'seg_fluc': "fluc",
}


dct_score={
    'single_1': score_single,
    'single_0': -score_single,
    'single_eq': 0,
    'first_1': score_first,
    'first_0': -score_first,
    'first_eq': 0,
    'last_1': score_last,
    'last_0': -score_last,
    'last_eq': False,
    'seg_eq': 0,
    'seg_1': score_seg,
    'seg_0': -score_seg,
    'seg_change_1': score_seg_change,
    'seg_change_0': -score_seg_change,
    'seg_fluc': 0,
}


print("finding log2 and group means...")


for c in g0+g1:
    colnamenew='log2tpm%s' %c
    # ignore warning of log2of 0 givin -inf, later replaced by None
    with np.errstate(divide='ignore'):
        df[colnamenew]=np.log2(df[c])


df=df.replace(-np.inf,np.nan)

mark_rows(df)

df['meantpm%s' %lb0]=df.loc[:,g0].mean(axis=1)
df['meantpm%s' %lb1]=df.loc[:,g1].mean(axis=1)
df['meanlog2tpm%s' %lb0]=df.loc[:,['log2tpm%s' %c for c in g0]].mean(axis=1)
df['meanlog2tpm%s' %lb1]=df.loc[:,['log2tpm%s' %c for c in g1]].mean(axis=1)
df['difflog2tpm']=df.loc[:,'meanlog2tpm%s' %lb1]-df.loc[:,'meanlog2tpm%s' %lb0]
#df['%s-%slogtpmmean' %(lb1,lb0)]=df['meanlog2tpm%s' %lb1]-df['meanlog2tpm%s' %lb0]


print("writing the count table before analysis...")
df.to_csv('%s_%s_vs_%s.csv' %(outfile2,lb1,lb0))


df_filtered=None



ii=df.query('is_first!=1 and meantpm%s<%f' %(lb0,mean_filter)).index
tmp=df.loc[ii,]
tmp['FILTER']=None
df_filtered=tmp
df_filtered.loc[ii,"FILTER"]="not first and mean TPM of %s<%f" %(lb0,mean_filter)
df.drop(ii,inplace=True)



ii=df.query('is_only==1 and meantpm%s<%f' %(lb0,mean_filter)).index
tmp=df.loc[ii,]
tmp['FILTER']=None
df_filtered=df_filtered.append(tmp)
df_filtered.loc[ii,"FILTER"]="single and mean TPM of %s<%f" %(lb0,mean_filter)
df.drop(ii,inplace=True)


print("filtering meanTPM, reads left=%s" %len(df))


vv=df.loc[df['is_first']!=1,['log2tpm%s' %c for c in g0]].isna().sum(axis=1)>=missing_filter
ii=vv[vv==True].index
tmp=df.loc[ii,]
tmp['FILTER']=None
df_filtered=df_filtered.append(tmp)
df_filtered.loc[ii,"FILTER"]="not first and number of missing in group %s >= %s" %(lb0,missing_filter)
df.drop(ii,inplace=True)


vv=df.loc[df['is_first']!=1,['log2tpm%s' %c for c in g1]].isna().sum(axis=1)>=missing_filter
ii=vv[vv==True].index
tmp=df.loc[ii,]
tmp['FILTER']=None
df_filtered=df_filtered.append(tmp)
df_filtered.loc[ii,"FILTER"]="not first and number of missing in group %s >= %s" %(lb1,missing_filter)
df.drop(ii,inplace=True)


print("filtering missing log2TPM, reads left=%s" %len(df))

# redo marking since some rows will be filtered
mark_rows(df)



# redo filtering for the segments that become single after the  remarking
ii=df.query('is_only==1 and meantpm%s<%f' %(lb0,mean_filter)).index
tmp=df.loc[ii,]
tmp['FILTER']=None
df_filtered=df_filtered.append(tmp)
df_filtered.loc[ii,"FILTER"]="after re-marking, single and mean TPM of %s<%f" %(lb0,mean_filter)
df.drop(ii,inplace=True)


print("filtering meanTPM after re-marking, reads left=%s" %len(df))



print("finding normalized differences...")

for gname in df.index.get_level_values('GeneID'):
    for c in g0+g1:
        # idx without second colon, it only returns second index, hides the given gene name
        for i in sorted(df.loc[idx[gname,],:].index):
            colname='log2tpm%s' %c
            colnamenew='%snormdiff' %c
            df.at[(gname,i),colnamenew]=df.at[(gname,i),colname]-df.at[(gname,1),colname]
           
        

print("finding segment differences...")

for gname in df.index.get_level_values('GeneID'):
    for c in g0+g1:
        tt=sorted(df.loc[idx[gname,],:].index)
        if len(tt)>1:
            for i,j in zip(tt,tt[1:]):
                colname='log2tpm%s' %c
                colnamenew='%ssegdiff' %c
                df.at[(gname,j),colnamenew]=df.at[(gname,j),colname]-df.at[(gname,i),colname] 
            # first row assign 0
            df.at[(gname,1),colnamenew]=0
        else:
            df.at[(gname,1),colnamenew]=0




df['%snormdiffmean' %lb0]=df.loc[:,['%snormdiff' %c for c in g0]].mean(axis=1)
df['%snormdiffmean' %lb1]=df.loc[:,['%snormdiff' %c for c in g1]].mean(axis=1)
df['%ssegdiffmean' %lb0]=df.loc[:,['%ssegdiff' %c for c in g0]].mean(axis=1)
df['%ssegdiffmean' %lb1]=df.loc[:,['%ssegdiff' %c for c in g1]].mean(axis=1)

df['%s-%snormdiff' %(lb1,lb0)]=df['%snormdiffmean' %lb1]-df['%snormdiffmean' %lb0]
df['%s-%ssegdiff' %(lb1,lb0)]=df['%ssegdiffmean' %lb1]-df['%ssegdiffmean' %lb0]





#col_logtpmmean='%s-%slogtpmmean' %(lb1,lb0)
col_logtpmmean='difflog2tpm'
col_normdiff='%s-%snormdiff' %(lb1,lb0)
col_segdiff='%s-%ssegdiff' %(lb1,lb0)


if len(g0)>2 and len(g1)>2:

    print("performing t-test...")

    for i,r in df.iterrows():


        z1=['log2tpm%s' %c for c in g0]
        z2=['log2tpm%s' %c for c in g1]
        df.loc[i,'pvallogtpm']=do_ttest(r[z1],r[z2])


        z1=['%snormdiff' %c for c in g0]
        z2=['%snormdiff' %c for c in g1]
        df.loc[i,'pvalnormdiff']=do_ttest(r[z1],r[z2])


        z1=['%ssegdiff' %c for c in g0]
        z2=['%ssegdiff' %c for c in g1]
        df.loc[i,'pvalsegdiff']=do_ttest(r[z1],r[z2])




    print("calculating the scores...")


    for i,r in df.iterrows():



        if pd.isna(r['pvallogtpm']):
            continue

        value_logtpm=r[col_logtpmmean]
        value_normdiff=r[col_normdiff]
        value_segdiff=r[col_segdiff]

        if r['is_only']:
            if r['pvallogtpm']<pval_filter:
                if value_logtpm>0:
                    df.loc[i,"eval"]=dct_catg['single_1']
                    df.loc[i,"score"]=dct_score['single_1']
                else:
                    df.loc[i,"eval"]=dct_catg['single_0']
                    df.loc[i,"score"]=dct_score['single_0']             
            else:
                df.loc[i,"eval"]=dct_catg['single_eq']
                df.loc[i,"score"]=dct_score['single_eq']
        elif r['is_first']:
            if r['pvallogtpm']<pval_filter:
                if value_logtpm>0:
                    df.loc[i,"eval"]=dct_catg['first_1']
                    df.loc[i,"score"]=dct_score['first_1']              
                else:
                    df.loc[i,"eval"]=dct_catg['first_0']
                    df.loc[i,"score"]=dct_score['first_0']              
            else:
                df.loc[i,"eval"]=dct_catg['first_eq']
                df.loc[i,"score"]=dct_score['first_eq']         
        elif r['is_last']:
            if pd.isna(r['pvalnormdiff']):
                continue
            if r['pvalnormdiff']<pval_filter:
                if value_normdiff>0:
                    df.loc[i,"eval"]=dct_catg['last_1']
                    df.loc[i,"score"]=dct_score['last_1']               
                else:
                    df.loc[i,"eval"]=dct_catg['last_0']
                    df.loc[i,"score"]=dct_score['last_0']
            else:
                df.loc[i,"eval"]=dct_catg['last_eq']        
        else:
            if pd.isna(r['pvalnormdiff']):
                continue
            if r['pvalnormdiff']<pval_filter:
                if value_normdiff >0:
                    if pd.isna(r['pvalsegdiff']):
                        continue
                    if r['pvalsegdiff']<pval_filter:
                        if value_segdiff >0:
                            df.loc[i,"eval"]=dct_catg['seg_change_1']
                            df.loc[i,"score"]=dct_score['seg_change_1']                     
                        else:
                            df.loc[i,"eval"]=dct_catg['seg_fluc']
                            df.loc[i,"score"]=dct_score['seg_fluc']                     
                    else:
                        df.loc[i,"eval"]=dct_catg['seg_1']
                        df.loc[i,"score"]=dct_score['seg_1']                    
                else:
                    if pd.isna(r['pvalsegdiff']):
                        continue
                    if r['pvalsegdiff']<pval_filter:
                        if value_segdiff <=0:
                            df.loc[i,"eval"]=dct_catg['seg_change_0']
                            df.loc[i,"score"]=dct_score['seg_change_0']                     
                        else:
                            df.loc[i,"eval"]=dct_catg['seg_fluc']
                            df.loc[i,"score"]=dct_score['seg_fluc']                     
                    else:
                        df.loc[i,"eval"]=dct_catg['seg_0']
                        df.loc[i,"score"]=dct_score['seg_0']                    
            else:
                df.loc[i,"eval"]=dct_catg['seg_eq']
                df.loc[i,"score"]=dct_score['seg_eq']
else:

    print("Number of samples less than 3, so skipping t-test, using cut-off")


    for i,r in df.iterrows():

        value_logtpm=r[col_logtpmmean]
        value_normdiff=r[col_normdiff]
        value_segdiff=r[col_segdiff]

        if r['is_only']:
            if value_logtpm >= cutoff:
                df.loc[i,"eval"]=dct_catg['single_1']
                df.loc[i,"score"]=dct_score['single_1']
            elif value_logtpm <= -cutoff:
                df.loc[i,"eval"]=dct_catg['single_0']
                df.loc[i,"score"]=dct_score['single_0'] 
            else:
                df.loc[i,"eval"]=dct_catg['single_eq']
                df.loc[i,"score"]=dct_score['single_eq']        
        elif r['is_first']:
            if value_logtpm >= cutoff:
                df.loc[i,"eval"]=dct_catg['first_1']
                df.loc[i,"score"]=dct_score['first_1']              
            elif value_logtpm <= -cutoff:
                df.loc[i,"eval"]=dct_catg['first_0']
                df.loc[i,"score"]=dct_score['first_0']              
            else:
                df.loc[i,"eval"]=dct_catg['first_eq']
                df.loc[i,"score"]=dct_score['first_eq']         
        elif r['is_last']:
            if value_normdiff >= cutoff:
                df.loc[i,"eval"]=dct_catg['last_1']
                df.loc[i,"score"]=dct_score['last_1']               
            elif value_normdiff <= -cutoff:
                df.loc[i,"eval"]=dct_catg['last_0']
                df.loc[i,"score"]=dct_score['last_0']
            else:
                df.loc[i,"eval"]=dct_catg['last_eq']        
        else:
            genename=i[0]
            curr_seg=i[1]
            pre_seg=get_pre_segment(df,genename,curr_seg)
            post_seg=get_post_segment(df,genename,curr_seg)
            try:
                if r['strand']=="-":
                    value_compare=df.at[(genename,pre_seg),col_segdiff]
                elif r['strand']=="+":
                    value_compare=df.at[(genename,post_seg),col_segdiff]
                else:
                    print("error, strand information missing")
                    print(r)

                if value_segdiff <0 and value_compare >= -value_segdiff:
                    df.loc[i,"eval"]=dct_catg['seg_fluc']
                    df.loc[i,"score"]=dct_score['seg_fluc'] 
                elif value_segdiff >0 and value_compare  <= -value_segdiff:
                    df.loc[i,"eval"]=dct_catg['seg_fluc']
                    df.loc[i,"score"]=dct_score['seg_fluc'] 
                elif value_normdiff <= -cutoff:
                    if value_segdiff <= -cutoff:
                        df.loc[i,"eval"]=dct_catg['seg_change_0']
                        df.loc[i,"score"]=dct_score['seg_change_0'] 
                    else :
                        df.loc[i,"eval"]=dct_catg['seg_0']
                        df.loc[i,"score"]=dct_score['seg_0'] 
                elif value_normdiff >= cutoff:
                    if value_segdiff >= cutoff:
                        df.loc[i,"eval"]=dct_catg['seg_change_1']
                        df.loc[i,"score"]=dct_score['seg_change_1'] 
                    else :
                        df.loc[i,"eval"]=dct_catg['seg_1']
                        df.loc[i,"score"]=dct_score['seg_1'] 
                else:
                        df.loc[i,"eval"]=dct_catg['seg_eq']
                        df.loc[i,"score"]=dct_score['seg_eq'] 
 
            except:
                print(sys.exc_info())
                #print(i,r,pre_seg,col_segdiff,post_seg)
                #print(f"current segment{curr_seg},values{value_compare}")
                sys.exit(0)


for i,r in df.iterrows():
    if r['eval']==dct_catg['last_eq']:
        df.loc[i,"score"]=get_score(df,i[0])



print("writing verdict ...")

for gname in df.index.get_level_values('GeneID'):

    if df.at[(gname,1),'is_only']==1:
        df.at[(gname,1),"verdict"]=df.at[(gname,1),"eval"]
        continue
    
    if df.at[(gname,1),'is_first']==1:  
        z=get_sum_verdict(df,gname)
        
        if z==0:
            df.at[(gname,1),"verdict"]="no APA, no differences"
        
        if z==-4:
            df.at[(gname,1),"verdict"]="equal CDS, longer 3’UTR in %s" %lb0
        if z==4:
            df.at[(gname,1),"verdict"]="equal CDS, longer 3’UTR in %s" %lb1
        
        if z==-3:
            df.at[(gname,1),"verdict"]="no APA, more abundant in %s" %lb0
        if z==3:
            df.at[(gname,1),"verdict"]="no APA, more abundant in %s" %lb1
        
        if z==-1:
            df.at[(gname,1),"verdict"]="less CDS but longer 3’UTR in %s" %lb0
        if z==1:
            df.at[(gname,1),"verdict"]="less CDS but longer 3’UTR in %s" %lb1

        if z==-7:
            df.at[(gname,1),"verdict"]="more CDS and longer 3’UTR in %s" %lb0
        if z==7:
            df.at[(gname,1),"verdict"]="more CDS and longer 3’UTR in %s" %lb1




for gname in df.index.get_level_values('GeneID'):
    vv=df.loc[(gname,),]
    vv=vv.sort_index(level='segment')
    if vv.iloc[0]['is_first']==1:
        df.at[(gname,1),'lastdiff']=vv.iloc[-1]['%s-%snormdiff' %(lb1,lb0)]
    else:
        df.at[(gname,1),'lastdiff']=0


for e in g0+g1:
    df.rename(columns={e:'TPM%s' %e})
    df_filtered.rename(columns={e:'TPM%s' %e})


df=df.append(df_filtered,sort=False)

tt=df.pop('FILTER')
df['FILTER']=tt

print("writing the results...")
df.to_csv('%s_%s_vs_%s.csv' %(outfile,lb1,lb0))


#print("writing the main table...")
#df.to_csv('%s_%s_vs_%s.csv' %(outfile2,lb1,lb0))





