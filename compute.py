def compute(r,cat_dict):
    import pandas as pd
    import numpy as np
    X=pd.read_csv('flaskexample/static/peplife_X.csv')
    X=X.drop('Unnamed: 0',axis=1)
    y=pd.read_csv('flaskexample/static/peplife_y.csv',header=None)
    y=y.drop(0,axis=1)
    y=np.ravel(y)

    from sklearn.preprocessing import StandardScaler
    scaler=StandardScaler()
    X_f=pd.DataFrame(scaler.fit_transform(X.iloc[:,5:]))
    X_c=pd.DataFrame(X.iloc[:,1:5])
    X=pd.concat([X_c,X_f],axis=1)
    X=pd.get_dummies(X)

    def generate_features(seq):
        """
        expect a list of sequences (a list of one for single sequence input)
        return pandas dataframe containing 20 unscaled features 10 from modlamp, 
        10 from custom feature generateion

        """
        from modlamp.descriptors import GlobalDescriptor
        custom_features = pd.Series(seq).apply(generate_custom_features)
        gdesc = GlobalDescriptor(seq)
        gdesc.calculate_all()
        modlamp_features = pd.DataFrame(gdesc.descriptor)
        modlamp_features.columns=gdesc.featurenames
        out = pd.concat([modlamp_features,custom_features],axis=1)
        return out
    def generate_custom_features(seq):
        def Asp_exists(seq):
            return seq.count('D')
        def Ser_exists(seq):
            return seq.count('S')
        def Cys_exists(seq):
            return seq.count('C')
        def Met_exists(seq):
            return seq.count('M')
        def DP_exists(seq):
            return seq.count('DP')
        def DG_exists(seq):
            return seq.count('DG')
        def NG_exists(seq):
            return seq.count('NG')
        def QG_exists(seq):
            return seq.count('QG')
        def Diketo_exists1(seq):
            return int(seq[2]=='G')
        def Diketo_exists2(seq):
            tmp=('P' in seq[:2]) or ('G' in seq[:2])
            return int(tmp and Diketo_exists1(seq))
        out=pd.Series([Asp_exists(seq),Ser_exists(seq),Cys_exists(seq),Met_exists(seq),
        DP_exists(seq),DG_exists(seq),NG_exists(seq),QG_exists(seq),Diketo_exists1(seq),
        Diketo_exists2(seq)],index=['Asp','Ser','Cys','Met','DP','DG','NG','QG','Diketo','Diketo_2'])
        return out
    def generate_cat(cat_dict):
        out=[]
        if cat_dict['in_vivo_in_vitro']=='in vitro':
            out+=[0,1]
        elif cat_dict['in_vivo_in_vitro']=='in vivo':
            out+=[0,1]
        if cat_dict['Linear_cyclic']=='Cyclic':
            out+=[1,0]
        elif cat_dict['Linear_cyclic']=='Linear':
            out+=[0,1]
        if cat_dict['N_ter_mod']=='Acetylation':
            out+=[1,0,0,0]
        elif cat_dict['N_ter_mod']=='Free':
            out+=[0,1,0,0]
        elif cat_dict['N_ter_mod']=='Glycosylation':
            out+=[0,0,1,0]
        elif cat_dict['N_ter_mod']=='Hydroxylation':
            out+=[0,0,0,1]
        if cat_dict['C_ter_mod']=='Amidation':
            out+=[1,0,0,0]
        elif cat_dict['C_ter_mod']=='Free':
            out+=[0,1,0,0]
        elif cat_dict['C_ter_mod']=='Pegylation':
            out+=[0,0,1,0]
        elif cat_dict['C_ter_mod']=='Propylamidation':
            out+=[0,0,0,1]
        return pd.Series(out)

    def generate_instance(seq,cat_dict):
        """
        Expect a list of sequences and a catorical dictionary, the categorical varibles should have four keys,
        they are 'in_vivo_in_vitro','Linear_cyclic','N_ter_mod','C_ter_mod'
        Return a dataframe with all continous varibles and categorical varibles generated in the same order of 
        training data.
        """
        df1=pd.DataFrame(scaler.transform(generate_features(seq)))
        df2=pd.concat([generate_cat(cat_dict)]*len(seq),axis=1)
        df2=df2.transpose()
        out=pd.concat([df1,df2],axis=1)
        out.columns=X.columns
        return out
    
    def generate_degenerate(seq):
        '''
        Expect individual sequence as string,allow for code degeneration
        '''
        seq=seq.upper()
        AA='ARNDCQEGHILKMFPSTWYVX[], '
        for i in seq:
            if i not in AA:
                raise Exception('Unrecognized Amino Acid code')
        seq=seq.upper()
        import re
        m=re.search('\[(.+)\]',seq)
        if not m:
            return [seq]
        elif m.group(1)=='X':
            seq_list=[]
            found=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
            seq1=seq[:seq.index('[')]
            seq2=seq[(seq.index(']')+1):]
            for i in found:
                seq_list.append(seq1+i+seq2)
            return seq_list
        else:
            seq_list=[]
            found=m.group(1).split(',')
            seq1=seq[:seq.index('[')]
            seq2=seq[(seq.index(']')+1):]
            for i in found:
                seq_list.append(seq1+i+seq2)
            return seq_list

    def generate_expanded_seq(seq_list):
        new_list=[]
        for seq in seq_list:
            new_list+=generate_degenerate(seq)
        return new_list

    def read_seq_string(string):
        """
        Expect a string of sequence return a list of sequences.
        """
        return [seq.strip() for seq in filter(None,string.split(';'))]

    def read_seq(string,cat_dict):
        """
        Expect a semicolumn seperated string of sequences,and a dictionary specifying sequence modifications
        convert to feature dataframe
        The sequence can use '[]' to create degenerate sequences
        return both extended sequence list and the feature dataframe
        """
        tmp1=read_seq_string(string)
        seq_list=generate_expanded_seq(tmp1)
        df=generate_instance(seq_list,cat_dict)
        return seq_list,df
    seq_list, X_test =read_seq(r,cat_dict)
    if 'clf' not in globals():
        global clf
        from sklearn.svm import SVC
        clf=SVC(kernel='linear')
    try:
        y_test=clf.predict(X_test)
    except:
        clf.fit(X,y)
        y_test=clf.predict(X_test)
    out_data=pd.DataFrame([seq_list,y_test.tolist()]).transpose()
    out_data.rename(columns={0:'Peptide Sequence',1:'Predicted Half Life'},inplace=True)
    def cleanup(df):
        def apply_color(val):
            if val== '>60min':
                color='#7BCCB5'
            elif val=='5~60min':
                color='#F5F5DC'
            elif val =='<5min':
                color = '#CD7F32'
            else:
                color='#d5d8dc'
            return 'color : %s'% color 
        df['Predicted Half Life']=pd.Categorical(df['Predicted Half Life'],['<5min','5~60min','>60min'])
        df_out =  df.sort_values('Predicted Half Life',ascending=False).reset_index(drop=True)
        s=df_out.style.applymap(apply_color)
        return s.render()
    return cleanup(out_data)
