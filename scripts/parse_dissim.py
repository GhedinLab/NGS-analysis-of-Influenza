import pandas as pd


filelist = ['rownames.H1N1.HA.kate.dissim_all.key.csv',
'rownames.H3N2.HA.kate.dissim_all.key.csv',
'rownames.H3N2.H3N2-HA.dissim_all.key.csv',
'rownames.H1N1.H1N1-HA.dissim_all.key.csv']

outfile = open('H1.H3.redemux.reseq.plotme.txt','w')
print>>outfile,'value,group,strain,type'

for filedissim in filelist:
    if 'kate' in filedissim:
        typeseq = 'redemux'
    else:
        typeseq = 'reseq'
    if 'H1N1' in filedissim:
        strain = 'H1N1'
    elif 'H3N2' in filedissim:
        strain = 'H3N2'
    df = pd.read_csv(filedissim,index_col=0)

    householdlist = []
    otherlist = []
    for i,row in df.iterrows():
        print row.name
        source_household = row.name.split('-')[0]
        # print row.index
        for sample in row.index:
            target_household = sample.split('-')[0]
            if target_household == source_household:
                householdlist.append(row[sample])
            else:
                otherlist.append(row[sample])
        #     print sample.name
        #     print sample
        # print row
    for i in householdlist:
        print>>outfile,str(i)+','+'household,'+strain+','+typeseq

    for i in otherlist:
        print>>outfile,str(i)+','+'other,'+strain+','+typeseq

