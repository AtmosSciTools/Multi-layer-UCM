####
####
import pandas as pd




df = pd.read_csv('Kugahara-full-datasets-ver1.1.csv', header=1)


'''
w = df.iloc[:24,2:13]
s1 = df.iloc[25:49,2:13]
s2 = df.iloc[50:74,2:13]

wd = pd.to_datetime('2005')+ pd.Timedelta(str(df.iloc[0,1]-1)+' days')
sd1 = pd.to_datetime('2005')+ pd.Timedelta(str(df.iloc[25,1]-1)+' days')
sd2 = pd.to_datetime('2005')+ pd.Timedelta(str(df.iloc[50,1]-1)+' days')


odir = "Kugahara-full-datasets"
for d, date in zip([w,s1,s2],[wd,sd1,sd2]):
    print wd, sd1, sd2
    d["date"] = [ date + pd.Timedelta(str(h)+' hours') for h in d[u'TIME(JST)'] ]
    d.set_index("date", inplace=True)
    d.drop([u'TIME(JST)'],axis=1,inplace=True)
    d.to_csv(date.strftime(odir+'/ts_top_allvar_%Y%m%d.csv'))
'''


w = df.iloc[:24,2:]
s1 = df.iloc[25:49,2:]
s2 = df.iloc[50:74,2:]

wd = pd.to_datetime('2005')+ pd.Timedelta(str(df.iloc[0,1]-1)+' days')
sd1 = pd.to_datetime('2005')+ pd.Timedelta(str(df.iloc[25,1]-1)+' days')
sd2 = pd.to_datetime('2005')+ pd.Timedelta(str(df.iloc[50,1]-1)+' days')



odir = "Kugahara-full-datasets"
for d, date in zip([w,s1,s2],[wd,sd1,sd2]):
    print wd, sd1, sd2
    time = [ date + pd.Timedelta(str(h)+' hours') for h in d[u'TIME(JST)'] ]
    d["date"] = time
    d.set_index("date", inplace=True)
    d.drop([u'TIME(JST)'],axis=1,inplace=True)
    #d.iloc[:,:10].to_csv(date.strftime(odir+'/ts_top_allvar_%Y%m%d.csv'))


    hum  = d.iloc[:,-11:-1]
    height = [float(h.split('m')[0][3:]) for h in hum.columns]

    temp = d.iloc[:,-23:-12]
    height = [float(h.split('m')[0][4:]) for h in temp.columns]

    wind = d.iloc[:,11:15]
    hum.to_csv(date.strftime(odir+'/hum_profile_%Y%m%d.csv'))
    temp.to_csv(date.strftime(odir+'/tmp_profile_%Y%m%d.csv'))
    wind.to_csv(date.strftime(odir+'/wps_profile_%Y%m%d.csv'))





