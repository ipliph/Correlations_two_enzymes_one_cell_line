import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import math

lim_A1 = 200
lim_A3 = 200
lim_canc = 500

step = 2.5

cancer_cell_line_no = 4
###Data extraction part

df_ic50 = pd.read_excel('PDI inhibitors data.xlsx', sheet_name='ic50',index_col=0)
df_ic50 = df_ic50[list(df_ic50.columns)[0:12]]
df_ic50 = df_ic50.rename(columns={
	'Unnamed: 1':'Scaffold',
	'Substituent':'R1',
	'Unnamed: 3':'R2',
	'Unnamed: 4':'R3',
	'Unnamed: 5':'R4',
	'Unnamed: 6':'R5',
	'Inhibiting activity':'PDI A1', 
	'Unnamed: 8':'PDI A3', 
	'Inhibiting activity, IC50 (µM)':'HT_1080',
	 'Unnamed: 10':'CaCo_2', 
	 'Unnamed: 11':'MDA_MB', 
	 'Unnamed: 12':'MCF_7'})

df_ic50 = df_ic50.iloc[2:]
df_ic50 = df_ic50.replace(['>200','n.t.'],[200,0])

df_ic50[['HT_1080av','HT_1080std']] = df_ic50.HT_1080.str.split('±',expand=True)
df_ic50[['CaCo_2av','CaCo_2std']] = df_ic50.CaCo_2.str.split('±',expand=True)
df_ic50[['MDA_MBav','MDA_MBstd']] = df_ic50.MDA_MB.str.split('±',expand=True)
df_ic50[['MCF_7av','MCF_7std']] = df_ic50.MCF_7.str.split('±',expand=True)

data_set_dict = {1:(df_ic50['HT_1080av'],'HT-1080'),2:(df_ic50['CaCo_2av'],'CaCo-2'),3:(df_ic50['MDA_MBav'],'MDA-MB'),4:(df_ic50['MCF_7av'],'MCF-7')}
data_set_to_process = data_set_dict[cancer_cell_line_no][0]

#Eliminating NaN
df_ic50_subtNaN = []
for num in range(len(data_set_to_process)):
	if math.isnan(float(data_set_to_process[num])):
		df_ic50_subtNaN.append(0)
	else:
		df_ic50_subtNaN.append(float(data_set_to_process[num]))

d_short = {'Comp_name':df_ic50.index,'Scaffold':df_ic50['Scaffold'],'R3':df_ic50['R3'],'R4':df_ic50['R4'],'R5':df_ic50['R5'],'PDI_A1':df_ic50['PDI A1'],'PDI_A3':df_ic50['PDI A3'],'IC50':df_ic50_subtNaN}
df_short = pd.DataFrame(d_short)

df_short_sorted = df_short.sort_values(by='IC50',ignore_index=True)



drop_list = []

for num in range(len(df_short_sorted['IC50'])):
	if df_short_sorted['IC50'][num] == 0.0:
		drop_list.append(num)

for num in range(len(df_short_sorted['IC50'])):
	if df_short_sorted['PDI_A1'][num] > lim_A1:
		drop_list.append(num)

for num in range(len(df_short_sorted['IC50'])):
	if df_short_sorted['PDI_A3'][num] > lim_A3:
		drop_list.append(num)

for num in range(len(df_short_sorted['IC50'])):
	if df_short_sorted['Scaffold'][num] == 'S1' and (df_short_sorted['R5'][num] == 'COOMe'):
		#pass
		drop_list.append(num)

for num in range(len(df_short_sorted['IC50'])):
	if df_short_sorted['Scaffold'][num] == 'S2' and (df_short_sorted['R4'][num] == 'COOMe'):
		#pass
		drop_list.append(num)

for num in range(len(df_short_sorted['IC50'])):
	if df_short_sorted['Scaffold'][num] == 'S3' and (df_short_sorted['R4'][num] == 'COOMe'):
		#pass
		drop_list.append(num)



df_droped = df_short_sorted.drop(index=drop_list)
df_droped = df_droped.reset_index()

###Corelation

##PDI A1
xA1 = df_droped['PDI_A1']
xA3 = df_droped['PDI_A3']

A1= np.vstack([xA1, np.ones(len(xA1))]).T
a_A1, b_A1 = np.linalg.lstsq(A1, df_droped['IC50'], rcond=None)[0]

cor_matrix = np.corrcoef(xA1, df_droped['IC50'])
cor_xy = cor_matrix[0,1]
r_squared_A1 = cor_xy**2

y_eq_A1 = a_A1 * df_droped['PDI_A1'] + b_A1

##PDI A3
A3 = np.vstack([xA3, np.ones(len(xA1))]).T
a_A3, b_A3 = np.linalg.lstsq(A3, df_droped['IC50'], rcond=None)[0]

cor_matrix = np.corrcoef(xA3, df_droped['IC50'])
cor_xy = cor_matrix[0,1]
r_squared_A3 = cor_xy**2

y_eq_A3 = a_A3 * df_droped['PDI_A3'] + b_A3


###Matplotlib
fig = plt.figure(figsize=(10,8))
gs = fig.add_gridspec(1, 2, hspace=0, wspace=0.01)
ax1, ax2 = gs.subplots(sharey='row',sharex='col')
suptitle = ('Correlation between antiproliferative effect against '+str(data_set_dict[cancer_cell_line_no][1])+
	'\nand inhibiting activity against PDI subclasses')
fig.suptitle(suptitle)

#Plot 1
for num in list(df_droped.index):
	if df_droped['Scaffold'][num] == 'S1':
		ax1.scatter(df_droped['PDI_A1'][num],df_droped['IC50'][num],c='b')
	elif df_droped['Scaffold'][num] == 'S2':
		ax1.scatter(df_droped['PDI_A1'][num],df_droped['IC50'][num],c='r')
	elif df_droped['Scaffold'][num] == 'S3':
		ax1.scatter(df_droped['PDI_A1'][num],df_droped['IC50'][num],c='g')

ax1.plot(xA1,y_eq_A1,c='r')
xlabel_A1 = 'Inhibiting activity against PDI A1, IC50 (µM)\na='+str(round(a_A1,5))+',b='+str(round(b_A1,5))+',r_sq='+str(round(r_squared_A1,4))
ax1.set(xlabel=xlabel_A1, ylabel = "Antiproliferative effect against cancer cell line, IC50 (µM)")

#Plot 2
for num in list(df_droped.index):
	if df_droped['Scaffold'][num] == 'S1':
		b = ax2.scatter(df_droped['PDI_A3'][num],df_droped['IC50'][num],c='b')
	elif df_droped['Scaffold'][num] == 'S2':
		r = ax2.scatter(df_droped['PDI_A3'][num],df_droped['IC50'][num],c='r')
	elif df_droped['Scaffold'][num] == 'S3':
		g = ax2.scatter(df_droped['PDI_A3'][num],df_droped['IC50'][num],c='g')

ax2.plot(xA3,y_eq_A3,c='r')
xlabel_A3 = 'Inhibiting activity against PDI A3, IC50 (µM)\na='+str(round(a_A3,5))+',b='+str(round(b_A3,5))+',r_sq='+str(round(r_squared_A3,4))
ax2.set(xlabel=xlabel_A3)

ax2.legend([b,r,g],['Scaffold 1','Scaffold 2','Scaffold 3'],loc=1)

file_name = data_set_dict[cancer_cell_line_no][1]+' amids_plot.png'
plt.savefig(file_name)

#plt.show()





















