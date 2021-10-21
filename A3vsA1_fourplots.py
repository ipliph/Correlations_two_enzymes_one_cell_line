import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

lim_x = 150
lim_y = 150

step = 2.5
###Data extraction part

df_ic50 = pd.read_excel('PDI inhibitors data.xlsx', sheet_name='ic50',index_col=0)
df_ic50 = df_ic50[list(df_ic50.columns)[6:12]]
df_ic50 = df_ic50.rename(columns={'Inhibiting activity':'PDI A1', 'Unnamed: 8':'PDI A3', 'Inhibiting activity, IC50 (µM)':'HT_1080', 'Unnamed: 10':'CaCo_2', 'Unnamed: 11':'MDA_MB', 'Unnamed: 12':'MCF_7'})
df_ic50 = df_ic50.iloc[2:]

df_ic50 = df_ic50.replace(['>200','n.t.'],[200,0])
###df[['First','Last']] = df.Name.str.split(expand=True)
###'A3HT_1080','CaCo_2','MDA_MB','MCF_7'

fig, ax = plt.subplots(2, 3, gridspec_kw={'width_ratios': [5, 5, 0.5]},figsize=(15,10))



#'PDI_A3HT_1080'
df_ic50[['HT-1080av','HT-1080std']] = df_ic50.HT_1080.str.split('±',expand=True)

for num in range(82):
	x = df_ic50['PDI A1'][num]
	y = df_ic50['PDI A3'][num]

	if float(df_ic50['HT-1080av'][num]) < 2*step:
		alpha = 1
	elif float(df_ic50['HT-1080av'][num]) < 3*step:
		alpha = 0.9
	elif float(df_ic50['HT-1080av'][num]) < 4*step:
		alpha = 0.8
	elif float(df_ic50['HT-1080av'][num]) < 5*step:
		alpha = 0.7
	elif float(df_ic50['HT-1080av'][num]) < 6*step:
		alpha = 0.6
	elif float(df_ic50['HT-1080av'][num]) < 7*step:
		alpha = 0.5
	elif float(df_ic50['HT-1080av'][num]) < 8*step:
		alpha = 0.4
	elif float(df_ic50['HT-1080av'][num]) < 9*step:
		alpha = 0.3
	elif float(df_ic50['HT-1080av'][num]) < 10*step:
		alpha = 0.2
	else:
		alpha = 0.1
	
	if x > lim_x or y > lim_y:
		pass
	else:
		ax[0,0].scatter(x,y,color='b',alpha=alpha)

ax[0,0].set_title("HT-1080")
ax[0,0].set(ylabel='PDI A3 - IC50 [µM]')

#'CaCo_2'
df_ic50[['CaCo_2av','CaCo_2std']] = df_ic50.CaCo_2.str.split('±',expand=True)

for num in range(82):
	x = df_ic50['PDI A1'][num]
	y = df_ic50['PDI A3'][num]
	alpha = 0

	if float(df_ic50['CaCo_2av'][num]) < 2*step:
		alpha = 1
	elif float(df_ic50['CaCo_2av'][num]) < 3*step:
		alpha = 0.9
	elif float(df_ic50['CaCo_2av'][num]) < 4*step:
		alpha = 0.8
	elif float(df_ic50['CaCo_2av'][num]) < 5*step:
		alpha = 0.7
	elif float(df_ic50['CaCo_2av'][num]) < 6*step:
		alpha = 0.6
	elif float(df_ic50['CaCo_2av'][num]) < 7*step:
		alpha = 0.5
	elif float(df_ic50['CaCo_2av'][num]) < 8*step:
		alpha = 0.4
	elif float(df_ic50['CaCo_2av'][num]) < 9*step:
		alpha = 0.3
	elif float(df_ic50['CaCo_2av'][num]) < 10*step:
		alpha = 0.2
	else:
		alpha = 0.1

	if x > lim_x or y > lim_y:
		pass
	else:
		ax[0,1].scatter(x,y,color='g',alpha=alpha)
ax[0,1].set_title("CaCo_2")



#'MDA_MB'
df_ic50[['MDA_MBav','MDA_MBstd']] = df_ic50.MDA_MB.str.split('±',expand=True)

for num in range(82):
	x = df_ic50['PDI A1'][num]
	y = df_ic50['PDI A3'][num]
	alpha = 0

	if float(df_ic50['MDA_MBav'][num]) < 2*step:
		alpha = 1
	elif float(df_ic50['MDA_MBav'][num]) < 3*step:
		alpha = 0.9
	elif float(df_ic50['MDA_MBav'][num]) < 4*step:
		alpha = 0.8
	elif float(df_ic50['MDA_MBav'][num]) < 5*step:
		alpha = 0.7
	elif float(df_ic50['MDA_MBav'][num]) < 6*step:
		alpha = 0.6
	elif float(df_ic50['MDA_MBav'][num]) < 7*step:
		alpha = 0.5
	elif float(df_ic50['MDA_MBav'][num]) < 8*step:
		alpha = 0.4
	elif float(df_ic50['MDA_MBav'][num]) < 9*step:
		alpha = 0.3
	elif float(df_ic50['MDA_MBav'][num]) < 10*step:
		alpha = 0.2
	else:
		alpha = 0.1

	if x > lim_x or y > lim_y:
		pass
	else:
		ax[1,0].scatter(x,y,color='r',alpha=alpha)

ax[1,0].set_title("MDA-MB")
ax[1,0].set(xlabel='PDI A1 - IC50 [µM]',ylabel='PDI A3 - IC50 [µM]')




#'MCF_7'
df_ic50[['MCF_7av','MCF_7std']] = df_ic50.MCF_7.str.split('±',expand=True)

for num in range(82):
	x = df_ic50['PDI A1'][num]
	y = df_ic50['PDI A3'][num]
	alpha = 0

	if float(df_ic50['MCF_7av'][num]) < 2*step:
		alpha = 1
	elif float(df_ic50['MCF_7av'][num]) < 3*step:
		alpha = 0.9
	elif float(df_ic50['MCF_7av'][num]) < 4*step:
		alpha = 0.8
	elif float(df_ic50['MCF_7av'][num]) < 5*step:
		alpha = 0.7
	elif float(df_ic50['MCF_7av'][num]) < 6*step:
		alpha = 0.6
	elif float(df_ic50['MCF_7av'][num]) < 7*step:
		alpha = 0.5
	elif float(df_ic50['MCF_7av'][num]) < 8*step:
		alpha = 0.4
	elif float(df_ic50['MCF_7av'][num]) < 9*step:
		alpha = 0.3
	elif float(df_ic50['MCF_7av'][num]) < 10*step:
		alpha = 0.2
	else:
		alpha = 0.1

	if x > lim_x or y > lim_y:
		pass
	else:
		ax[1,1].scatter(x,y,color='brown',alpha=alpha)


ax[1,1].set_title("MCF-7")
ax[1,1].set(xlabel='PDI A1 - IC50 [µM]')

ax[0,2].set_visible(False)

x = [0 for num in range(10)]
y = [num*step for num in range(1,11)]
al = np.linspace(0.1,1,10)
al = np.flip(al)
for num in range(10):
	ax[1,2].scatter(x[num], y[num], c='b',alpha=al[num])
ax[1,2].set_ylabel('Antiproliferative effect against cancer cell line\nIC50 [µM]')
ax[1,2].set_title("Legend")
plt.xticks([])

fig.suptitle('Novel PDI inhibitors'+
	' - relationship between A1 and A3 subtypes IC50 values\n'+
	'combined with antiproliferative effect\n'+
	'(points with value IC50 > 150 - excluded)')

plt.savefig('Combined 4 plots_150.png')
plt.show()
