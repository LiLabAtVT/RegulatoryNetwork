import numpy as np
import csv
import scipy.stats
import os

<<<<<<< HEAD
os.chdir("/Users/Craftsman/Documents/Wei He/Blacksburg/Study/research/")

=======
>>>>>>> 3aa85e2a277dbd026babddeb57c6d2f8c89b5e9a
filename = 'GSE5623_ave_salt_root.csv'
data_mat = []
with open(filename, 'rt') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=' ')
    for row in spamreader:
        data_mat.append(row)
csvfile.close()

time_step = data_mat[0]
time_step = time_step[0]
time_step = time_step.strip()
time_step = time_step.split(',')
del time_step[0]
np.array(time_step)

time_stepuse = list()
for time1 in time_step:
    ind1 = time1.index('m')
    timeuse = float(time1[0:ind1])
    time_stepuse.append(timeuse)
np.asarray(time_stepuse)

data_mat_allgenename = []
data_mat_data = []
for i in range(1, len(data_mat)):
        data_mati = data_mat[i]
        data_mati = data_mati[0]
        data_mati = data_mati.strip()
        data_mati = data_mati.split(',')
        genename1 = data_mati[0]
        del data_mati[0]
        data_mati = [float(ii) for ii in data_mati]
        data_mat_allgenename.insert(i - 1, genename1)
        data_mat_data.append(data_mati)

mRNAtime = 30 # set the mRNA time to 30 mins
a_value = mRNAtime/np.log(2)
dataoutput = list()
for data1 in data_mat_data:
    np.array(data1)
    data1diff = np.diff(data1)
    timediff = np.diff(time_stepuse)
    result1 = np.divide(data1diff,timediff)
    result1 = 1/a_value * result1
    result1 = result1 + data1[0:-1]
    dataoutput.append(result1)

pathsave = 'processed'
if not os.path.exists(pathsave):
    os.makedirs(pathsave)  # create folder
os.chdir(pathsave)


exp_num = len(data1diff)
exp_name = ['  ']
for i in range(0, exp_num):
    exp_name.append('exp'+str(i+1))

with open(filename[0: -4] + 'process.csv', 'w') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_NONE)
    wr.writerow(exp_name)
    for i in range(0, len(dataoutput)):
        rel1 = dataoutput[i]
        writeresult = [data_mat_allgenename[i]]
        for j in range(0, len(rel1)):
            writeresult.append(str(rel1[j]))
        wr.writerow(writeresult)
myfile.close()












