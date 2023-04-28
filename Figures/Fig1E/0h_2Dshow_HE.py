import matplotlib.pyplot as plt
import matplotlib.image as imgplt



##################################

spot = "spot1"
imgnum = "1"
slide = "A1.txt"
crop_size = 30
path = "./data/"

################raw spot#########################

data = []

with open(path+"prepare/"+spot+".txt") as f:
    for line in f.readlines():
        temp = line.split()
        temp[0] = int(temp[0])
        temp[1] = int(temp[1])
        data.append(temp)

data.sort(key=lambda x: (x[1], x[0]))

all_s = []
sess = []
s = int(data[0][1])

for i in range(len(data)):
    if data[i][1] == s:
        sess.append(data[i])
    else:
        if int(data[i][1]) < s + crop_size:
            sess.append(data[i])
        else:
            maxv = int(data[i-1][1])
            averv = (s + maxv)/2
            for j in range(len(sess)):
                sess[j][1] = averv
            all_s.append(sess)
            sess = []
            s = data[i][1]
            sess.append(data[i])

maxv = int(data[-1][1])
averv = (s + maxv)/2
for j in range(len(sess)):
    sess[j][1] = averv
all_s.append(sess)

alls = []
for i in range(len(all_s)):
    for j in range(len(all_s[i])):
        alls.append(all_s[i][j])

##

Data = {}
for each in alls:
    temp = each[2].split("_")[0]
    Data[temp] = [each[0], each[1]]

lab = []
with open(path + "cell2location/" + slide) as f:
    for line in f.readlines():
        temp = line.split()
        lab.append(temp)

selectall = []
for each in lab:
    if each[0].split("-c")[0] in Data.keys():
        temp = each
        temp.extend(Data[each[0].split("-c")[0]])
        selectall.append(temp)


## Plot

fig, ax = plt.subplots()

# The range of x and y
ax.set_xlim(left=0, right=8000)
ax.set_ylim(bottom=6000, top=0) 
plt.axis('off')


Image = imgplt.imread(path+"prepare/image-"+imgnum+".png")
plt.imshow(Image)

for each in selectall:
    value = [float(k) for k in each[1:-2]]
    celltype = value.index(max(value))+1
    if int(celltype) == 1:  # Epidermal
        c = 'blue'
    if int(celltype) == 2:  # Gut
        c = 'green'
    if int(celltype) == 3:  # Muscle
        c = 'red'
    if int(celltype) == 4:  # Neoblast
        c = 'grey'
    if int(celltype) == 5:  # Neuronal
        c = 'orange'
    if int(celltype) == 6:  # Parenchymal
        c = 'pink'
    if int(celltype) == 7:  # Secretory
        c = 'purple'
    plt.scatter(float(each[-2]), float(each[-1]), c=c, s=2, marker='o')

plt.savefig("0h_"+spot+"_rawdata_HE.pdf", dpi=300)


################merge spot#########################

fill = []

with open(path+"predict/predict_"+spot+".txt") as f:
    for line in f.readlines():
        temp = line.split()
        temp[3] = float(temp[3])
        temp[4] = float(temp[4])
        fill.append(temp)

Colors = {"Neuronal": 'orange', "Neoblast": 'grey', "Parenchymal": 'pink', "Gut": 'green', "Muscle": 'red', "Secretory": 'purple', "Epidermal": 'blue'}

for each in fill:
    c = Colors[str(each[2])]
    plt.scatter(float(each[3]), float(each[4]), c=c, s=2, marker='o')

plt.savefig("0h_"+spot+"_mergedata_HE.pdf", dpi=300)
#plt.clf()

# ## Plot
#
# fig, ax = plt.subplots()
#
# ax.set_xlim(left=0, right=8000)
# ax.set_ylim(bottom=6000, top=0)
# plt.axis('off')
#
#
# Image = imgplt.imread(path+"prepare/image-"+imgnum+".png")
# plt.imshow(Image)
#
# for each in fill:
#     c = Colors[str(each[2])]
#     plt.scatter(float(each[3]), float(each[4]), c=c, s=2, marker='o')
#
# plt.savefig("0h_"+spot+"_filldata_HE.pdf", dpi=300)