import matplotlib.pyplot as plt
import numpy as np
np.set_printoptions(suppress=True)
import copy
import time
import math



############################

celltype = "Secretory"

##
inputpath_raw = "./data/cell2location/"
inputpath_rawxy = "./data/prepare/"
inputpath_fill = "./data/predict/"
inputpath_rest = "./data/rest//"

##
lab = {"Epidermal": [0, "Blues"], "Gut": [1, "Greens"], "Muscle": [2, "Reds"], "Neoblast": [3, "binary"],
       "Neuronal": [4, "Oranges"],  "Parenchymal": [5, "RdPu"], "Secretory": [6, "Purples"], }

##
crop_size = 30

##
filter_value = 0

##
spot_NUM = [1, 2, 3, 5, 6, 7, 8, 9, 10, 11]

##
Spot_Slide = {"spot1": "A1", "spot2": "B1", "spot3": "C1", "spot5": "A1", "spot6": "B1",
              "spot7": "C1", "spot8": "D1", "spot9": "A1", "spot10": "B1", "spot11": "C1"}

##
resetx, resety = 1151, 1019
imgx, imgy = 5913, 2039

Transmessage = {"spot2": [1222, 972, 7087, 3146, "N"], "spot3": [1183, 1074, 7150, 2948, "S"],
                "spot4": [1183, 990, 7150, 2864, "S"], "spot5": [1246, 1064, 7207, 2958, "S"],
                "spot6": [1367, 900, 7217, 3113, "N"], "spot7": [1050, 1061, 7007, 2966, "S"],
                "spot8": [1306, 1028, 7267, 2922, "S"], "spot9": [1034, 977, 6996, 2869, "S"],
                "spot10": [1136, 922, 7032, 3010, "N"], "spot11": [1233, 1049, 7163, 3037, "N"]}

##
img_spot = {"img7": "spot1", "img9": "spot2", "img12": "spot3", "img14": "spot4", "img17": "spot5",
             "img20": "spot6", "img22": "spot7", "img24": "spot8", "img27": "spot9", "img30": "spot10"}


############################

##
def Srotate(angle, valuex, valuey, pointx, pointy):

    sRotatex = (valuex-pointx)*math.cos(angle) + (valuey-pointy)*math.sin(angle) + pointx
    sRotatey = (valuey-pointy)*math.cos(angle) - (valuex-pointx)*math.sin(angle) + pointy

    return sRotatex, sRotatey


##
def Nrotate(angle, valuex, valuey, pointx, pointy):

    nRotatex = (valuex-pointx)*math.cos(angle) - (valuey-pointy)*math.sin(angle) + pointx
    nRotatey = (valuex-pointx)*math.sin(angle) + (valuey-pointy)*math.cos(angle) + pointy

    return nRotatex, nRotatey


##
def get_rawspot(fname, fslide, rawindex, spotnum):

    data = []

    with open(fname, 'r') as f:
        for line in f.readlines():
            temp = line.split()
            temp[0] = float(temp[0])
            temp[1] = float(temp[1])
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
                maxv = int(data[i - 1][1])
                averv = (s + maxv) / 2
                for j in range(len(sess)):
                    sess[j][1] = averv
                all_s.append(sess)
                sess = []
                s = data[i][1]
                sess.append(data[i])
    maxv = int(data[-1][1])
    averv = (s + maxv) / 2
    for j in range(len(sess)):
        sess[j][1] = averv
    all_s.append(sess)

    alls = []
    for i in range(len(all_s)):
        for j in range(len(all_s[i])):
            alls.append(all_s[i][j])

    Data = {}
    for each in alls:
        temp = each[2].split("_")[0]
        Data[temp] = [each[0], each[1]]

    ## cell2location

    slide_C = []
    with open(fslide, 'r') as f:
        for line in f.readlines():
            temp = line.split()
            slide_C.append(temp)

    spot_C = []
    for each in slide_C:
        if each[0].split("-c")[0] in Data.keys():
            temp = each[1:]
            temp.extend(Data[each[0].split("-c")[0]])
            spot_C.append(temp)

    ##

    m1 = np.array(spot_C, dtype="float")
    m2 = m1.reshape(len(spot_C), 9)
    m3 = []
    for each in m2:
        total_value = [float(k) for k in each[0:-2]]
        if each[rawindex]/sum(total_value) >= filter_value:
            m3.append(each[-2])
            m3.append(each[-1])
            m3.append(float(each[rawindex])/sum(total_value))
    m4 = np.array(m3, dtype='float')
    m5 = m4.reshape(int(len(m4)/3), 3)
    # x = [k[0] for k in m5]
    # y = [k[1] for k in m5]
    # alpha = [k[2] for k in m5]

    ##

    if spotnum == "spot1":
        finaldata = m5.tolist()

    else:
        xleft, yleft, xright, yright, direction = Transmessage[spotnum]

        x3 = float(xleft) - resetx  
        y3 = float(yleft) - resety
        basex1 = 0
        basey1 = 0
        x4 = float(xright) - resetx
        y4 = float(yright) - resety
        basex2 = float(imgx)  
        basey2 = float(imgy) 

        x_shift = x3 - basex1
        y_shift = y3 - basey1
        a = math.sqrt((x4 - (basex2 + x_shift)) ** 2 + (y4 - (basey2 + y_shift)) ** 2)
        b = math.sqrt((x4 - x3) ** 2 + (y4 - y3) ** 2)
        c = math.sqrt(((basex2 + x_shift) - x3) ** 2 + (((basey2 + y_shift)) - y3) ** 2)
        A = math.acos((a ** 2 - b ** 2 - c ** 2) / (-2 * b * c))
        # print(A)

        finaldata = []
        direction = direction

        for each in m5:
            each_line = list(each)
            each_line[0] = float(each[0]) + x_shift
            each_line[1] = float(each[1]) + y_shift
            if direction == "S":
                each_line[0], each_line[1] = Srotate(A, each_line[0], each_line[1], x3, y3)
            if direction == "N":
                each_line[0], each_line[1] = Nrotate(A, each_line[0], each_line[1], x3, y3)
            finaldata.append(each_line)

    return finaldata


##
def get_spotfill(fname, fillindex, spotnum):

    xy = open(fname, 'r') 
    point = xy.read()
    xy.close()
    l1 = point.replace('\n', '\t')
    l2 = l1.split('\t')
    l2.pop()
    m1 = np.array(l2[0:len(l2)])
    m2 = m1.reshape(int(len(l2)/10), 10)

    m3 = []
    for each in m2:
        value = [float(k) for k in each[1:-2]]
        if value[fillindex] >= filter_value:
            tmp1 = float(each[-2])
            tmp2 = float(each[-1])
            tmp3 = float(value[fillindex])
            tmp = [tmp1, tmp2, tmp3]
            m3.append(tmp)

    ##

    if spotnum == "spot1":
        finaldata = m3
    else:
        xleft, yleft, xright, yright, direction = Transmessage[spotnum]

        x3 = float(xleft) - resetx 
        y3 = float(yleft) - resety
        basex1 = 0
        basey1 = 0
        x4 = float(xright) - resetx
        y4 = float(yright) - resety
        basex2 = float(imgx)  
        basey2 = float(imgy) 

        x_shift = x3 - basex1
        y_shift = y3 - basey1
        a = math.sqrt((x4 - (basex2 + x_shift)) ** 2 + (y4 - (basey2 + y_shift)) ** 2)
        b = math.sqrt((x4 - x3) ** 2 + (y4 - y3) ** 2)
        c = math.sqrt(((basex2 + x_shift) - x3) ** 2 + (((basey2 + y_shift)) - y3) ** 2)
        A = math.acos((a ** 2 - b ** 2 - c ** 2) / (-2 * b * c))
        # print(A)

        finaldata = []
        direction = direction

        for each in m3:
            each_line = list(each)
            each_line[0] = float(each[0]) + x_shift
            each_line[1] = float(each[1]) + y_shift
            if direction == "S":
                each_line[0], each_line[1] = Srotate(A, each_line[0], each_line[1], x3, y3)
            if direction == "N":
                each_line[0], each_line[1] = Nrotate(A, each_line[0], each_line[1], x3, y3)
            finaldata.append(each_line)

    return finaldata


##
def get_neighbor(spotlist, x0, y0, minv):
    for each in spotlist:
        x = float(each[0])
        y = float(each[1])
        tmp = ((x - x0)**2 + (y - y0)**2)**0.5
        if tmp < minv:
            minv = tmp
            n = spotlist.index(each)
        else:
            continue
    return minv, n


##
def get_basespot(fname1, fname2):

    BaseData = []

    ##fname1

    data = []

    with open(fname1, 'r') as f:
        for line in f.readlines():
            temp = line.split()
            temp[0] = float(temp[0])
            temp[1] = float(temp[1])
            data.append(temp[0:-1])

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
                maxv = int(data[i - 1][1])
                averv = (s + maxv) / 2
                for j in range(len(sess)):
                    sess[j][1] = averv
                all_s.append(sess)
                sess = []
                s = data[i][1]
                sess.append(data[i])
    maxv = int(data[-1][1])
    averv = (s + maxv) / 2
    for j in range(len(sess)):
        sess[j][1] = averv
    all_s.append(sess)

    alls = []
    for i in range(len(all_s)):
        for j in range(len(all_s[i])):
            alls.append(all_s[i][j])

    BaseData.extend(alls)

    ##fname2

    xy = open(fname2, 'r') 
    point = xy.read()
    xy.close()
    l1 = point.replace('\n', '\t')
    l2 = l1.split('\t')
    l2.pop()
    m1 = np.array(l2[0:len(l2)])
    m2 = m1.reshape(int(len(l2) / 10), 10)

    m3 = []
    for each in m2:
        tmp1 = float(each[-2])
        tmp2 = float(each[-1])
        tmp = [tmp1, tmp2]
        m3.append(tmp)

    BaseData.extend(m3)

    return BaseData



localtime_start = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
print("start:"+localtime_start)


############################

##raw
List = []

for i in spot_NUM:
    spotname = "spot" + str(i)
    fname = inputpath_rawxy + spotname + "/" + spotname + ".txt"
    fslide = inputpath_raw + Spot_Slide[spotname] + ".txt"
    spot = get_rawspot(fname, fslide, lab[celltype][0], spotname)
    List.append(spot)

q = 0
for each in List:
    q = q + len(each)
print("select raw spot: " + str(q))

##fill

for i in spot_NUM:
    spotname = "spot" + str(i)
    fname = inputpath_fill + "predict_" + spotname + "_detail.txt"
    spot = get_spotfill(fname, lab[celltype][0], spotname)
    List.append(spot)

q = 0
for each in List:
    q = q + len(each)
print("select raw and fill spot: " + str(q))

##rest

for each in img_spot.keys():
    spotname = img_spot[each]
    fname = inputpath_rest + "predict_" + each + "_detail.txt"
    spot = get_spotfill(fname, lab[celltype][0], spotname)
    List.append(spot)

q = 0
for each in List:
    q = q + len(each)
print("select raw and fill and rest spot: " + str(q))

##basesopt

nol = get_basespot(inputpath_rawxy + "spot1/spot1.txt", inputpath_fill + "predict_spot1_detail.txt")


############验证################

# for i in range(len(List)):
#
#     fig, ax = plt.subplots()
#     ax.set_xlim(left=0, right=8000)
#     ax.set_ylim(bottom=6000, top=0)
#     plt.axis('off')
#
#     x = [k[0] for k in List[i]]
#     y = [k[1] for k in List[i]]
#
#     plt.scatter(x, y, c='r', marker='o', s=3)
#
#     plt.savefig("./verification/"+str(i)+".png")
#
#     plt.clf()


# fig, ax = plt.subplots()
# ax.set_xlim(left=0, right=8000)
# ax.set_ylim(bottom=6000, top=0)
# plt.axis('off')
#
# basespotx = []
# basespoty = []
# for each in nol:
#     basespotx.append(each[0])
#     basespoty.append(each[1])
#
# plt.scatter(basespotx, basespoty, c='r', marker='o', s=3)
#
# plt.savefig("./verification/base.png")
# plt.clf()


############################

SpotL1 = copy.deepcopy(nol)


##
AA = []
BB = []

for i in List:
    aa = []
    bb = []
    for each in i:
        minv, n = get_neighbor(SpotL1, each[0], each[1], 6400)
        aa.append(minv)
        bb.append(n)
        alpha = float(each[2])
        nol[n].append(alpha)
    AA.append(aa)
    BB.append(bb)


############################

newnol = []

for each in nol:
    #Alpha = float(np.mean(each[2:]))
    Alpha = float(sum(each[2:]))
    tmpspot = each[0:2]
    tmpspot.append(Alpha)
    newnol.append(tmpspot)


tilingspot = []
for each in nol:
    tilingspot.append(int(len(each)-2))
print(max(tilingspot), min(tilingspot))


size = []
for each in newnol:
    size.append(float(each[2]))
print(max(size), min(size))
size.sort()

# for each in newnol:
#     each[2] = each[2]/maxsize
#
# size = []
# for each in newnol:
#     size.append(each[2])
# maxsize = max(size)
# print(maxsize)

for i in range(len(AA)):
    if len(AA[i]) != len(List[i]):
        print("Wrong!")


############画布设置###########

fig, ax = plt.subplots()
#fig, ax = plt.subplots(facecolor='gray')

ax.set_xlim(left=0, right=8000)
ax.set_ylim(bottom=6000, top=0) 

plt.axis('off')

x = [k[0] for k in newnol]
y = [k[1] for k in newnol]
alpha = [k[2] for k in newnol]

#plt.scatter(x, y, c=alpha, cmap=lab[celltype][1], marker='o', s=6, linewidths=0, vminv=0, vmax=1)
plt.scatter(x, y, c=alpha, cmap=lab[celltype][1], marker='o', s=6, linewidths=0, vmax=int(size[-100])+1)
#plt.show()



#######################

localtime_end = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
print("end:"+localtime_end)

plt.savefig("./" + celltype + "_2D.pdf")
