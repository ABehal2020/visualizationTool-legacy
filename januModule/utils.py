import numpy as np

def getParseArray(experumentName, chipse_array, labels_list):

    array = np.array(chipse_array)
    rows = array.shape[0]
    cols = array.shape[1]
    rowTable = []
    colTable = []
    corrTable = []
    rowLabelTable = []
    colLabelTable = []

    for i in range(0, len(labels_list)):
        rowLabelTable.append(labels_list[i])
        colLabelTable.append(labels_list[i])


    for i in range( 0, rows):
        for j in range( 0, cols):
            if j == 0:
                rowTable.append(array[i, j])
            if j == 1:
                colTable.append(array[i, j])

            if j == 2:
                corrTable.append(array[i, j])

    return rowTable, colTable, corrTable, rowLabelTable, colLabelTable


