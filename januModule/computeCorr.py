# note that there is a problem with all 3 bigwig files for four hours treatment of A549 cells with dexamethasone
# this is the link to the bad dataset: https://www.encodeproject.org/experiments/ENCSR748MVX/
import numpy as np
import sys
import os
import deeptools.getScorePerBigWigBin as score_bw
from deeptools.correlation import Correlation

labels_list = []
bw_labels = []
np_array = None
listUnnested = []
BASE_DIR = os.path.dirname(os.path.abspath('../tableRun.py'))
BASE_DIR = os.path.join(BASE_DIR, 'code/januModule')

# useful for quick run with docker compose command
def print_files():
    print('BASE_DIR...')
    print(BASE_DIR)
    a = os.path.join( BASE_DIR, 'treatment-data/files/0h/R1.bigwig' )
    print(a)

    # Generate the local bigWig file paths and cell line labels
    files0 = [ os.path.join(BASE_DIR,'treatment-data/files/0h/R1.bigwig'),
              os.path.join(BASE_DIR,'treatment-data/files/0h/R2.bigwig'),
              os.path.join(BASE_DIR,'treatment-data/files/0h/R3.bigwig')]
    files05 = [os.path.join(BASE_DIR,'treatment-data/files/0.5h/R1.bigwig'),
               os.path.join(BASE_DIR,'treatment-data/files/0.5h/R2.bigwig'),
               os.path.join(BASE_DIR, 'treatment-data/files/0.5h/R3.bigwig')]
    files1 = [os.path.join(BASE_DIR, 'treatment-data/files/1h/R1.bigwig'),
              os.path.join(BASE_DIR, 'treatment-data/files/1h/R2.bigwig'),
              os.path.join(BASE_DIR, 'treatment-data/files/1h/R3.bigwig')]
    files2 = [os.path.join(BASE_DIR, 'treatment-data/files/2h/R1.bigwig'),
              os.path.join(BASE_DIR, 'treatment-data/files/2h/R2.bigwig'),
              os.path.join(BASE_DIR, 'treatment-data/files/2h/R3.bigwig')]
    files3 = [os.path.join(BASE_DIR, 'treatment-data/files/3h/R1.bigwig'),
              os.path.join(BASE_DIR, 'treatment-data/files/3h/R2.bigwig'),
              os.path.join(BASE_DIR, 'treatment-data/files/3h/R3.bigwig')]

    # 4 hours bigwig files are corrupted
    # files4 = [os.path.join(BASE_DIR, 'treatment-data/files/4h/R1.bigwig'),os.path.join(BASE_DIR, 'treatment-data/files/4h/R2.bigwig'),os.path.join(BASE_DIR, 'treatment-data/files/4h/R3.bigwig')]

    files5 = [os.path.join(BASE_DIR, 'treatment-data/files/5h/R1.bigwig'),
              os.path.join(BASE_DIR, 'treatment-data/files/5h/R2.bigwig'),
              os.path.join(BASE_DIR, 'treatment-data/files/5h/R3.bigwig')]
    files6 = [os.path.join(BASE_DIR, 'treatment-data/files/6h/R1.bigwig'),
              os.path.join(BASE_DIR, 'treatment-data/files/6h/R2.bigwig'),
              os.path.join(BASE_DIR, 'treatment-data/files/6h/R3.bigwig')]
    files7 = [os.path.join(BASE_DIR, 'treatment-data/files/7h/R1.bigwig'),
              os.path.join(BASE_DIR, 'treatment-data/files/7h/R2.bigwig'),
              os.path.join(BASE_DIR, 'treatment-data/files/7h/R3.bigwig')]
    files8 = [os.path.join(BASE_DIR, 'treatment-data/files/8h/R1.bigwig'),
              os.path.join(BASE_DIR, 'treatment-data/files/8h/R2.bigwig'),
              os.path.join(BASE_DIR, 'treatment-data/files/8h/R3.bigwig')]
    files10 = [os.path.join(BASE_DIR, 'treatment-data/files/10h/R1.bigwig'),
               os.path.join(BASE_DIR, 'treatment-data/files/10h/R2.bigwig'),
               os.path.join(BASE_DIR, 'treatment-data/files/10h/R3.bigwig')]
    files12 = [os.path.join(BASE_DIR, 'treatment-data/files/12h/R1.bigwig'),
               os.path.join(BASE_DIR, 'treatment-data/files/12h/R2.bigwig'),
               os.path.join(BASE_DIR, 'treatment-data/files/12h/R3.bigwig')]
    files2only = [os.path.join(BASE_DIR, 'treatment-data/files/12h/R1.bigwig'),
                  os.path.join(BASE_DIR, 'treatment-data/files/12h/R2.bigwig')]

def get_labels_and_correlation(
        bw_files,
        # chrs_to_skip,
        bin_size=10000,
        method='pearson',
        fileset_name='result',
        blacklist=None,
        labels=bw_labels,
        output_dir='/Users/baditya02/Downloads/treatment-data/graphs/test/'
):
    my_listUnnested = []
    my_labels_list = []
    assert method in ('pearson', 'spearman'), 'Invalid correlation method'
    # Autogenerate labels from filename if not provided
    if not labels:
        labels = [filename.split( '/' )[-1].split( '.' )[0] for filename in bw_files]
    # Generate a name for the unique combination
    test_name = fileset_name + '_' + method
    if blacklist:
        blacklist_title = 'Blacklisted '
        test_name += '_blacklisted'
    else:
        blacklist_title = ''
    image_name = test_name + '.png'
    # Bin the bigwig data in 10kb increments
    num_reads_per_bin = score_bw.getScorePerBin(
        bw_files,
        bin_size,
        # chrsToSkip=chrs_to_skip,
        blackListFileName=blacklist
    )
    # Write to npz file
    filename = output_dir + test_name + '.npz'
    with open(filename, "wb") as f:
        np.savez_compressed(f, matrix=num_reads_per_bin, labels=labels)
    # Compute the correlations
    corr = Correlation(filename, method, labels=labels)
    np_array = corr.compute_correlation()
    listNested = np_array.tolist()

    def removeNestings(listNest):
        for i in listNest:
            if type(i) == list:
                removeNestings(i)
            else:
                my_listUnnested.append(i)

    removeNestings(listNested)

    plot_title = '{}{} Correlation of {}'.format(
        blacklist_title,
        method.capitalize(),
        fileset_name
    )
    # Create a png file of correlation heatmap
    image_path = output_dir + image_name
    corr.plot_correlation( image_path, plot_title=plot_title )

    # return np_ar
    my_labels_list = labels
    return image_path, my_labels_list, my_listUnnested

def make_corr_plot(
        bw_files,
        # chrs_to_skip,
        bin_size=10000,
        method='pearson',
        fileset_name='result',
        blacklist=None,
        labels=bw_labels,
        output_dir=os.path.join(BASE_DIR, 'treatment-data/graphs/test/')
):

    assert method in ('pearson', 'spearman'), 'Invalid correlation method'
    # Autogenerate labels from filename if not provided
    if not labels:
        labels = [filename.split( '/' )[-1].split( '.' )[0] for filename in bw_files]
    # Generate a name for the unique combination
    test_name = fileset_name + '_' + method
    if blacklist:
        blacklist_title = 'Blacklisted '
        test_name += '_blacklisted'
    else:
        blacklist_title = ''
    image_name = test_name + '.png'
    # Bin the bigwig data in 10kb increments
    num_reads_per_bin = score_bw.getScorePerBin(
        bw_files,
        bin_size,
        # chrsToSkip=chrs_to_skip,
        blackListFileName=blacklist
    )
    # Write to npz file
    filename = output_dir + test_name + '.npz'

    with open( filename, "wb" ) as f:
        np.savez_compressed(f, matrix=num_reads_per_bin, labels=labels)

    # Compute the correlations
    corr = Correlation( filename, method, labels=labels )
    np_array = corr.compute_correlation()

    listNested = np_array.tolist()

    def removeNestings(listNest):
        for i in listNest:
            if type(i)  == list:
                removeNestings(i)
            else:
                listUnnested.append(i)

    removeNestings(listNested)

    plot_title = '{}{} Correlation of {}'.format(
        blacklist_title,
        method.capitalize(),
        fileset_name
    )
    # Create a png file of correlation heatmap
    image_path = output_dir + image_name
    corr.plot_correlation( image_path, plot_title=plot_title )

    # return np_ar
    return image_path

def process_files():
    # Use the Spearman correlation for comparison, it is less susceptible to
    # outliers than Pearson correlation
    method = 'spearman'


    print('BASE_DIR...')
    print(BASE_DIR)
    a = os.path.join(BASE_DIR, 'treatment-data/files/0h/R1.bigwig')

    # Generate the local bigWig file paths and cell line labels
    files0 = [ os.path.join(BASE_DIR,'treatment-data/files/0h/R1.bigwig'),
              os.path.join(BASE_DIR,'treatment-data/files/0h/R2.bigwig'),
              os.path.join(BASE_DIR,'treatment-data/files/0h/R3.bigwig')]
    files05 = [os.path.join(BASE_DIR,'treatment-data/files/0.5h/R1.bigwig'),
               os.path.join(BASE_DIR,'treatment-data/files/0.5h/R2.bigwig'),
               os.path.join(BASE_DIR, 'treatment-data/files/0.5h/R3.bigwig')]
    files1 = [os.path.join(BASE_DIR, 'treatment-data/files/1h/R1.bigwig'),
              os.path.join(BASE_DIR, 'treatment-data/files/1h/R2.bigwig'),
              os.path.join(BASE_DIR, 'treatment-data/files/1h/R3.bigwig')]
    files2 = [os.path.join(BASE_DIR, 'treatment-data/files/2h/R1.bigwig'),
              os.path.join(BASE_DIR, 'treatment-data/files/2h/R2.bigwig'),
              os.path.join(BASE_DIR, 'treatment-data/files/2h/R3.bigwig')]
    files3 = [os.path.join(BASE_DIR, 'treatment-data/files/3h/R1.bigwig'),
              os.path.join(BASE_DIR, 'treatment-data/files/3h/R2.bigwig'),
              os.path.join(BASE_DIR, 'treatment-data/files/3h/R3.bigwig')]

    # 4 hours bigwig files are corrupted
    # files4 = [os.path.join(BASE_DIR, 'treatment-data/files/4h/R1.bigwig'),os.path.join(BASE_DIR, 'treatment-data/files/4h/R2.bigwig'),os.path.join(BASE_DIR, 'treatment-data/files/4h/R3.bigwig')]

    files5 = [os.path.join(BASE_DIR, 'treatment-data/files/5h/R1.bigwig'),
              os.path.join(BASE_DIR, 'treatment-data/files/5h/R2.bigwig'),
              os.path.join(BASE_DIR, 'treatment-data/files/5h/R3.bigwig')]
    files6 = [os.path.join(BASE_DIR, 'treatment-data/files/6h/R1.bigwig'),
              os.path.join(BASE_DIR, 'treatment-data/files/6h/R2.bigwig'),
              os.path.join(BASE_DIR, 'treatment-data/files/6h/R3.bigwig')]
    files7 = [os.path.join(BASE_DIR, 'treatment-data/files/7h/R1.bigwig'),
              os.path.join(BASE_DIR, 'treatment-data/files/7h/R2.bigwig'),
              os.path.join(BASE_DIR, 'treatment-data/files/7h/R3.bigwig')]
    files8 = [os.path.join(BASE_DIR, 'treatment-data/files/8h/R1.bigwig'),
              os.path.join(BASE_DIR, 'treatment-data/files/8h/R2.bigwig'),
              os.path.join(BASE_DIR, 'treatment-data/files/8h/R3.bigwig')]
    files10 = [os.path.join(BASE_DIR, 'treatment-data/files/10h/R1.bigwig'),
               os.path.join(BASE_DIR, 'treatment-data/files/10h/R2.bigwig'),
               os.path.join(BASE_DIR, 'treatment-data/files/10h/R3.bigwig')]
    files12 = [os.path.join(BASE_DIR, 'treatment-data/files/12h/R1.bigwig'),
               os.path.join(BASE_DIR, 'treatment-data/files/12h/R2.bigwig'),
               os.path.join(BASE_DIR, 'treatment-data/files/12h/R3.bigwig')]
    files2only = [os.path.join(BASE_DIR, 'treatment-data/files/12h/R1.bigwig'),
                  os.path.join(BASE_DIR, 'treatment-data/files/12h/R2.bigwig')]
    print(files10)
    # bw_files = files0 + files05 + files1 + files2 + files3 + files5 + files6 + files7 + files8 + files10 + files12
    bw_files = files0

    labels0 = ['R1 - 0 h', 'R2 - 0 h', 'R3 - 0 h']
    labels05 = ['R1 - 0.5 h', 'R2 - 0.5 h', 'R3 - 0.5 h']
    labels1 = ['R1 - 1 h', 'R2 - 1 h', 'R3 - 1 h']
    labels2 = ['R1 - 2 h', 'R2 - 2 h', 'R3 - 2 h']
    labels3 = ['R1 - 3 h', 'R2 - 3 h', 'R3 - 3 h']

    # 4 hour bigwig files are corrupted
    # labels4 = ['R1 - 4 h', 'R2 - 4 h', 'R3 - 4 h']

    labels5 = ['R1 - 5 h', 'R2 - 5 h', 'R3 - 5 h']
    labels6 = ['R1 - 6 h', 'R2 - 6 h', 'R3 - 6 h']
    labels7 = ['R1 - 7 h', 'R2 - 7 h', 'R3 - 7 h']
    labels8 = ['R1 - 8 h', 'R2 - 8 h', 'R3 - 8 h']
    labels10 = ['R1 - 10 h', 'R2 - 10 h', 'R3 - 10 h']
    labels12 = ['R1 - 12 h', 'R2 - 12 h', 'R3 - 12 h']
    labels2only = ['R1 - 12 h', 'R2 - 12 h']

    # bw_labels = labels0 + labels05 + labels1 + labels2 + labels3 + labels5 + labels6 + labels7 + labels8 + labels10 + labels12
    bw_labels = labels0

    imagename = make_corr_plot(
        bw_files,
        # chrs_to_skip,
        method=method,
        fileset_name='A549 cells treated with 100 nM dexamethasone',
        labels=bw_labels,
        output_dir=os.path.join(BASE_DIR, 'treatment-data/graphs/test')
    )

def process_labels_and_corr():
    # Use the Spearman correlation for comparison, it is less susceptible to
    # outliers than Pearson correlation
    method = 'spearman'

    # Generate the local bigWig file paths and cell line labels
    files0 = ['/Users/baditya02/Documents/sites/svpython/pyconnectWdl/januModule/treatment-data/files/0h/R1.bigwig',
              '/Users/baditya02/Documents/sites/svpython/pyconnectWdl/januModule/treatment-data/files/0h/R2.bigwig',
              '/Users/baditya02/Documents/sites/svpython/pyconnectWdl/januModule/treatment-data/files/0h/R3.bigwig']
    files05 = ['/Users/baditya02/Documents/sites/svpython/pyconnectWdl/januModule/treatment-data/files/0.5h/R1.bigwig',
               '/Users/baditya02/Documents/sites/svpython/pyconnectWdl/januModule/treatment-data/files/0.5h/R2.bigwig',
               '/Users/baditya02/Documents/sites/svpython/pyconnectWdl/januModule/treatment-data/files/0.5h/R3.bigwig']
    files1 = ['/Users/baditya02/Documents/sites/svpython/pyconnectWdl/januModule/treatment-data/files/1h/R1.bigwig',
              '/Users/baditya02/Documents/sites/svpython/pyconnectWdl/januModule/treatment-data/files/1h/R2.bigwig',
              '/Users/baditya02/Documents/sites/svpython/pyconnectWdl/januModule/treatment-data/files/1h/R3.bigwig']
    files2 = ['/Users/baditya02/Documents/sites/svpython/pyconnectWdl/januModule/treatment-data/files/2h/R1.bigwig',
              '/Users/baditya02/Documents/sites/svpython/pyconnectWdl/januModule/treatment-data/files/2h/R2.bigwig',
              '/Users/baditya02/Documents/sites/svpython/pyconnectWdl/januModule/treatment-data/files/2h/R3.bigwig']
    files3 = ['/Users/baditya02/Documents/sites/svpython/pyconnectWdl/januModule/treatment-data/files/3h/R1.bigwig',
              '/Users/baditya02/Documents/sites/svpython/pyconnectWdl/januModule/treatment-data/files/3h/R2.bigwig',
              '/Users/baditya02/Documents/sites/svpython/pyconnectWdl/januModule/treatment-data/files/3h/R3.bigwig']

    # 4 hours bigwig files are corrupted
    # files4 = ['/Users/baditya02/Documents/sites/svpython/pyconnectWdl/januModule/treatment-data/files/4h/R1.bigwig','/Users/baditya02/Documents/sites/svpython/pyconnectWdl/januModule/treatment-data/files/4h/R2.bigwig','/Users/baditya02/Documents/sites/svpython/pyconnectWdl/januModule/treatment-data/files/4h/R3.bigwig']

    files5 = ['/Users/baditya02/Documents/sites/svpython/pyconnectWdl/januModule/treatment-data/files/5h/R1.bigwig',
              '/Users/baditya02/Documents/sites/svpython/pyconnectWdl/januModule/treatment-data/files/5h/R2.bigwig',
              '/Users/baditya02/Documents/sites/svpython/pyconnectWdl/januModule/treatment-data/files/5h/R3.bigwig']
    files6 = ['/Users/baditya02/Documents/sites/svpython/pyconnectWdl/januModule/treatment-data/files/6h/R1.bigwig',
              '/Users/baditya02/Documents/sites/svpython/pyconnectWdl/januModule/treatment-data/files/6h/R2.bigwig',
              '/Users/baditya02/Documents/sites/svpython/pyconnectWdl/januModule/treatment-data/files/6h/R3.bigwig']
    files7 = ['/Users/baditya02/Documents/sites/svpython/pyconnectWdl/januModule/treatment-data/files/7h/R1.bigwig',
              '/Users/baditya02/Documents/sites/svpython/pyconnectWdl/januModule/treatment-data/files/7h/R2.bigwig',
              '/Users/baditya02/Documents/sites/svpython/pyconnectWdl/januModule/treatment-data/files/7h/R3.bigwig']
    files8 = ['/Users/baditya02/Documents/sites/svpython/pyconnectWdl/januModule/treatment-data/files/8h/R1.bigwig',
              '/Users/baditya02/Documents/sites/svpython/pyconnectWdl/januModule/treatment-data/files/8h/R2.bigwig',
              '/Users/baditya02/Documents/sites/svpython/pyconnectWdl/januModule/treatment-data/files/8h/R3.bigwig']
    files10 = ['/Users/baditya02/Documents/sites/svpython/pyconnectWdl/januModule/treatment-data/files/10h/R1.bigwig',
               '/Users/baditya02/Documents/sites/svpython/pyconnectWdl/januModule/treatment-data/files/10h/R2.bigwig',
               '/Users/baditya02/Documents/sites/svpython/pyconnectWdl/januModule/treatment-data/files/10h/R3.bigwig']
    files12 = ['/Users/baditya02/Documents/sites/svpython/pyconnectWdl/januModule/treatment-data/files/12h/R1.bigwig',
               '/Users/baditya02/Documents/sites/svpython/pyconnectWdl/januModule/treatment-data/files/12h/R2.bigwig',
               '/Users/baditya02/Documents/sites/svpython/pyconnectWdl/januModule/treatment-data/files/12h/R3.bigwig']
    files2only = ['/Users/baditya02/Documents/sites/svpython/pyconnectWdl/januModule/treatment-data/files/12h/R1.bigwig',
                  '/Users/baditya02/Documents/sites/svpython/pyconnectWdl/januModule/treatment-data/files/12h/R2.bigwig']

    # bw_files = files0 + files05 + files1 + files2 + files3 + files5 + files6 + files7 + files8 + files10 + files12
    bw_files = files0

    labels0 = ['R1 - 0 h', 'R2 - 0 h', 'R3 - 0 h']
    labels05 = ['R1 - 0.5 h', 'R2 - 0.5 h', 'R3 - 0.5 h']
    labels1 = ['R1 - 1 h', 'R2 - 1 h', 'R3 - 1 h']
    labels2 = ['R1 - 2 h', 'R2 - 2 h', 'R3 - 2 h']
    labels3 = ['R1 - 3 h', 'R2 - 3 h', 'R3 - 3 h']

    # 4 hour bigwig files are corrupted
    # labels4 = ['R1 - 4 h', 'R2 - 4 h', 'R3 - 4 h']

    labels5 = ['R1 - 5 h', 'R2 - 5 h', 'R3 - 5 h']
    labels6 = ['R1 - 6 h', 'R2 - 6 h', 'R3 - 6 h']
    labels7 = ['R1 - 7 h', 'R2 - 7 h', 'R3 - 7 h']
    labels8 = ['R1 - 8 h', 'R2 - 8 h', 'R3 - 8 h']
    labels10 = ['R1 - 10 h', 'R2 - 10 h', 'R3 - 10 h']
    labels12 = ['R1 - 12 h', 'R2 - 12 h', 'R3 - 12 h']
    labels2only = ['R1 - 12 h', 'R2 - 12 h']

    #bw_labels = labels0 + labels05 + labels1 + labels2 + labels3 + labels5 + labels6 + labels7 + labels8 + labels10 + labels12
    bw_labels = labels0

    imagename, labels_list_final, listUnnested_final = get_labels_and_correlation(
        bw_files,
        # chrs_to_skip,
        method=method,
        fileset_name='A549 cells treated with 100 nM dexamethasone',
        labels=bw_labels,
        output_dir='/Users/baditya02/Documents/sites/svpython/pyconnectWdl/januModule/treatment-data/graphs/test/'
    )

    return imagename, labels_list_final, listUnnested_final


def get_chips_generated_array(experimentName = 'Replicate5',  corrList = [1, 2, 3, 4, 5, 6, 7, 8, 9], chps_generated_labels_list = ["R0", "R1", "R2"]):

    length_labels_list = len(chps_generated_labels_list)
    chps_generated_array = []
    rowsNumber = length_labels_list
    colsNumber = length_labels_list
    labelRow = chps_generated_labels_list
    labelCol = chps_generated_labels_list

    for i in range(len(corrList)):
        for j in range(rowsNumber):
            x = j + 1
            for k in range(colsNumber):
                y = k + 1
                corrInterList = [x, y, corrList[i], experimentName, labelRow[j], labelCol[k]]
                if i != len( corrList ) - 1:
                    i = i + 1
                if x == rowsNumber and y == colsNumber:
                    break
            else:
                continue
            break
        else:
            continue
        break

    return chps_generated_array

def main(args):
    process_labels_and_corr()

if __name__ == '__main__':
    main(sys.argv[1:] )

