from januModule.utils import multiply, parseArray, compute, getParseArray
from januModule.computeCorr import get_chips_generated_array, process_labels_and_corr, process_files, print_files
import sys
import pymysql

def insert_rows_dna_final():

    imageFinal, labelsFinal, corrFinal = process_labels_and_corr()
    chps_generated_labels_list = labelsFinal
    corrList = corrFinal
    chps_generated_array = corrList
    experiment_name = 'replicateSV21'
    list_rows = get_chips_generated_array(experiment_name,chps_generated_array, chps_generated_labels_list)
    cell_tuple_list = []

    query = "INSERT INTO db.cs_correlation(row_num, col_num, cor, experiment_name, row_label, col_label)  VALUES( %s, %s, %s, %s, %s, %s)"

    for item in list_rows:
        print(item)
        cell = (item[0],item[1], item[2], item[3], item[4], item[5])
        cell_tuple_list.append(tuple(cell))

    dbServerName = "127.0.0.1"
    dbUser = "root"
    dbPassword = "Shiva$123"
    dbName = "db"
    charSet = "utf8mb4"
    cusrorType = pymysql.cursors.DictCursor
    connection = pymysql.connect(host=dbServerName, user=dbUser, password=dbPassword,
                                  db=dbName, charset=charSet, cursorclass=cusrorType)

    try:
        cursor = connection.cursor()
        cursor.executemany( query, cell_tuple_list )
        connection.commit()
    except Exception as e:
        print(e)
    finally:
        connection.close()

def get_data_cs():
    con = pymysql.connect(host='localhost',
                          user='root',
                          password='Shiva$123',
                          db='db',
                          charset='utf8mb4',
                          cursorclass=pymysql.cursors.DictCursor)

    ar_rows = list()

    try:
        with con.cursor() as cursor:
            sql = "select `id`, `row_num`,`col_num`, `cor`, `experiment_name`, `row_label`, `col_label`  from db.cs_correlation"
            cursor.execute(sql)
            rows = cursor.fetchall()
            for row in rows:
                t_row = (row['id'], row['row_num'], row['col_num'], row['cor'], row['experiment_name'], row['row_label'], row['col_label'])
                ar_rows.append(t_row)
                print(row['id'])
    except Exception as error:
        print(error)
    finally:
        con.close()

    return ar_rows

def main(args):
    process_files()

if __name__ == '__main__':
    main(sys.argv[1:])

