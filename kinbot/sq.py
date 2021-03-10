import sqlite3
from sqlite3 import Error
import numpy as np
import io
import logging

"""
create kinbot table with sqlite3
"""


def create_connection(db_file):
    """
    create db connetion for reading/writing
    """
    conn = None
    try:
        conn = sqlite3.connect(db_file)
    except Error as err:
        print(err)

    return conn


def preamble():
    # Converts np.array to TEXT when inserting
    sqlite3.register_adapter(np.ndarray, adapt_array)

    # Converts TEXT to np.array when selecting
    sqlite3.register_converter("array", convert_array)


def adapt_array(arr):
    """ Convert a numpy array to a binary object (BLOB) for sqlite table storage
    : param arr: numpy array
    : return: sqlite blob
    """
    out = io.BytesIO()
    np.save(out, arr, allow_pickle=True)
    out.seek(0)
    return sqlite3.Binary(out.read())


def convert_array(text):

    """
    Convert an sqlite blob to a numpy array
    : param: sqlite blob
    : return: numpy array
    """
    out = io.BytesIO(text)
    out.seek(0)
    return np.load(out, allow_pickle=True)


def create_kinbot_table(conn):

    """
    database elements (
              chemid
              fragment chemids
              l1 energy
              l2 energy
              l3 energy
              zpe
              atoms
              frag xyzs
              symm factor
              frag freqs
              frag red. freqs
              #rotors
              rotor groups
              rotor axis
              rotor symm
              hir potentials
    """

    sql = ''' CREATE TABLE if not exists kinbot (

                 id text PRIMARY KEY,
                 wellorts integer NOT NULL,
                 l1e REAL NOT NULL,
                 l2e REAL,
                 l3e REAL,
                 l1_zpe REAL NOT NULL,
                 l2_zpe REAL,
                 atoms array,
                 l1_xyz array,
                 l2_xyz array,
                 l1_hess array,
                 L2_hess array,
                 L1_freq array,
                 l2_freq array,
                 l1_red_freq array,
                 l2_red_freq array,
                 hir_potentials array
                 ); '''

    try:
        dbc = conn.cursor()
        dbc.execute(sql)
    except Error as err:
        print(err)


def create_kinbot(conn, kinbot_table):

    sql = ''' INSERT OR REPLACE INTO kinbot (id, wellorts, l1e, l2e, l3e, l1_zpe, l2_zpe, atoms, l1_xyz, l2_xyz, l1_hess, l2_hess, l1_freq, l2_freq, l1_red_freq, l2_red_freq, hir_potentials)
              VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?) '''
    cur = conn.cursor()
    cur.execute(sql, kinbot_table)
    conn.commit()

    return cur.lastrowid


def select_all_values(conn, table):

    cur = conn.cursor()
    cur.execute("SELECT * FROM {}".format(table))
    rows = cur.fetchall()
    return rows


def get_sql_data(parent, species, level):
    """
    Collect data from sqlite db for mess file
    sqlite db for kinbot holds information as follows
    col       0: rxn/species name
    col       1: well_prod_ts value
    col 2, 3, 4: L1, L2, L3 energies respectively
    col    5, 6: L1, L2 zpe respectively
    col       7: list of atoms
    col    8, 9: L1, L2 geometries respectively
    col  10, 11: L1, L2 hessian respectively
    col  12, 13: L1, L2 freq respectively
    col  14, 15: L1, L2 reduced freq
    col      16: hir potentials
    """
    db = str(parent) + "/" + str(parent) + "_sql.db"
    conn = create_connection(db)
    data = select_all_values(conn, 'kinbot')
    list_data = []
    for i, row in enumerate(data):
        if row[0] == species:
            for j, col in enumerate(row):
                if type(col) == bytes:
                    array_col = convert_array(col)
                else:
                    array_col = col
                list_data.append(array_col)
    if len(list_data) > 0:
        species = list_data[0]
        if level == 1:
            energy = list_data[2]
            zpe = list_data[5]
            geom = list_data[8]
            hess = list_data[10]
            freq = list_data[12]
            red_freq = list_data[14]
        elif level == 2:
            energy = list_data[3]
            zpe = list_data[6]
            geom = list_data[9]
            hess = list_data[11]
            freq = list_data[13]
            red_freq = list_data[15]
        elif level == 3:
            energy = list_data[4]
            zpe = list_data[6]
            geom = list_data[9]
            hess = list_data[11]
            freq = list_data[13]
            red_freq = list_data[15]
        else:
            logging.error("Error reading SQLite Database")
    else:
        print("{}/{} EMPTY".format(db, row[0]))
        return 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

    hir = list_data[16]
    atoms = list_data[7]

    return list_data, species, energy, zpe, geom, hess, freq, red_freq, hir, atoms


def get_sql_mess_data(data_array, species, level, parameter):
    """
    Pull specific value from sql data array
    levels = L1, L2, L3
    parameter = energy, zpe, atoms, geoms, red_freq, hir
    """
    for row in data_array:
        if species == row[0]:
            if parameter == 'hir':
                return row[16]
            elif parameter == 'red_freq' and level == 'L1':
                return row[14]
            elif parameter == 'red_freq' and level == 'L2':
                return row[15]
            elif parameter == 'geom' and level == 'L1':
                return row[8]
            elif parameter == 'geom' and level == 'L2':
                return row[9]
            elif parameter == 'atoms':
                return row[7]
            elif parameter == 'zpe' and level == 'L1':
                return row[5]
            elif parameter == 'zpe' and level == 'L2':
                return row[6]
            elif parameter == 'energy' and level == 'L1':
                return row[2]
            elif parameter == 'energy' and level == 'L2':
                return row[3]
            elif parameter == 'energy' and level == 'L3':
                return row[4]
        else:
            pass


"""
def main():

    database = 'test.db'
    conn = create_connection(database)
    create_kinbot_table(conn)
    kinbot_table = ('well1', 0, 1.025)
    kinbot = create_kinbot(conn, kinbot_table)

if __name__ == '__main__':
     main()
"""
